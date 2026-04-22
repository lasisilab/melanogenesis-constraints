"""
07_compute_pi.py

Computes per-gene nucleotide diversity (π) for all 5 populations.

π_gene = (1 / L) × Σ_i  2 · p_i · (1 − p_i)

where p_i = ALT allele frequency at biallelic SNP i (computed from AC/AN
            after bcftools +fill-tags, so it reflects the subsetted samples),
      L   = region length in base pairs (gene body ± 10 kb from BED file).

Only sites with AN > 0 contribute. Sites with missing genotypes reduce the
effective sample size at that position but do not affect π by default because
we normalise by L (region length), not by n_sites. This matches the standard
convention used in population genetics studies comparing diversity across loci
of different sizes.

Output: output/pi_per_gene.csv
  gene, chrom, start, end, region_bp,
  n_sites_african, pi_african,
  n_sites_melanesian, pi_melanesian,
  n_sites_eastasian, pi_eastasian,
  n_sites_southasian, pi_southasian,
  n_sites_european, pi_european
"""

import argparse
import os
import subprocess
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--base",
                    default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints",
                    help="Cluster base directory")
args = parser.parse_args()

BASE       = args.base
BED        = os.path.join(BASE, "data", "gene_regions.bed")
VCF_DIR    = os.path.join(BASE, "vcf", "final")
OUT_DIR    = os.path.join(BASE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

POPULATIONS = ["african", "melanesian", "eastasian", "southasian", "european"]

# ── Load gene regions ──────────────────────────────────────────────────────
bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene"])
print(f"Loaded {len(bed)} gene regions from {BED}")


def query_ac_an(vcf_path: str, region: str) -> pd.DataFrame:
    """
    Returns a DataFrame with columns [pos, ac, an] for all biallelic SNPs
    in the given region.  Uses bcftools +fill-tags to recompute AC/AN from
    the subsetted genotypes (the INFO fields in the final VCF still reflect
    the full cohort).
    """
    cmd = (
        f"bcftools view -r {region} {vcf_path} "
        f"| bcftools +fill-tags -- -t AC,AN "
        f"| bcftools query -f '%POS\\t%INFO/AC\\t%INFO/AN\\n'"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        # Empty region or no variants — not an error
        return pd.DataFrame(columns=["pos", "ac", "an"])

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        try:
            # AC can be comma-separated for multiallelics; we filtered to
            # biallelic in step 05 so this should always be a single int.
            pos = int(parts[0])
            ac  = int(parts[1].split(",")[0])
            an  = int(parts[2])
            if an > 0:
                rows.append((pos, ac, an))
        except ValueError:
            continue

    return pd.DataFrame(rows, columns=["pos", "ac", "an"])


def compute_pi(ac: pd.Series, an: pd.Series, region_bp: int) -> float:
    """
    π per base pair = (1/L) × Σ 2·p·(1−p)
    Returns NaN if no callable sites.
    """
    valid = an > 0
    if valid.sum() == 0 or region_bp == 0:
        return np.nan
    p = ac[valid] / an[valid]
    return float((2 * p * (1 - p)).sum() / region_bp)


# ── Main loop ──────────────────────────────────────────────────────────────
records = []
n_genes = len(bed)

for idx, row in bed.iterrows():
    gene       = row["gene"]
    chrom      = row["chrom"]
    start      = int(row["start"])
    end        = int(row["end"])
    region_bp  = end - start
    region_str = f"{chrom}:{start}-{end}"

    rec = {
        "gene":      gene,
        "chrom":     chrom,
        "start":     start,
        "end":       end,
        "region_bp": region_bp,
    }

    for pop in POPULATIONS:
        vcf = os.path.join(VCF_DIR, f"{pop}.vcf.gz")
        if not os.path.exists(vcf):
            print(f"  WARNING: {vcf} not found — skipping {pop}")
            rec[f"n_sites_{pop}"] = np.nan
            rec[f"pi_{pop}"]      = np.nan
            continue

        df_sites = query_ac_an(vcf, region_str)
        n_sites  = len(df_sites)
        pi_val   = compute_pi(df_sites["ac"], df_sites["an"], region_bp)

        rec[f"n_sites_{pop}"] = n_sites
        rec[f"pi_{pop}"]      = pi_val

    records.append(rec)

    if (idx + 1) % 10 == 0 or (idx + 1) == n_genes:
        print(f"  [{idx+1}/{n_genes}] {gene}  "
              f"π_afr={rec.get('pi_african', float('nan')):.5f}  "
              f"π_mel={rec.get('pi_melanesian', float('nan')):.5f}")

# ── Write output ───────────────────────────────────────────────────────────
out_df = pd.DataFrame(records)
out_csv = os.path.join(OUT_DIR, "pi_per_gene.csv")
out_df.to_csv(out_csv, index=False)
print(f"\nSaved → {out_csv}  ({len(out_df)} genes)")

# ── Quick summary ──────────────────────────────────────────────────────────
print("\nMedian π by population:")
for pop in POPULATIONS:
    col = f"pi_{pop}"
    if col in out_df:
        med = out_df[col].median()
        n   = out_df[col].notna().sum()
        print(f"  {pop:12s}: median π = {med:.6f}  (n={n} genes)")

print("\nDone!")
