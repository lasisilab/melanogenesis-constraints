"""
08_compute_pbs.py

Computes per-gene PBS (Population Branch Statistic) for the 4 approved scans
using the Hudson FST estimator (Bhatia et al. 2013, Genome Research).

Hudson FST (ratio-of-means, not mean-of-ratios):
  For populations X and Y at site i with frequencies p_X, p_Y and
  chromosome counts n_X, n_Y:

    num_i   = (p_X - p_Y)^2
              - p_X(1-p_X)/(n_X - 1)
              - p_Y(1-p_Y)/(n_Y - 1)

    denom_i = p_X(1-p_Y) + p_Y(1-p_X)

  FST_XY = Σ num_i / Σ denom_i   (sum over all sites with data in both pops)

PBS transform to branch length:
  T_XY = -ln(1 - FST_XY)         (FST clipped to [0, 0.9999] before log)

PBS for target population X, outgroup Y, distant outgroup Z:
  PBS_X = (T_XY + T_XZ - T_YZ) / 2

Approved scans (Tina Lasisi):
  PBS-1: Target=African,    Outgroup=SouthAsian, Distant=Papuan/Melanesian
  PBS-2: Target=African,    Outgroup=European,   Distant=Papuan/Melanesian
  PBS-3: Target=Papuan,     Outgroup=SouthAsian, Distant=African
  PBS-4: Target=Papuan,     Outgroup=European,   Distant=African

Reference:
  Bhatia G et al. (2013) Estimating and interpreting FST: The impact of rare
  variants. Genome Research 23(9):1514–1521.

Output: output/pbs_per_gene.csv
  gene, chrom, start, end, n_sites_shared,
  fst_afr_sas, fst_afr_eur, fst_afr_mel,
  fst_sas_mel, fst_eur_mel, fst_sas_eur,
  pbs1_african, pbs2_african,
  pbs3_melanesian, pbs4_melanesian
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

BASE    = args.base
BED     = os.path.join(BASE, "data", "gene_regions.bed")
VCF_DIR = os.path.join(BASE, "vcf", "final")
OUT_DIR = os.path.join(BASE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

POPULATIONS = ["african", "melanesian", "eastasian", "southasian", "european"]

# PBS scan definitions: (target, outgroup1, distant_outgroup)
PBS_SCANS = [
    ("pbs1_african",    "african",    "southasian", "melanesian"),
    ("pbs2_african",    "african",    "european",   "melanesian"),
    ("pbs3_melanesian", "melanesian", "southasian", "african"),
    ("pbs4_melanesian", "melanesian", "european",   "african"),
]

# All FST pairs needed (derive from scans)
FST_PAIRS = set()
for _, a, b, c in PBS_SCANS:
    FST_PAIRS.add(tuple(sorted([a, b])))
    FST_PAIRS.add(tuple(sorted([a, c])))
    FST_PAIRS.add(tuple(sorted([b, c])))
FST_PAIRS = sorted(FST_PAIRS)

# ── Load gene regions ──────────────────────────────────────────────────────
bed = pd.read_csv(BED, sep="\t", header=None,
                  names=["chrom", "start", "end", "gene"])
print(f"Loaded {len(bed)} gene regions from {BED}")


def query_ac_an(vcf_path: str, region: str) -> pd.DataFrame:
    """
    Returns DataFrame [pos, ac, an] for all biallelic SNPs in region.
    Uses bcftools +fill-tags to recompute AC/AN from subsetted genotypes.
    """
    cmd = (
        f"bcftools view -r {region} {vcf_path} "
        f"| bcftools +fill-tags -- -t AC,AN "
        f"| bcftools query -f '%POS\\t%INFO/AC\\t%INFO/AN\\n'"
    )
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        return pd.DataFrame(columns=["pos", "ac", "an"])

    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        try:
            pos = int(parts[0])
            ac  = int(parts[1].split(",")[0])
            an  = int(parts[2])
            if an > 0:
                rows.append((pos, ac, an))
        except ValueError:
            continue

    return pd.DataFrame(rows, columns=["pos", "ac", "an"])


def hudson_fst(merged: pd.DataFrame, pop_a: str, pop_b: str) -> float:
    """
    Computes Hudson FST between two populations using the ratio-of-means
    estimator (Bhatia et al. 2013).

    merged must have columns: ac_{pop_a}, an_{pop_a}, ac_{pop_b}, an_{pop_b}
    Returns NaN if fewer than 2 shared sites.
    """
    a_ac = merged[f"ac_{pop_a}"].values.astype(float)
    a_an = merged[f"an_{pop_a}"].values.astype(float)
    b_ac = merged[f"ac_{pop_b}"].values.astype(float)
    b_an = merged[f"an_{pop_b}"].values.astype(float)

    # Require data in both populations
    valid = (a_an > 1) & (b_an > 1)
    if valid.sum() < 2:
        return np.nan

    p1 = a_ac[valid] / a_an[valid]
    p2 = b_ac[valid] / b_an[valid]
    n1 = a_an[valid]
    n2 = b_an[valid]

    num   = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    denom = p1 * (1 - p2) + p2 * (1 - p1)

    sum_denom = denom.sum()
    if sum_denom <= 0:
        return np.nan

    return float(num.sum() / sum_denom)


def fst_to_t(fst: float) -> float:
    """Branch-length transform: T = -ln(1 - FST). FST clipped to [0, 0.9999]."""
    if np.isnan(fst):
        return np.nan
    fst_clipped = max(0.0, min(fst, 0.9999))
    return float(-np.log(1.0 - fst_clipped))


def compute_pbs(t_xy: float, t_xz: float, t_yz: float) -> float:
    """PBS_X = (T_XY + T_XZ - T_YZ) / 2."""
    if any(np.isnan(v) for v in [t_xy, t_xz, t_yz]):
        return np.nan
    return (t_xy + t_xz - t_yz) / 2.0


# ── Main loop ──────────────────────────────────────────────────────────────
records = []
n_genes = len(bed)

# Determine which populations are actually needed (skip unused ones)
pops_needed = set()
for _, a, b, c in PBS_SCANS:
    pops_needed.update([a, b, c])

for idx, row in bed.iterrows():
    gene       = row["gene"]
    chrom      = row["chrom"]
    start      = int(row["start"])
    end        = int(row["end"])
    region_str = f"{chrom}:{start}-{end}"

    # Query AC/AN for each needed population
    pop_data: dict[str, pd.DataFrame] = {}
    for pop in pops_needed:
        vcf = os.path.join(VCF_DIR, f"{pop}.vcf.gz")
        if not os.path.exists(vcf):
            print(f"  WARNING: {vcf} not found")
            pop_data[pop] = pd.DataFrame(columns=["pos", "ac", "an"])
        else:
            pop_data[pop] = query_ac_an(vcf, region_str)

    # Merge all populations on POS (inner join — only sites with data in all)
    # Use outer merge so each FST pair can use its own shared sites
    all_pos_frames = {}
    for pop, df in pop_data.items():
        if len(df) > 0:
            all_pos_frames[pop] = df.rename(columns={
                "ac": f"ac_{pop}", "an": f"an_{pop}"
            }).set_index("pos")

    if not all_pos_frames:
        rec = {"gene": gene, "chrom": chrom, "start": start, "end": end,
               "n_sites_shared": 0}
        for _, *_ in PBS_SCANS:
            pass
        records.append(rec)
        continue

    # Full outer merge across all populations
    merged = None
    for pop, frame in all_pos_frames.items():
        if merged is None:
            merged = frame
        else:
            merged = merged.join(frame, how="outer")
    merged = merged.fillna(0).reset_index()

    n_shared = len(merged)

    rec = {
        "gene":           gene,
        "chrom":          chrom,
        "start":          start,
        "end":            end,
        "n_sites_shared": n_shared,
    }

    # Compute Hudson FST for each needed pair
    fst_cache: dict[tuple, float] = {}
    for pop_a, pop_b in FST_PAIRS:
        if f"ac_{pop_a}" not in merged.columns or f"ac_{pop_b}" not in merged.columns:
            fst_cache[(pop_a, pop_b)] = np.nan
        else:
            fst_cache[(pop_a, pop_b)] = hudson_fst(merged, pop_a, pop_b)
        key = f"fst_{pop_a}_{pop_b}"
        rec[key] = fst_cache[(pop_a, pop_b)]

    # Compute PBS for each scan
    for scan_name, target, outgroup, distant in PBS_SCANS:
        pair_to  = tuple(sorted([target,   outgroup]))
        pair_td  = tuple(sorted([target,   distant]))
        pair_od  = tuple(sorted([outgroup, distant]))

        fst_to = fst_cache.get(pair_to, np.nan)
        fst_td = fst_cache.get(pair_td, np.nan)
        fst_od = fst_cache.get(pair_od, np.nan)

        t_to = fst_to_t(fst_to)
        t_td = fst_to_t(fst_td)
        t_od = fst_to_t(fst_od)

        pbs_val = compute_pbs(t_to, t_td, t_od)
        # Floor at 0: negative PBS means no excess population-specific
        # differentiation — biologically equivalent to zero signal.
        if not np.isnan(pbs_val):
            pbs_val = max(0.0, pbs_val)
        rec[scan_name] = pbs_val

    records.append(rec)

    if (idx + 1) % 10 == 0 or (idx + 1) == n_genes:
        print(f"  [{idx+1}/{n_genes}] {gene}  n_sites={n_shared}  "
              f"PBS1_afr={rec.get('pbs1_african', float('nan')):.4f}  "
              f"PBS3_mel={rec.get('pbs3_melanesian', float('nan')):.4f}")

# ── Write output ───────────────────────────────────────────────────────────
# Reorder columns: identifiers → FST pairs → PBS scans
fst_cols = [f"fst_{a}_{b}" for a, b in FST_PAIRS]
pbs_cols = [name for name, *_ in PBS_SCANS]
id_cols  = ["gene", "chrom", "start", "end", "n_sites_shared"]

out_df  = pd.DataFrame(records)
col_order = id_cols + [c for c in fst_cols if c in out_df.columns] + pbs_cols
out_df  = out_df[[c for c in col_order if c in out_df.columns]]

out_csv = os.path.join(OUT_DIR, "pbs_per_gene.csv")
out_df.to_csv(out_csv, index=False)
print(f"\nSaved → {out_csv}  ({len(out_df)} genes)")

# ── Quick summary ──────────────────────────────────────────────────────────
print("\nTop 10 genes by PBS-1 (African, SouthAsian outgroup, Melanesian distant):")
if "pbs1_african" in out_df.columns:
    top = out_df[["gene", "pbs1_african", "pbs2_african",
                  "pbs3_melanesian", "pbs4_melanesian"]].dropna(
        subset=["pbs1_african"]).nlargest(10, "pbs1_african")
    print(top.to_string(index=False))

print("\nTop 10 genes by PBS-3 (Melanesian, SouthAsian outgroup, African distant):")
if "pbs3_melanesian" in out_df.columns:
    top = out_df[["gene", "pbs3_melanesian", "pbs4_melanesian",
                  "pbs1_african", "pbs2_african"]].dropna(
        subset=["pbs3_melanesian"]).nlargest(10, "pbs3_melanesian")
    print(top.to_string(index=False))

print("\nDone!")
