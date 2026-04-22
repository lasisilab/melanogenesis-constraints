"""
11_compute_pbs_genomewide.py

Computes per-gene PBS for all protein-coding genes on one chromosome,
using the same Hudson FST estimator as 08_compute_pbs.py.

Called by 11_compute_pbs_genomewide.sh (SLURM array, one job per chromosome).

Inputs (all in GW_DIR = vcf/genomewide/):
    {pop}.chr{CHR}.vcf.gz  — per-population whole-chromosome VCFs from step 10
    data/gencode_v38_protein_coding_10kb.bed  — gene windows from step 09

Output:
    output/pbs_genomewide_chr{CHR}.csv
        gene_name, gene_id, chrom, start, end, n_snps,
        fst_afr_sas, fst_afr_mel, fst_mel_sas,
        fst_afr_eur, fst_eur_mel,
        pbs1_african, pbs2_african, pbs3_melanesian, pbs4_melanesian

PBS scan design (mirrors 08_compute_pbs.py):
    PBS-1: African target,    South Asian outgroup, Melanesian distant
    PBS-2: African target,    European outgroup,    Melanesian distant
    PBS-3: Melanesian target, South Asian outgroup, African distant
    PBS-4: Melanesian target, European outgroup,    African distant

Usage:
    python 11_compute_pbs_genomewide.py --chr 1 \
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
"""

import argparse
import os
import subprocess
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--chr",  required=True, help="Chromosome number (1–22)")
    p.add_argument("--base",
                   default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints")
    return p.parse_args()


# ── Hudson FST (ratio-of-means) ────────────────────────────────────────────────
def hudson_fst(ac1, an1, ac2, an2):
    """
    Vectorised Hudson FST (Bhatia et al. 2013, ratio-of-means estimator).
    ac/an are 1-D arrays across SNPs.  Returns a scalar FST for the window.
    Negative FST floored to 0.
    """
    # Drop sites monomorphic in both pops or missing in either
    mask = (an1 > 0) & (an2 > 0)
    ac1, an1, ac2, an2 = ac1[mask], an1[mask], ac2[mask], an2[mask]
    if len(ac1) == 0:
        return np.nan

    p1 = ac1 / an1
    p2 = ac2 / an2

    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (an1 - 1) - p2 * (1 - p2) / (an2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)

    denom_sum = den.sum()
    if denom_sum == 0:
        return np.nan

    fst = num.sum() / denom_sum
    return max(0.0, fst)


def fst_to_t(fst):
    """T = -ln(1 - FST), clipped to avoid log(0)."""
    fst_clipped = min(fst, 0.9999)
    return -np.log(1.0 - fst_clipped)


def pbs(t_target_out, t_target_dist, t_out_dist):
    """PBS formula.  Negative values floored to 0."""
    val = (t_target_out + t_target_dist - t_out_dist) / 2.0
    return max(0.0, val)


# ── Query AC/AN for a region from a VCF ───────────────────────────────────────
def query_ac_an(vcf_path, region):
    """
    Returns (ac_array, an_array) for all biallelic SNPs in region.
    Uses bcftools +fill-tags to recompute AC/AN from subsetted genotypes.
    """
    cmd = (
        f"bcftools view -r {region} {vcf_path} "
        f"| bcftools +fill-tags -- -t AC,AN "
        f"| bcftools query -f '%INFO/AC\\t%INFO/AN\\n'"
    )
    try:
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=120
        )
        if result.returncode != 0 or not result.stdout.strip():
            return np.array([]), np.array([])
        rows = [line.split("\t") for line in result.stdout.strip().split("\n")]
        ac = np.array([int(r[0]) for r in rows if len(r) == 2])
        an = np.array([int(r[1]) for r in rows if len(r) == 2])
        return ac, an
    except (subprocess.TimeoutExpired, Exception):
        return np.array([]), np.array([])


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    chrom_num = args.chr
    chrom = f"chr{chrom_num}"

    gw_dir  = os.path.join(args.base, "vcf", "genomewide")
    out_dir = os.path.join(args.base, "output")
    data_dir = os.path.join(args.base, "data")
    os.makedirs(out_dir, exist_ok=True)

    out_csv = os.path.join(out_dir, f"pbs_genomewide_chr{chrom_num}.csv")
    if os.path.exists(out_csv):
        print(f"Output already exists — {out_csv}.  Delete to rerun.")
        return

    # ── VCF paths ─────────────────────────────────────────────────────────────
    vcfs = {
        pop: os.path.join(gw_dir, f"{pop}.chr{chrom_num}.vcf.gz")
        for pop in ("african", "melanesian", "southasian", "european")
    }
    for pop, path in vcfs.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing VCF: {path}")

    # ── Load gene annotation for this chromosome ───────────────────────────────
    bed_path = os.path.join(data_dir, "gencode_v38_protein_coding_10kb.bed")
    if not os.path.exists(bed_path):
        raise FileNotFoundError(
            f"{bed_path} not found.  Run 09_fetch_gene_annotation.py first."
        )

    genes = pd.read_csv(bed_path, sep="\t", comment="#",
                        names=["chrom", "start", "end",
                               "gene_name", "gene_id", "strand"])
    genes = genes[genes["chrom"] == chrom].reset_index(drop=True)
    print(f"chr{chrom_num}: {len(genes)} genes to process")

    if len(genes) == 0:
        print(f"No genes found for {chrom} — exiting.")
        return

    # ── Compute PBS per gene ───────────────────────────────────────────────────
    results = []
    for i, row in genes.iterrows():
        region = f"{row['chrom']}:{row['start']}-{row['end']}"

        ac_afr, an_afr = query_ac_an(vcfs["african"],    region)
        ac_mel, an_mel = query_ac_an(vcfs["melanesian"],  region)
        ac_sas, an_sas = query_ac_an(vcfs["southasian"],  region)
        ac_eur, an_eur = query_ac_an(vcfs["european"],    region)

        n_snps = len(ac_afr)

        # Pairwise FST (only defined when both pops have data for the same sites)
        # We use the per-population AC/AN independently (ratio-of-means is
        # computed per pair, not requiring identical site sets — monomorphic
        # sites are dropped within each pair by hudson_fst()).
        def pair_fst(ac_a, an_a, ac_b, an_b):
            # Align to shorter array if lengths differ (can happen at edges)
            n = min(len(ac_a), len(ac_b))
            if n == 0:
                return np.nan
            return hudson_fst(ac_a[:n], an_a[:n], ac_b[:n], an_b[:n])

        f_afr_sas = pair_fst(ac_afr, an_afr, ac_sas, an_sas)
        f_afr_mel = pair_fst(ac_afr, an_afr, ac_mel, an_mel)
        f_mel_sas = pair_fst(ac_mel, an_mel, ac_sas, an_sas)
        f_afr_eur = pair_fst(ac_afr, an_afr, ac_eur, an_eur)
        f_eur_mel = pair_fst(ac_eur, an_eur, ac_mel, an_mel)

        def safe_pbs(f_to, f_td, f_od):
            if any(np.isnan(x) for x in (f_to, f_td, f_od)):
                return np.nan
            return pbs(fst_to_t(f_to), fst_to_t(f_td), fst_to_t(f_od))

        pbs1 = safe_pbs(f_afr_sas, f_afr_mel, f_mel_sas)  # AFR/SAS/MEL
        pbs2 = safe_pbs(f_afr_eur, f_afr_mel, f_eur_mel)  # AFR/EUR/MEL
        pbs3 = safe_pbs(f_mel_sas, f_afr_mel, f_afr_sas)  # MEL/SAS/AFR
        pbs4 = safe_pbs(f_eur_mel, f_afr_mel, f_afr_eur)  # MEL/EUR/AFR

        results.append({
            "gene_name":       row["gene_name"],
            "gene_id":         row["gene_id"],
            "chrom":           row["chrom"],
            "start":           row["start"],
            "end":             row["end"],
            "n_snps":          n_snps,
            "fst_afr_sas":     f_afr_sas,
            "fst_afr_mel":     f_afr_mel,
            "fst_mel_sas":     f_mel_sas,
            "fst_afr_eur":     f_afr_eur,
            "fst_eur_mel":     f_eur_mel,
            "pbs1_african":    pbs1,
            "pbs2_african":    pbs2,
            "pbs3_melanesian": pbs3,
            "pbs4_melanesian": pbs4,
        })

        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(genes)} genes processed...")

    out_df = pd.DataFrame(results)
    out_df.to_csv(out_csv, index=False)
    print(f"\nSaved → {out_csv}  ({len(out_df)} genes)")


if __name__ == "__main__":
    main()
