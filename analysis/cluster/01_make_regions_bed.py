"""
01_make_regions_bed.py

Builds a sorted, merged BED file of gene regions (gene body ± FLANK bp)
from data/phylop_scores.csv, which already has hg38 coordinates.

Output: data/gene_regions.bed (chrom, start, end, gene)
        also uploaded/copied to the cluster data/ directory.

Usage:
    python analysis/cluster/01_make_regions_bed.py
    python analysis/cluster/01_make_regions_bed.py --flank 10000   # default
"""

import argparse
import os
import pandas as pd

FLANK = 10_000  # bp flanking either side of gene body

parser = argparse.ArgumentParser()
parser.add_argument("--flank", type=int, default=FLANK,
                    help="Flanking region in bp (default: 10000)")
args = parser.parse_args()

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PHYLOP_CSV  = os.path.join(PROJECT_DIR, "data", "phylop_scores.csv")
OUT_BED     = os.path.join(PROJECT_DIR, "data", "gene_regions.bed")

df = pd.read_csv(PHYLOP_CSV)
df = df.dropna(subset=["chrom", "start", "end"])

# Ensure chrom has "chr" prefix (bcftools expects it for hg38)
df["chrom"] = df["chrom"].astype(str)
df.loc[~df["chrom"].str.startswith("chr"), "chrom"] = "chr" + df["chrom"]

# Add flanking
df["bed_start"] = (df["start"] - args.flank).clip(lower=0).astype(int)
df["bed_end"]   = (df["end"]   + args.flank).astype(int)

bed = df[["chrom", "bed_start", "bed_end", "gene"]].copy()
bed.columns = ["chrom", "start", "end", "gene"]

# Sort by chrom then start (numeric chrom order)
chrom_order = {f"chr{i}": i for i in range(1, 23)}
chrom_order.update({"chrX": 23, "chrY": 24, "chrM": 25})
bed["sort_key"] = bed["chrom"].map(chrom_order).fillna(99).astype(int)
bed = bed.sort_values(["sort_key", "start"]).drop(columns="sort_key")

bed.to_csv(OUT_BED, sep="\t", index=False, header=False)
print(f"Wrote {len(bed)} regions to {OUT_BED}")
print(f"Flanking: ±{args.flank:,} bp")
print(f"Chromosomes covered: {sorted(bed['chrom'].unique(), key=lambda c: chrom_order.get(c, 99))}")
