"""
03_make_sample_lists.py

Downloads gnomAD HGDP+1KGP sample metadata and writes per-population
sample list files to $BASE/data/.

Also writes a list of SGDP Melanesian samples (Papuan populations)
based on the SGDP metadata table.

Outputs (all in BASE/data/):
  samples_african.txt          — 1KGP + HGDP continental African, admixed excluded
  samples_melanesian_hgdp.txt  — HGDP Papuan + Bougainville
  samples_melanesian_sgdp.txt  — SGDP Papuan populations
  samples_eastasian.txt        — 1KGP + HGDP East Asian
  samples_southasian.txt       — 1KGP + HGDP South Asian
  samples_all_keep.txt         — union of all above (for VCF filtering)
  samples_exclude.txt          — ASW + ACB (admixed African-descent, excluded)

Usage:
    python analysis/cluster/03_make_sample_lists.py \
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
"""

import argparse
import os
import pandas as pd
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--base",
                    default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints",
                    help="Cluster working directory")
args = parser.parse_args()

OUT_DIR = os.path.join(args.base, "data")
os.makedirs(OUT_DIR, exist_ok=True)

# ─── gnomAD HGDP+1KGP sample metadata ────────────────────────────────────────
META_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
            "release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/"
            "gnomad_meta_updated.tsv")

print("Downloading gnomAD HGDP+1KGP sample metadata...")
meta_path = os.path.join(OUT_DIR, "gnomad_hgdp_1kg_metadata.tsv")
if not os.path.exists(meta_path):
    r = requests.get(META_URL, timeout=60)
    r.raise_for_status()
    with open(meta_path, "wb") as f:
        f.write(r.content)
    print(f"  Saved → {meta_path}")
else:
    print(f"  Already exists → {meta_path}")

meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
print(f"  {len(meta)} samples loaded")

# Column names may vary — inspect available columns
print("\nAvailable columns (first 20):", list(meta.columns[:20]))

# Typical columns: 's' (sample ID), 'population' or 'pop', 'project_pop',
# 'hgdp_tgp_meta.Population', 'hgdp_tgp_meta.Genetic region', etc.
# Adjust the column names below if they differ in the actual metadata file.

# Try to find the population column
pop_col = None
for candidate in ["population", "pop", "project_pop",
                   "hgdp_tgp_meta.Population", "Population"]:
    if candidate in meta.columns:
        pop_col = candidate
        break

sample_col = "s" if "s" in meta.columns else meta.columns[0]

if pop_col is None:
    print("\nERROR: Could not detect population column. "
          "Inspect gnomad_hgdp_1kg_metadata.tsv and set pop_col manually.")
    raise SystemExit(1)

print(f"\nUsing sample column: '{sample_col}', population column: '{pop_col}'")
print("Unique populations:", sorted(meta[pop_col].dropna().unique()))

# ─── Population assignments ───────────────────────────────────────────────────

# 1KGP population codes → role
AFRICAN_1KGP     = {"YRI", "LWK", "ESN", "GWD", "MSL"}          # keep
ADMIXED_EXCLUDE  = {"ASW", "ACB"}                                  # exclude (admixed)
EAST_ASIAN_1KGP  = {"CHB", "JPT", "CHS", "CDX", "KHV"}
SOUTH_ASIAN_1KGP = {"GIH", "PJL", "BEB", "STU", "ITU"}
EUROPEAN_1KGP    = {"CEU", "TSI", "FIN", "GBR", "IBS"}

# HGDP population labels (check against your metadata)
AFRICAN_HGDP     = {"Yoruba", "Mandenka"}
MELANESIAN_HGDP  = {"Papuan", "PapuanHighlands", "PapuanSepik", "Bougainville"}
EAST_ASIAN_HGDP  = {"Han", "Japanese", "She", "Tujia", "Naxi", "Lahu", "Dai",
                     "Cambodians", "Yakut"}
SOUTH_ASIAN_HGDP = {"Balochi", "Brahui", "Makrani", "Pathan", "Sindhi",
                     "Burusho", "Hazara", "Kalash"}
EUROPEAN_HGDP    = {"French", "Sardinian", "Basque", "Orcadian", "BergamoItalian",
                     "Tuscan", "Russian", "Adygei"}

pops = meta[pop_col].fillna("")

african     = meta[pops.isin(AFRICAN_1KGP | AFRICAN_HGDP)][sample_col]
melanesian  = meta[pops.isin(MELANESIAN_HGDP)][sample_col]
east_asian  = meta[pops.isin(EAST_ASIAN_1KGP | EAST_ASIAN_HGDP)][sample_col]
south_asian = meta[pops.isin(SOUTH_ASIAN_1KGP | SOUTH_ASIAN_HGDP)][sample_col]
european    = meta[pops.isin(EUROPEAN_1KGP | EUROPEAN_HGDP)][sample_col]
excluded    = meta[pops.isin(ADMIXED_EXCLUDE)][sample_col]

print(f"\nSample counts:")
print(f"  African (1KGP + HGDP, excl. admixed): {len(african)}")
print(f"  Melanesian (HGDP only):               {len(melanesian)}")
print(f"  East Asian (1KGP + HGDP):             {len(east_asian)}")
print(f"  South Asian (1KGP + HGDP):            {len(south_asian)}")
print(f"  European (1KGP + HGDP):               {len(european)}")
print(f"  Excluded (ASW + ACB):                  {len(excluded)}")

def write_list(samples, path):
    samples.reset_index(drop=True).to_csv(path, index=False, header=False)
    print(f"  → {path}  ({len(samples)} samples)")

print("\nWriting sample lists...")
write_list(african,     os.path.join(OUT_DIR, "samples_african.txt"))
write_list(melanesian,  os.path.join(OUT_DIR, "samples_melanesian_hgdp.txt"))
write_list(east_asian,  os.path.join(OUT_DIR, "samples_eastasian.txt"))
write_list(south_asian, os.path.join(OUT_DIR, "samples_southasian.txt"))
write_list(european,    os.path.join(OUT_DIR, "samples_european.txt"))
write_list(excluded,    os.path.join(OUT_DIR, "samples_exclude.txt"))

# ─── SGDP Melanesian samples ──────────────────────────────────────────────────
# Sample IDs confirmed directly from Simons.vcf.gz (data/sgdp/Simons.vcf.gz).
# 14 × S_Papuan, 1 × B_Papuan-15, 2 × S_Bougainville = 17 total.
# NOTE: must be defined BEFORE all_keep so SGDP samples are included.

SGDP_MELANESIAN_SAMPLES = [
    "S_Papuan-1",  "S_Papuan-2",  "S_Papuan-3",  "S_Papuan-4",
    "S_Papuan-5",  "S_Papuan-6",  "S_Papuan-7",  "S_Papuan-8",
    "S_Papuan-9",  "S_Papuan-10", "S_Papuan-11", "S_Papuan-12",
    "S_Papuan-13", "S_Papuan-14",
    "B_Papuan-15",
    "S_Bougainville-1", "S_Bougainville-2",
]

sgdp_melanesian = pd.Series(SGDP_MELANESIAN_SAMPLES)
print(f"\n  SGDP Melanesian (Papuan + Bougainville): {len(sgdp_melanesian)} samples")
write_list(sgdp_melanesian, os.path.join(OUT_DIR, "samples_melanesian_sgdp.txt"))

# all_keep includes SGDP so 05_merge_filter_vcfs.sh doesn't silently drop them
all_keep = pd.concat([african, melanesian, sgdp_melanesian,
                      east_asian, south_asian, european]).drop_duplicates()
write_list(all_keep,    os.path.join(OUT_DIR, "samples_all_keep.txt"))

print("\nDone. Copy data/ directory to the cluster before running 02_extract_vcfs.sh.")
