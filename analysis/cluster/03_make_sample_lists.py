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
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints
"""

import argparse
import os
import pandas as pd
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--base",
                    default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints",
                    help="Cluster working directory")
args = parser.parse_args()

OUT_DIR = os.path.join(args.base, "data")
os.makedirs(OUT_DIR, exist_ok=True)

# ─── gnomAD HGDP+1KGP sample metadata ────────────────────────────────────────
META_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/"
            "release/3.1/secondary_analyses/hgdp_1kg/"
            "gnomad.hgdp_1kg.sample_qc_metadata.tsv")

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

# HGDP population labels (check against your metadata)
AFRICAN_HGDP     = {"Yoruba", "Mandenka"}
MELANESIAN_HGDP  = {"Papuan", "Bougainville"}
EAST_ASIAN_HGDP  = {"Han", "Japanese", "She", "Tujia", "Naxi", "Lahu", "Dai",
                     "Cambodians", "Yakut"}
SOUTH_ASIAN_HGDP = {"Balochi", "Brahui", "Makrani", "Pathan", "Sindhi",
                     "Burusho", "Hazara", "Kalash"}

pops = meta[pop_col].fillna("")

african     = meta[pops.isin(AFRICAN_1KGP | AFRICAN_HGDP)][sample_col]
melanesian  = meta[pops.isin(MELANESIAN_HGDP)][sample_col]
east_asian  = meta[pops.isin(EAST_ASIAN_1KGP | EAST_ASIAN_HGDP)][sample_col]
south_asian = meta[pops.isin(SOUTH_ASIAN_1KGP | SOUTH_ASIAN_HGDP)][sample_col]
excluded    = meta[pops.isin(ADMIXED_EXCLUDE)][sample_col]

print(f"\nSample counts:")
print(f"  African (1KGP + HGDP, excl. admixed): {len(african)}")
print(f"  Melanesian (HGDP only):               {len(melanesian)}")
print(f"  East Asian (1KGP + HGDP):             {len(east_asian)}")
print(f"  South Asian (1KGP + HGDP):            {len(south_asian)}")
print(f"  Excluded (ASW + ACB):                  {len(excluded)}")

def write_list(samples, path):
    samples.reset_index(drop=True).to_csv(path, index=False, header=False)
    print(f"  → {path}  ({len(samples)} samples)")

print("\nWriting sample lists...")
write_list(african,     os.path.join(OUT_DIR, "samples_african.txt"))
write_list(melanesian,  os.path.join(OUT_DIR, "samples_melanesian_hgdp.txt"))
write_list(east_asian,  os.path.join(OUT_DIR, "samples_eastasian.txt"))
write_list(south_asian, os.path.join(OUT_DIR, "samples_southasian.txt"))
write_list(excluded,    os.path.join(OUT_DIR, "samples_exclude.txt"))

all_keep = pd.concat([african, melanesian, east_asian, south_asian]).drop_duplicates()
write_list(all_keep,    os.path.join(OUT_DIR, "samples_all_keep.txt"))

# ─── SGDP Melanesian samples ──────────────────────────────────────────────────
# SGDP metadata from: https://reichdata.hms.harvard.edu/pub/datasets/sgdp/
# Papuan populations in SGDP: "Papuan" entries (various subgroups)
# The SGDP VCF uses sample IDs of the form: S_<population>-<N>
# Manually curated list of Papuan-group populations in SGDP:

SGDP_PAPUAN_POPS = [
    "Papuan",
    "New_Guinea",
]

SGDP_META_URL = ("https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
                 "SGDP_metadata.279public.21signedLetter.samples.txt")

sgdp_meta_path = os.path.join(OUT_DIR, "sgdp_metadata.txt")
if not os.path.exists(sgdp_meta_path):
    try:
        print("\nDownloading SGDP metadata...")
        r = requests.get(SGDP_META_URL, timeout=30)
        r.raise_for_status()
        with open(sgdp_meta_path, "wb") as f:
            f.write(r.content)
        print(f"  Saved → {sgdp_meta_path}")
    except Exception as e:
        print(f"  WARNING: Could not download SGDP metadata: {e}")
        print("  Download manually from: https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/")
        sgdp_meta_path = None

if sgdp_meta_path and os.path.exists(sgdp_meta_path):
    sgdp = pd.read_csv(sgdp_meta_path, sep="\t", low_memory=False)
    print("SGDP columns:", list(sgdp.columns[:10]))
    # Column may be 'SGDP-Population ID', 'Population', etc.
    # Filter for Papuan-group populations
    for col in ["Population ID", "SGDP-Population ID", "Population", "population"]:
        if col in sgdp.columns:
            papuan_mask = sgdp[col].str.contains(
                "|".join(SGDP_PAPUAN_POPS), case=False, na=False)
            sgdp_sample_col = [c for c in ["Sample ID", "SGDP-ID", "sample_id"]
                                if c in sgdp.columns][0]
            sgdp_melanesian = sgdp.loc[papuan_mask, sgdp_sample_col]
            print(f"  SGDP Melanesian (Papuan): {len(sgdp_melanesian)} samples")
            write_list(sgdp_melanesian,
                       os.path.join(OUT_DIR, "samples_melanesian_sgdp.txt"))
            break

print("\nDone. Copy data/ directory to the cluster before running 02_extract_vcfs.sh.")
