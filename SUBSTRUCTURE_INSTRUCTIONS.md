# Population Substructure Fix — Instructions for Claude Code

## Context

This is the `melanogenesis-constraints` project. It analyzes evolutionary constraint and
population-specific adaptation across the 129-gene Raghunath et al. (2015) melanocyte
signaling network.

The current PBS pipeline pools all African samples into one branch and all South Asian
samples into another. A reviewer flagged this as insensitive to population substructure —
particularly the distinction between West African and East African populations, which have
different evolutionary histories and selection pressures. East African populations are also
the primary target for the planned hypoxia-constraints extension of this project.

**Goal:** Split the African PBS branch into West African and East African sub-branches,
and split the Melanesian branch into Highland vs. Lowland, without breaking any existing
analyses. Existing PBS columns in outputs should be preserved.

---

## Pipeline Overview

The cluster pipeline lives in `analysis/cluster/`. Scripts run on the University of Michigan
Great Lakes HPC cluster. Key scripts for this task:

| Script | Role |
|---|---|
| `03_make_sample_lists.py` | Pulls gnomAD HGDP+1KGP metadata, writes per-population sample `.txt` files to `data/` |
| `13_make_sgdp_only_sample_lists.py` | Writes SGDP sample lists to `data/sgdp_only/` |
| `08_compute_pbs.py` | Reads sample lists + merged VCFs, computes per-gene PBS for approved scan configurations, writes `output/pbs_per_gene.csv` |
| `11_compute_pbs_genomewide.py` / `15_compute_pbs_sgdp_genomewide.py` | Genome-wide PBS (same logic, per-chromosome) |

Existing approved PBS scans (defined in `08_compute_pbs.py`):
- PBS-1: Target=African, Outgroup=SouthAsian, Distant=Melanesian
- PBS-2: Target=African, Outgroup=European, Distant=Melanesian
- PBS-3: Target=Melanesian, Outgroup=SouthAsian, Distant=African
- PBS-4: Target=Melanesian, Outgroup=European, Distant=African

---

## Changes Required

### 1. `03_make_sample_lists.py` — Split African and Melanesian sample lists

**Current African populations (pooled):**
```python
AFRICAN_1KGP  = {"YRI", "LWK", "ESN", "GWD", "MSL"}
AFRICAN_HGDP  = {"Yoruba", "Mandenka"}
```

**Replace with sub-group definitions:**
```python
# West African — Nigeria, Sierra Leone, Gambia
AFRICAN_WEST_1KGP  = {"YRI", "ESN", "GWD", "MSL"}
AFRICAN_WEST_HGDP  = {"Yoruba", "Mandenka"}

# East African — Kenya (Luhya is the only East African population in HGDP+1KGP)
# NOTE: LWK (Luhya in Webuye, Kenya) are Bantu-speaking — not altitude-adapted.
# True East African signal (Somali, Masai, Dinka) requires SGDP (see script 13 changes).
AFRICAN_EAST_1KGP  = {"LWK"}
AFRICAN_EAST_HGDP  = set()   # no East African HGDP populations in current pull
```

**Keep the existing combined `samples_african.txt` for backward compatibility.**
Add two new output files:
```python
african_west = meta[pops.isin(AFRICAN_WEST_1KGP | AFRICAN_WEST_HGDP)][sample_col]
african_east = meta[pops.isin(AFRICAN_EAST_1KGP | AFRICAN_EAST_HGDP)][sample_col]

write_list(african_west, os.path.join(OUT_DIR, "samples_african_west.txt"))
write_list(african_east, os.path.join(OUT_DIR, "samples_african_east.txt"))
```

**Melanesian split — already in the HGDP labels:**
```python
MELANESIAN_HIGHLANDS_HGDP = {"PapuanHighlands"}
MELANESIAN_LOWLANDS_HGDP  = {"PapuanSepik", "Papuan", "Bougainville"}
# Keep existing MELANESIAN_HGDP combined list for backward compatibility
```

Add two new output files:
```python
melanesian_highlands = meta[pops.isin(MELANESIAN_HIGHLANDS_HGDP)][sample_col]
melanesian_lowlands  = meta[pops.isin(MELANESIAN_LOWLANDS_HGDP)][sample_col]

write_list(melanesian_highlands, os.path.join(OUT_DIR, "samples_melanesian_highlands.txt"))
write_list(melanesian_lowlands,  os.path.join(OUT_DIR, "samples_melanesian_lowlands.txt"))
```

---

### 2. `13_make_sgdp_only_sample_lists.py` — Split SGDP African into East/West

The SGDP African set (in the existing `AFRICAN` dict) includes East African populations
not present in HGDP+1KGP: **Somali, Masai, Dinka, Luo, BantuKenya**. These are the
populations most relevant for altitude adaptation and the hypoxia follow-up.

**Replace the single `AFRICAN` set with sub-group sets:**
```python
AFRICAN_WEST = {
    "Yoruba", "Mandenka", "Biaka", "Esan", "Gambian", "Mende",
    "Mbuti", "Ju_hoan_North", "Khomani_San",
}
AFRICAN_EAST = {
    "Somali", "Masai", "Dinka", "Luo", "BantuKenya",
}
AFRICAN_SOUTH = {
    "BantuSouthAfrica", "BantuTutsi", "BantuHerero",
}
# Keep combined AFRICAN set for backward compatibility
AFRICAN = AFRICAN_WEST | AFRICAN_EAST | AFRICAN_SOUTH
```

Update `POP_MAP` to map each sub-group:
```python
for pop in AFRICAN_WEST:  POP_MAP[pop] = "african_west"
for pop in AFRICAN_EAST:  POP_MAP[pop] = "african_east"
for pop in AFRICAN_SOUTH: POP_MAP[pop] = "african_south"
# Keep legacy "african" mapping for backward compat — map combined set
```

Write additional output files in `data/sgdp_only/`:
- `samples_sgdp_african_west.txt`
- `samples_sgdp_african_east.txt`
- `samples_sgdp_african_south.txt` (small n, probably combine with east or exclude)

**Important:** Check n for each SGDP sub-group before proceeding — if
`african_east` has fewer than ~10 samples, PBS estimates will be unreliable.
Print counts and flag any group with n < 15.

---

### 3. `08_compute_pbs.py` — Add new PBS scan configurations

**Do not modify existing PBS-1 through PBS-4.** Add new scans using the new sample list
files. Append new columns to `output/pbs_per_gene.csv`.

New scans to add:
```
PBS-5: Target=AfricanWest,  Outgroup=SouthAsian, Distant=Melanesian
PBS-6: Target=AfricanEast,  Outgroup=SouthAsian, Distant=Melanesian
       (Use SGDP East African samples: Somali, Masai, Dinka, Luo, BantuKenya)
PBS-7: Target=MelanesianHighlands, Outgroup=SouthAsian, Distant=AfricanWest
PBS-8: Target=MelanesianLowlands,  Outgroup=SouthAsian, Distant=AfricanWest
```

New output columns (append to existing CSV, do not drop existing columns):
```
pbs5_african_west
pbs6_african_east
pbs7_mel_highlands
pbs8_mel_lowlands
```

The SGDP East African samples live in a separate VCF from the gnomAD HGDP+1KGP samples.
PBS-6 will require merging SGDP East African samples with the existing gnomAD Melanesian
and SouthAsian samples for FST computation — see how `08_compute_pbs.py` currently handles
the SGDP Melanesian merge for reference.

---

### 4. Downstream: `analysis/phase2_pbs_within_network.py` and `phase2_pbs_within_network_combined.py`

Once `output/pbs_per_gene.csv` has the new PBS columns, re-run the within-network
regression analysis for the new branches. The regression pattern is identical to the
existing AFR/MEL branches — Spearman ρ + adjusted linear model + permutation test
for each predictor (τ, betweenness, KEGG count, LOEUF, tissue breadth).

Add branches `AFR_WEST`, `AFR_EAST`, `MEL_HIGHLANDS`, `MEL_LOWLANDS` to whatever
branch loop currently iterates over `["AFR", "MEL", "EUR", "SAS"]`.

Update `output/phase2_pbs_within_network_combined.txt` to include results for new branches.

---

## What NOT to touch

- `output/pbs_per_gene.csv` existing columns — append only, do not drop
- `output/phase2_pbs_within_network_combined.txt` — add new branch sections, preserve existing
- Any figure-generation scripts — leave those for a separate pass once results are confirmed
- The genome-wide PBS scripts (11, 15) — defer substructure changes there until the per-gene
  results look sensible

---

## Practical notes

- Scripts run on Great Lakes HPC. Local paths use `/Users/ypryor/GitHub/melanogenesis-constraints/`;
  cluster paths use `/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/`.
- The merged VCFs in `vcf/final/` should already contain all HGDP+1KGP samples — no
  re-extraction needed for steps 1–3. The SGDP East African samples for PBS-6 may require
  re-running `05_merge_filter_vcfs.sh` if SGDP non-Melanesian samples were excluded earlier.
  Check whether SGDP samples beyond the 17 Melanesian ones are present in the merged VCF
  before assuming a re-merge is needed.
- `push_to_greatlakes.sh` handles syncing code to the cluster.

---

## Priority order

1. Modify `03_make_sample_lists.py` — add West/East African and Highland/Lowland Melanesian
   sample list outputs (backward-compatible)
2. Modify `13_make_sgdp_only_sample_lists.py` — split SGDP African into East/West/South,
   check n per group
3. Modify `08_compute_pbs.py` — add PBS-5 through PBS-8, append new columns
4. Rerun within-network regression for new branches
5. Decide whether SGDP East African re-merge is needed (check VCF sample headers first)
