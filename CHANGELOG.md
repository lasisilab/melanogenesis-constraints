# Changelog

All notable changes to the melanogenesis-constraints project are documented here.

---

## [Unreleased] — Phase 2 Population Genomics Pipeline

### Fixed
- `analysis/cluster/03_make_sample_lists.py` — `all_keep` was built before `sgdp_melanesian` was defined, causing `--force-samples` in 05 to silently drop all 17 SGDP Papuan/Bougainville samples from `melanesian.vcf.gz`; fixed by moving `sgdp_melanesian` definition above `all_keep` and including it in the union; verified: `melanesian.vcf.gz` now has exactly 47 samples (30 HGDP + 17 SGDP)
- Phase 2 VCF prep complete: all 5 population VCFs verified in `vcf/final/` (african 184,859 vars/747 samples; southasian 167,302/790; eastasian 137,663/718; european 135,432/788; melanesian 47 samples confirmed)

### Added (Quarto website updates)
- `phase2.qmd` — Phase 2 population genomics page: π stats + tables, PBS by category, top-10 gene tables per scan, MC1R case study with auto-generated interpretation, PBS tree figure, data provenance table
- `_quarto.yml` — Phase 2 added to navbar
- `index.qmd` — updated with analysis phase summary table, key findings section, corrected gene count and scope

### Added (Phase 2 — PBS population branch trees)
- `analysis/phase2_pbs_trees.py` — draws unrooted 3-leaf PBS star trees for top-3 genes per scan + MC1R; 4×4 grid (rows = scans, cols = top-1/2/3/MC1R); branch lengths = PBS values for all 3 populations derived from pairwise FST; target arm colored by functional category; output: `output/figure_phase2_pbs_trees.png/pdf`

### Added (Phase 2 — π and PBS visualisation)
- `analysis/phase2_pi_pbs_analysis.py` — visualises π and PBS outputs:
  - `figure_phase2_pi.png/pdf`: 2×3 layout — π by category for all 5 populations (A–E) + MC1R bar chart (F); MC1R starred in each panel
  - `figure_phase2_pbs.png/pdf`: 2×3 layout — all 4 PBS scans as category boxplots (A–D) + lollipop of top genes across all scans (E)
  - `output/table_pi_pigment_specific.csv` — raw π for all pigment-specific genes across 5 populations
  - `output/table_pi_top_divergence.csv` — top 10 genes by max−min π divergence across populations
  - Statistics: Kruskal-Wallis + Dunn's posthoc on π by category, Spearman ρ(π, LOEUF), MC1R percentile ranks; negative PBS floored to 0

### Added (Phase 2 — π and PBS analysis scripts)
- `analysis/cluster/07_compute_pi.sh` — SLURM single job; calls 07_compute_pi.py
- `analysis/cluster/07_compute_pi.py` — per-gene nucleotide diversity π = (1/L)Σ 2p(1−p) for all 5 populations; outputs `output/pi_per_gene.csv`
- `analysis/cluster/08_compute_pbs.sh` — SLURM single job; calls 08_compute_pbs.py; can run concurrently with 07
- `analysis/cluster/08_compute_pbs.py` — per-gene PBS using Hudson FST estimator (Bhatia et al. 2013) for all 4 approved scans (PBS-1/2 African target, PBS-3/4 Melanesian target); outputs `output/pbs_per_gene.csv`

### Added (Phase 0 — Network Position & KEGG Connectivity)
- `analysis/fetch_kegg_pathways.py` — fetches number of KEGG pathways per gene via KEGG REST API (2 bulk requests); saves to `data/kegg_pathway_counts.csv`
- `analysis/phase0_network_constraint.py` — replicates Tina's original figure layout (1×3): Panel A betweenness scatter, Panel B KEGG pathway scatter, Panel C LOEUF-by-category boxplot with white-circle medians; degrades to 1×2 (A+C) if KEGG data not yet fetched
- `phase0.qmd` — Quarto page for Phase 0; displays figure, summary tables, per-category statistics, and auto-generated interpretation text
- Phase 0 added to `_quarto.yml` nav

### Added
- `analysis/cluster/00_setup_dirs.sh` — creates working directory tree on cluster
- `analysis/cluster/01_make_regions_bed.py` — builds hg38 BED file of 128 gene regions ±10 kb
- `analysis/cluster/02_extract_vcfs.sh` — SLURM array: extracts gene regions from gnomAD HGDP+1KGP v3.1.2 (HTTP stream) and local SGDP Simons.vcf.gz; skip-if-exists for gnomAD
- `analysis/cluster/03_make_sample_lists.py` — generates per-population sample list files from gnomAD metadata; SGDP Melanesian sample IDs hardcoded from confirmed VCF header
- `analysis/cluster/04_liftover_sgdp.sh` — SLURM array: renames SGDP chromosomes (1→chr1), lifts hg19→hg38 with bcftools +liftover, sorts, normalizes
- `analysis/cluster/05_merge_filter_vcfs.sh` — SLURM array: merges gnomAD + lifted SGDP, filters to biallelic SNPs, writes per-population subsets
- `analysis/cluster/06_concat_vcfs.sh` — concatenates per-chromosome VCFs into one file per population (african, melanesian, eastasian, southasian, european)
- `CLUSTER_REQUEST.md` — compute and storage resource documentation for UMich Great Lakes

### Pipeline configuration
- Cluster: UMich Great Lakes (SLURM), account `tlasisi0`
- Modules: `Bioinformatics bcftools/1.21 htslib`
- 1 CPU per job, `%6` concurrency, email notifications to `ypryor@umich.edu`
- Base directory: `/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/`

### Populations
| Population | N | Source |
|---|---|---|
| African | 747 | gnomAD HGDP+1KGP (YRI, LWK, ESN, GWD, MSL, Yoruba, Mandenka) |
| Melanesian | 47 | gnomAD HGDP (PapuanHighlands, PapuanSepik, Bougainville) + SGDP (Papuan) |
| East Asian | 718 | gnomAD HGDP+1KGP |
| South Asian | 790 | gnomAD HGDP+1KGP |
| European | 788 | gnomAD HGDP+1KGP (CEU, TSI, FIN, GBR, IBS + HGDP European) |

### PBS scan design (approved by T. Lasisi)
| Scan | Target | Outgroup | Distant outgroup |
|---|---|---|---|
| PBS-1 | African | South Asian | Papuan |
| PBS-2 | African | European | Papuan |
| PBS-3 | Papuan | South Asian | African |
| PBS-4 | Papuan | European | African |

### Added MC1R case study subsection to phase1_ancestry.qmd
- Documents the paradox: highest LOEUF (1.967) + 53/54 tissue breadth + African purifying selection signal invisible to LoF metrics
- Explains why LOEUF misses African purifying selection (missense, not LoF)
- Planned follow-up: Missense Z + PBS to capture population-specific constraint
- Marked as "in progress" for PEQG poster

### Decided: iHS out of scope for PEQG 2026
iHS requires phased whole-chromosome data and genome-wide normalization — incompatible with targeted extraction. Deferred to future work. PEQG deliverables are π and PBS only.

---

## [2026-04-07] — Phase 1 Analyses Complete

### Added
- `analysis/phase1_gnomad_v4_ancestry.py` — ancestry-stratified LOEUF using Han et al. 2025 (Nat Commun) Supplementary Data 11 (MOESM11); NFE, AFR, EAS, SAS, ASJ
- `analysis/phase1_gtex_analysis.py` — GTEx v8 tissue expression breadth (54 tissues, TPM > 1)
- `analysis/phase1_phylop_analysis.py` — PhyloP 100-way vertebrate conservation scores
- `phase1_ancestry.qmd` — Quarto page for ancestry-stratified LOEUF analysis
- `data/network_constraint_ancestry_loeuf.csv` — master dataset with v4 ancestry LOEUF columns
- `data/network_constraint_gtex.csv` — master dataset with LOEUF + tissue breadth
- `data/network_constraint_phylop.csv` — master dataset with LOEUF + PhyloP
- `output/figure_phase1_ancestry_loeuf.png/pdf`
- `output/figure_phase1_gtex.png/pdf`
- `output/figure_phase1_phylop.png/pdf`

### Key findings
- Pigment-specific genes are significantly more LoF-tolerant than generic signaling genes (KW p = 1.70e-08)
- Constraint architecture holds in non-European ancestry data (NFE/AFR LOEUF concordant)
- MC1R: highest LOEUF in dataset (1.967), expressed in 53/54 tissues; European red-hair alleles may inflate NFE estimate
- Tissue breadth positively correlated with LOEUF (Spearman ρ > 0)

---

## [2026-04-01] — Initial Setup

### Added
- 129-gene set from Raghunath et al. (2015) melanocyte signaling network
- LOEUF from gnomAD v2.1.1, annotated by functional category (6 categories) and disease class
- Network betweenness centrality from Raghunath et al.
- `LOEUF_by_functional_category.xlsx` — source data
- Quarto website (`index.qmd`, `phase1.qmd`, `_quarto.yml`)
