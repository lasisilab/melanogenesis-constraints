# Melanogenesis Constraints

Analysis pipeline for the PEQG 2026 poster:
**"Testing Network-Based Predictions of Genetic Constraint in Melanogenesis."**

We test whether a gene's position in the melanogenesis signaling network predicts
its evolutionary constraint, combining LoF intolerance (LOEUF), phylogenetic
conservation (PhyloP), tissue expression breadth (GTEx), and population genomics
(nucleotide diversity π, PBS) across 5 populations.

**Gene set:** Raghunath et al. 2015 — 129-gene melanocyte signaling network,
classified into 6 functional categories.

**Live site:** built and deployed by GitHub Actions on every push to `main`.
See `.github/workflows/publish.yml`.

---

## Repository layout

```
.
├── README.md                    This file
├── _quarto.yml                  Quarto site config (navbar, theme, render list)
├── styles.css                   Site styles
├── index.qmd                    Site landing page (overview + phase index)
├── requirements.txt             Python dependencies for analysis + rendering
│
├── pages/                       Site pages (rendered into _site/pages/)
│   ├── initial_loeuf.qmd        Baseline LOEUF by functional category & disease class
│   ├── kegg.qmd                 Phase 0: betweenness centrality + KEGG vs. LOEUF
│   ├── gtex.qmd                 Phase 1: PhyloP conservation + GTEx tissue breadth
│   ├── ancestry_loeuf.qmd       Phase 1.3: gnomAD v4 ancestry-stratified LOEUF
│   └── popgen_selection.qmd     Phase 2: π and PBS across 5 populations
│
├── analysis/                    Production analysis scripts
│   ├── fetch_kegg_pathways.py   Fetches KEGG pathway membership per gene
│   ├── fetch_phylop_scores.py   Fetches PhyloP 100-way conservation scores
│   ├── phase0_*.py              Network position + KEGG analysis
│   ├── phase1_*.py              PhyloP, GTEx, gnomAD ancestry analyses
│   ├── phase2_*.py              π, PBS, network selection figures
│   ├── gene_list_overlap_analysis.py  Cross-source pigmentation gene overlap
│   ├── generate_figures.py      Helper for batch figure generation
│   ├── notebooks/               Exploratory Jupyter notebooks (not run by deploy)
│   └── cluster/                 SLURM/HPC scripts for population genomics on Great Lakes
│
├── data/                        Input data + processed CSVs (tracked)
│   ├── 13104_2015_1128_MOESM2_ESM.xlsx   Raghunath 2015 supp data
│   ├── baxter2018_*.xlsx                 Baxter 2018 pigmentation gene list
│   ├── bajpai2023_*.xlsx                 Bajpai 2023 CRISPR screen
│   ├── gnomad_*.tsv                      gnomAD constraint metrics
│   ├── LOEUF_by_functional_category.xlsx Manually curated LOEUF + categories
│   ├── kegg_pathway_counts.csv           Output of fetch_kegg_pathways.py
│   ├── phylop_scores.csv                 Output of fetch_phylop_scores.py
│   ├── network_constraint_*.csv          Merged datasets per phase
│   ├── gtex_tissue_breadth.csv           GTEx-derived tissue expression breadth
│   ├── gwas_*.csv                        GWAS catalog pigmentation associations
│   └── gene_regions.bed                  Gene coordinates for population genomics
│
├── output/                      Generated figures + summary tables (tracked)
│   ├── figure_phase0_*.{png,pdf}         Network/KEGG vs. LOEUF
│   ├── figure_phase1_*.{png,pdf}         PhyloP, GTEx, ancestry LOEUF
│   ├── figure_phase2_*.{png,pdf}         π, PBS, network selection
│   └── table_*.csv                       Summary tables embedded in site
│
├── docs/                        Project documentation (not part of site)
│   ├── plan.md                  Master analysis plan
│   ├── todo.md                  Active task tracker
│   ├── cluster_request.md       UMich Great Lakes resource request
│   ├── aim1_grant_language.docx Grant language reference
│   └── prompts/                 Saved Claude Code prompts for repeatable tasks
│
└── .github/workflows/
    └── publish.yml              GitHub Actions: render Quarto + deploy to Pages
```

---

## Building the site

### Local preview

```bash
quarto render        # builds _site/ from index.qmd + pages/*.qmd
quarto preview       # live-reload server on http://localhost:port
```

### Deploy

Pushing to `main` triggers `.github/workflows/publish.yml`, which renders the
site on a fresh Ubuntu runner and deploys to GitHub Pages. No manual publish
step required.

---

## Regenerating figures

Figures in `output/` are committed (so the deployed site has them without
running heavy analysis on CI). To regenerate from scratch:

```bash
# Phase 0 — KEGG + network
python analysis/fetch_kegg_pathways.py
python analysis/phase0_network_constraint.py

# Phase 1 — conservation + expression
python analysis/fetch_phylop_scores.py
python analysis/phase1_phylop_analysis.py
python analysis/phase1_gtex_analysis.py
python analysis/phase1_gnomad_v4_ancestry.py

# Phase 2 — population genomics (requires cluster outputs)
python analysis/phase2_pi_pbs_analysis.py
python analysis/phase2_pbs_summary.py
python analysis/phase2_pbs_trees.py
python analysis/phase2_network_selection.py
python analysis/phase2_network_constraint_loeuf.py
```

Phase 2 depends on per-gene π and PBS computed on the UMich Great Lakes cluster
(see `analysis/cluster/`). The resulting `output/pi_per_gene.csv` and
`output/pbs_per_gene.csv` are pulled back to the local repo.

---

## Key findings

- **Pigment-specific genes are significantly more LoF-tolerant** than generic
  signaling genes (Kruskal-Wallis p = 1.70×10⁻⁸), consistent with the recessive
  biology hypothesis: heterozygous LoF carries no fitness penalty.
- **Tolerance ordering is preserved across ancestry groups** (NFE/AFR LOEUF
  concordant), ruling out European ascertainment bias as the sole driver.
- **Tissue expression breadth strongly predicts constraint** (Spearman
  ρ > 0, p < 0.001): broadly expressed genes are more LoF-intolerant
  regardless of pigmentation function.
- **MC1R paradox:** highest LOEUF in the network yet expressed in 53/54 tissues;
  African PBS elevated despite high π — consistent with African purifying
  selection acting on *missense* variants, invisible to LoF metrics.
