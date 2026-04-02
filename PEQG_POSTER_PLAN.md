# PEQG Poster Plan: Testing Network-Based Predictions of Genetic Constraint in Melanogenesis

**Conference:** PEQG, June 2026
**Last updated:** 2026-04-02

---

## What the abstract commits us to

The abstract names four constraint/selection metrics (LOEUF, PhyloP, π, PBS/iHS), two population resources (1000 Genomes, HGDP), a population comparison (African vs. Melanesian), a tissue expression analysis, and three explicit hypotheses. Below is the full inventory of what exists in the repo and what remains, organized as a work plan.

---

## PHASE 0 — Already complete

These analyses exist in the repo and can go directly onto the poster.

### 0.1 LOEUF × functional category
- **Notebook:** `melanogenesis_network_constraint_v2.ipynb`, `generate_figures.py`
- **What it does:** 130 Raghunath network genes classified into 5–6 functional categories; LOEUF compared across categories with Kruskal-Wallis + Dunn's posthoc
- **Key result:** Pigment-specific enzymes (TYR, TYRP1, DCT, etc.) have median LOEUF = 1.89; developmental/NC and signaling genes are constrained (median ~0.30–0.36)
- **Output:** `network_constraint_data.csv`, 3-panel composite figure (centrality scatter, KEGG pathway count scatter, category boxplot)
- **Status:** ✅ Done

### 0.2 LOEUF × disease class
- **Notebook:** `melanogenesis_network_constraint_v2.ipynb`
- **What it does:** Classifies genes by clinical phenotype (syndromic vs. isolated pigment); Mann-Whitney U test
- **Key result:** Complete separation — syndromic max LOEUF = 0.66, isolated pigment min LOEUF = 0.86 (p = 6.78e-05)
- **Status:** ✅ Done

### 0.3 Betweenness centrality × LOEUF
- **Notebook:** `melanogenesis_network_constraint_v2.ipynb`, `generate_figures.py`
- **What it does:** Scatter of network betweenness centrality vs. LOEUF with Spearman correlation
- **Key result:** Topologically central genes are generic signaling hubs, not melanogenesis effectors
- **Status:** ✅ Done

### 0.4 GWAS coding vs. regulatory constraint
- **Notebook:** `pigmentation_constraint_analysis.ipynb`
- **What it does:** Pulls pigmentation GWAS associations from GWAS Catalog API; classifies variants as coding vs. regulatory; compares LOEUF
- **Key result:** Coding GWAS variants sit in more tolerant genes (median LOEUF 0.86) than regulatory variants (0.60), p = 0.02
- **Output:** `pigmentation_gwas_constraint.csv`, figure
- **Status:** ✅ Done

### 0.5 LOEUF by clinical phenotype classification
- **Notebook:** `melanogenesis_constraint_analysis.ipynb`
- **What it does:** Curated pigmentation-restricted vs. pleiotropic/syndromic gene list; LOEUF comparison + GWAS overlay
- **Key result:** Pigmentation-restricted genes are significantly more LoF-tolerant than pleiotropic genes
- **Output:** `supplementary_table_genes.csv`, 2-panel figure
- **Status:** ✅ Done

### 0.6 Gene list overlap
- **Script:** `gene_list_overlap_analysis.py`
- **What it does:** Overlap of Raghunath (134 genes), Baxter (635 genes), Bajpai CRISPR (169 hits) — pairwise and three-way
- **Key result:** Only TYR, DCT, OCA2 in all three; union = 873 unique genes
- **Output:** `gene_list_overlap_results.json`, `gene_list_union_annotated.csv`, Venn diagram
- **Status:** ✅ Done

---

## PHASE 1 — Low-hanging fruit (days of work each)

### 1.1 PhyloP evolutionary conservation scores
- [ ] **Get data:** Download per-gene PhyloP scores (UCSC PhyloP 100-way or 470-way vertebrate conservation). Can use precomputed gene-level summaries (e.g., mean PhyloP across coding exons) from sources like Ensembl BioMart, or compute from bigWig + BED intervals.
- [ ] **Merge:** Join PhyloP scores to the existing 130-gene `network_constraint_data.csv` by gene symbol or Ensembl ID.
- [ ] **Analyze:** Correlate PhyloP with LOEUF (expect strong correlation). Compare PhyloP across functional categories (same Kruskal-Wallis framework). Scatter: PhyloP vs. LOEUF to show convergent evidence.
- [ ] **Figure:** Add PhyloP panel to poster (boxplot by category, or overlay on existing scatter).
- **Why it matters:** Abstract explicitly lists PhyloP as a constraint metric. Quick win that strengthens the "multiple lines of evidence" story.
- **Depends on:** Nothing — independent of population data.

### 1.2 GTEx tissue expression breadth
- [ ] **Get data:** Download GTEx v8 gene-level median TPM per tissue from GTEx Portal (file: `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz`). This is a single downloadable file.
- [ ] **Compute tissue breadth:** For each gene, count number of tissues with TPM > 1 (or similar threshold). This is a proxy for pleiotropy.
- [ ] **Merge:** Join to the 130-gene network data.
- [ ] **Analyze:** Correlate tissue breadth with LOEUF (expect: broader expression → lower LOEUF). Compare tissue breadth across functional categories.
- [ ] **Hypothesis 2 test:** Do genes with broader expression show distinct constraint? Regression: LOEUF ~ tissue_breadth + functional_category.
- **Why it matters:** Hypothesis 2 is about tissue expression × constraint. Without this, you can't address Hypothesis 2 at all.
- **Depends on:** Nothing — independent of population data.

---

## PHASE 2 — The population genomics layer (the heavy lift)

This is the bottleneck. Everything in the second half of the abstract depends on having population-level data for African and Melanesian populations.

### 2.1 Obtain and process 1000 Genomes / HGDP data
- [ ] **Identify populations:**
  - African: YRI (Yoruba), LWK (Luhya), ESN (Esan), GWD (Gambian), MSL (Mende) from 1KGP; and/or African populations from HGDP
  - Melanesian: Bougainville and/or Papuan from HGDP (not in 1KGP)
  - Outgroup for PBS: e.g., CEU or CHB from 1KGP
- [ ] **Download VCFs:** 1KGP phase 3 VCFs are available from the IGSR FTP. HGDP data available from gnomAD (HGDP+1KGP callset) or the Bergström et al. 2020 release.
  - *Decision needed:* Use the gnomAD HGDP+1KGP joint callset? This would simplify getting both African and Melanesian populations in one dataset with consistent variant calling.
- [ ] **Define gene regions:** For each of the 130 network genes, define genomic intervals (gene body ± some flanking region, e.g., ±10kb). Use Ensembl/UCSC gene coordinates.
- [ ] **Extract per-gene VCFs:** Subset the full VCFs to just the gene regions of interest. Tools: `bcftools view -R regions.bed`.
- **Why it matters:** This is the prerequisite for π, PBS, and iHS.
- **Estimated effort:** 1–3 days depending on data access and pipeline familiarity.

### 2.2 Nucleotide diversity (π) per gene per population
- [ ] **Compute π:** For each gene region, compute per-site nucleotide diversity for African and Melanesian populations separately. Tools: `pixy`, `scikit-allel`, or `vcftools --site-pi`.
- [ ] **Normalize:** Consider using π/divergence or π relative to genome-wide expectation to control for mutation rate variation.
- [ ] **Merge:** Join π values to the 130-gene dataset, per population.
- [ ] **Analyze:**
  - Compare π across functional categories within each population
  - Test African vs. Melanesian π (paired by gene) — do the same genes show reduced diversity in both?
  - Correlate π with LOEUF and with network centrality
- [ ] **Figure:** Panel showing π by functional category, split by population.
- **Depends on:** 2.1

### 2.3 PBS (Population Branch Statistic)
- [ ] **Compute allele frequencies:** Per-SNP allele frequencies for African, Melanesian, and outgroup (e.g., CEU) populations.
- [ ] **Compute PBS per gene:** PBS requires F_ST between three population pairs. Summarize per gene (e.g., max PBS SNP per gene, or mean PBS).
  - Tools: Custom script using Hudson's F_ST estimator, or `scikit-allel`.
- [ ] **Merge:** Add per-gene PBS to the dataset.
- [ ] **Analyze:**
  - Which genes show elevated PBS in African? In Melanesian?
  - Are PBS outliers concentrated at peripheral (pigment-specific) network positions? (Hypothesis 3)
  - Scatter: PBS vs. network centrality, colored by functional category
- [ ] **Figure:** PBS Manhattan-style or scatter panel for poster.
- **Depends on:** 2.1

### 2.4 iHS (integrated haplotype score)
- [ ] **Prepare phased data:** Need phased haplotypes. 1KGP phase 3 data is already phased. HGDP may need phasing (SHAPEIT/Eagle).
- [ ] **Compute iHS:** Use `selscan` or `hapbin` to compute iHS for African and Melanesian populations.
- [ ] **Summarize per gene:** Max |iHS| per gene, or proportion of SNPs with |iHS| > 2.
- [ ] **Merge and analyze:** Same framework as PBS — overlay on network, test Hypothesis 3.
- **Note:** This is the most computationally intensive metric. If time is short, PBS alone may suffice for the poster, with iHS noted as "in progress."
- **Depends on:** 2.1 + phased data

---

## PHASE 3 — Hypothesis testing

### 3.1 Hypothesis 1: Centrality predicts constraint regardless of population
- [ ] Test: Does network centrality (betweenness) correlate with constraint (LOEUF, PhyloP, π) consistently across African and Melanesian populations?
- [ ] Method: Spearman correlations of centrality vs. each constraint metric, per population. Compare correlation coefficients.
- [ ] This is partially addressed by existing LOEUF × centrality analysis (Phase 0.3), but the abstract promises population-level evidence.
- **Depends on:** Phase 1.1 (PhyloP) + Phase 2.2 (π)

### 3.2 Hypothesis 2: Tissue breadth predicts distinct constraint patterns
- [ ] Test: Do broadly expressed genes show different constraint than tissue-specific genes, after controlling for functional category?
- [ ] Method: Multiple regression — LOEUF ~ tissue_breadth + functional_category + centrality. Or split genes by tissue breadth tertile and compare LOEUF.
- **Depends on:** Phase 1.2 (GTEx)

### 3.3 Hypothesis 3: Population-specific selection at peripheral positions
- [ ] Test: Do PBS/iHS outliers concentrate at peripheral (low-centrality, pigment-specific) network positions while core pathway genes remain universally constrained?
- [ ] Method: Compare network centrality of selection outlier genes vs. non-outliers. Fisher's exact test: are outlier genes enriched among peripheral genes?
- [ ] Compare African vs. Melanesian outlier gene sets — do they overlap (convergent targets) or differ (independent paths)?
- **Depends on:** Phase 2.3 (PBS) and/or 2.4 (iHS)

---

## PHASE 4 — Poster creation

### 4.1 Figure layout
- [ ] Design poster figure panels (likely 5–7 panels total)
  - Panel 1: Network schematic or gene overlap Venn (context-setting)
  - Panel 2: LOEUF × functional category boxplot (existing)
  - Panel 3: PhyloP × functional category (or LOEUF vs. PhyloP scatter)
  - Panel 4: π and/or PBS by functional category, split by population
  - Panel 5: Tissue breadth × constraint (Hypothesis 2)
  - Panel 6: Selection outliers × network position (Hypothesis 3)
- [ ] Unify color scheme across all panels (use existing palette from `generate_figures.py`)

### 4.2 Draft poster
- [ ] Write intro/motivation text
- [ ] Write methods summary
- [ ] Assemble figures
- [ ] Write conclusions

---

## Priority order for bare-minimum poster

| Priority | Task | Effort | Why essential |
|----------|------|--------|---------------|
| 1 | PhyloP scores (1.1) | ~1 day | Named in abstract, easy, second constraint axis |
| 2 | GTEx tissue breadth (1.2) | ~1 day | Required for Hypothesis 2 |
| 3 | Get population data (2.1) | ~2–3 days | Prerequisite for everything population-level |
| 4 | Nucleotide diversity π (2.2) | ~2 days | Easiest population metric, strong story |
| 5 | PBS (2.3) | ~2–3 days | Core of African vs. Melanesian comparison |
| 6 | Hypothesis tests (3.1–3.3) | ~2 days | Abstract promises three hypotheses |
| 7 | iHS (2.4) | ~3–5 days | Can be marked "in progress" if time is short |
| 8 | Poster assembly (4.1–4.2) | ~3 days | The actual deliverable |

**Critical path:** Phases 1.1 and 1.2 can start immediately and in parallel. Phase 2.1 (data acquisition) should start ASAP since it gates everything population-level.

---

## Decisions needed

1. **Which population dataset?** Use the gnomAD HGDP+1KGP joint callset (simplest), or download 1KGP and HGDP separately?
2. **Which African populations?** YRI only, or pool multiple West African groups?
3. **Which Melanesian populations?** Bougainville, Papuan, or both?
4. **PBS outgroup:** CEU? CHB? Both?
5. **Gene region definition:** Gene body only, or ± flanking? How much?
6. **iHS priority:** Include on poster or mark as ongoing?
