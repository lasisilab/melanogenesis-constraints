# Within-network PBS association analysis — execution plan & progress

**Goal:** Among the 129 melanogenesis-network genes, test whether gene-level
features (tissue specificity, network centrality, pathway involvement, LOEUF,
expression breadth) predict **branch-specific genome-wide-percentile PBS** for
each focal branch (AFR, MEL, EUR, SAS).

This is a **within-network association** analysis (per the steering note), not a
network-vs-genome test. Genome-wide data is used only to convert raw PBS into
branch-specific percentiles. The `PBS ~ feature * network` interaction model is
explicitly shelved.

**Defaults chosen (2026-05-27):**
- Include `tissue_breadth` as a 5th predictor (continuity with constraint work).
- Centrality = `betweenness_centrality` only, chosen a priori (no graph recompute).
- Gene length = unpadded body length = `(end − start) − 20000` (strip ±10 kb pad).

---

## RESULT (2026-05-27): null across the board — and a dataset-robustness flag

On the **SGDP-only genome-wide percentile** outcome, **no predictor is
significant in any branch** (permutation p > 0.13 everywhere; all 95% CIs cross
zero). This **conflicts with** the earlier combined-dataset results:

| | combined HGDP+1KGP+SGDP | SGDP-only genome-wide |
|---|---|---|
| AFR τ vs PBS (Spearman) | ρ=−0.234, p=0.009 | ρ=+0.075, p=0.41 |
| AFR breadth vs PBS | ρ=+0.273, p=0.002 | ρ=−0.075, p=0.41 |

Spearman is rank-invariant, so this is **not** the percentile transform — it's
the **dataset**. Combined African PBS uses gnomAD's hundreds of African/S.Asian
samples; SGDP-only uses ~18 SGDP Africans + 17 Melanesians. The African τ signal
(which drove the "adaptation exploits tissue specificity" claim, Fisher p=0.046
and raw-PBS regression p=0.014 — both on the combined dataset) **does not
replicate, and sign-flips, in SGDP-only.**

**Decisive next step:** rerun this exact percentile analysis on the **HGDP
genome-wide** scan (larger samples, being built on the cluster):
`--prefix pbs_genomewide_chr`. If HGDP recovers the African signal → SGDP-only
was underpowered and the combined result stands. If HGDP is also null → the
original signal is dataset-specific and should be reported cautiously.

## Progress checklist

- [x] **Step 1 — Assemble data.** Concatenate `pbs_sgdp_genomewide_chr{1..22}.csv`
      (~19 k genes); merge network features (`network_constraint_categorized.csv`,
      `kegg_pathway_counts.csv`).
- [x] **Step 2 — Gene-level PBS, 4 branches, genome-wide.** AFR=`pbs1_african`,
      MEL=`pbs3_melanesian`; EUR/SAS computed from FST (Afr–Mel reference).
- [x] **Step 3 — Branch-specific genome-wide percentiles.** Rank each branch's
      gene PBS (finite, n_snps>0) → 0–100 percentile.
- [x] **Step 4 — Subset to 129 network genes.**
- [x] **Step 5 — Predictors & covariates.** Predictors: τ, betweenness,
      log1p(KEGG), LOEUF, tissue_breadth. Covariates: log1p(n_snps),
      log10(body_length). All z-scored.
- [x] **Step 6 — Three tests per branch × predictor.** (a) Spearman ρ;
      (b) adjusted single-predictor OLS (β, 95% CI, p, VIF);
      (c) within-network permutation (10 k shuffles → empirical p).
- [x] **Step 7 — Outputs.** `phase2_pbs_within_network.txt`,
      `phase2_pbs_within_network.csv`,
      `figure_phase2_pbs_within_network.{png,pdf}` (2×2 forest, one panel/branch).
- [x] **Step 8 — Caveats baked into text output.** Afr–Mel reference for EUR/SAS;
      n=17 MEL; PBS = differentiation not selection; EUR = contrast/positive
      control; MEL admixture/drift; omitted covariates (recomb/GC/mappability/B).
- [x] **Review** figure + table; report findings.

## Notes / deviations
(to be filled during execution)

## Omitted-by-design
- Recombination rate, GC content, mappability, background-selection score
  (not on disk; future controls).
- Interaction model `PBS ~ feature * network` (shelved per steering note).
- degree / closeness / eigenvector centrality (would require recomputing the
  Raghunath graph; betweenness used a priori instead).
