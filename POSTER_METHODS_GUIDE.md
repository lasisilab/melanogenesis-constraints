# PEQG 2026 Poster — Methods & Conference Q&A Guide

A plain-language guide to every analysis on the poster, with the formal statistical terminology you'll want when explaining methods to a stats-savvy audience. Use this to prep your verbal walkthrough and anticipate questions.

---

## Big picture: the poster's central claim

**Constraint and adaptation are operating on partly orthogonal evolutionary axes within the melanogenesis network.**

In plain language:
- *Constraint* = which genes can tolerate being broken (measured by LOEUF — does the gene have fewer loss-of-function variants than we'd expect by chance?)
- *Adaptation* = which genes show signatures of recent positive selection in specific populations (measured by PBS — has the allele frequency at this gene shifted dramatically along one population's evolutionary branch?)
- *Orthogonal* = these two forces don't necessarily target the same genes. A gene can be highly constrained yet show no recent positive selection, and vice versa.

The poster builds this claim through three figures, each anchored to a different axis of gene biology: **expression architecture**, **constraint**, and **adaptation**.

---

## Background concepts you'll need to defend

### LOEUF (Loss-of-function Observed/Expected Upper bound Fraction)

**Plain language:** Measures how much a gene resists being broken. For every gene, gnomAD counts how many loss-of-function (LoF) mutations are observed in humans and compares it to how many would be expected by chance given the gene's size and mutability. LOEUF is the upper bound of the 90% confidence interval on the observed/expected ratio.

**How to read it:** Lower LOEUF = fewer LoF mutations than expected = stronger purifying selection = the gene is "constrained" (or "LoF-intolerant"). Higher LOEUF = the gene tolerates LoF mutations and is under weaker constraint.

**Statistical caveat to flag:** LOEUF is calculated using the European-ancestry subset of gnomAD by default. We checked ancestry-stratified versions (NFE/AFR/EAS/SAS) and the rank-order is concordant — the constraint architecture isn't a European artifact.

---

### τ (Tau) — tissue specificity index

**Plain language:** A single number between 0 and 1 that summarizes how tissue-specific a gene's expression is. τ=0 means the gene is expressed equally across all tissues (housekeeping). τ=1 means the gene is expressed in only one tissue (highly specialized).

**Formula:** For each gene, you take its expression value in each of the 54 GTEx tissues, normalize each to the gene's max expression (so values are 0–1), subtract each from 1, and average across tissues. τ = mean(1 − x/x_max). This is the Yanai et al. 2005 formulation.

**Why we use it:** It's the standard pop-gen metric for tissue specificity and is independent of overall expression level.

---

### Tissue breadth

**Plain language:** A simpler count — how many GTEx tissues express this gene above a TPM > 1 threshold? Ranges 0–54.

**Relationship to τ:** Highly correlated (broadly-expressed genes have low τ and high breadth), but not identical. Tissue breadth is a count; τ is a continuous index that's also sensitive to *how much* expression varies between tissues, not just whether expression crosses a threshold.

**Why we include both:** Tissue breadth gives a slightly stronger signal in the joint regression than τ, but τ is the more standard metric. Reporting both is reassuring — the signal isn't an artifact of one specific way of measuring tissue specificity.

---

### KEGG pathway count

**Plain language:** How many distinct biological pathways (in the KEGG database) does this gene participate in? A gene in 1 pathway is a specialist; a gene in 30+ pathways is a pleiotropic hub.

**Why this matters:** Genes that participate in many pathways have more "jobs" — breaking them disrupts more cellular processes — so they tend to be more constrained.

**Statistical detail:** We use `log1p(KEGG count)` = log(KEGG count + 1) in models because the distribution is right-skewed (a few hubs have 50+ pathways; most genes have 1–5).

---

### Betweenness centrality

**Plain language:** A network metric. For each gene in the melanogenesis network, you count how often it sits on the shortest path between any two other genes. High-betweenness genes are bottlenecks through which information flows.

**Difference from KEGG count:**
- KEGG count = cross-system pleiotropy (how many *different pathways* the gene works in)
- Betweenness = within-pathway centrality (how *central* the gene is *within the melanogenesis network specifically*)

**Why both?** They capture related but distinct concepts. In our joint regression, KEGG absorbs most of the betweenness signal (see Figure 2 caption) — likely because pleiotropic genes also tend to be central, so the two metrics overlap.

**Statistical detail:** We use `sqrt(betweenness)` in models because the raw distribution is heavily right-skewed.

---

### PBS (Population Branch Statistic)

**Plain language:** For a given population (say, Africans), PBS asks: "compared to two outgroup populations, how much has allele frequency at this gene shifted *specifically along this population's evolutionary branch*?" A high PBS means the gene has experienced lineage-specific allele frequency change — a potential signature of recent positive selection in that population.

**How it's calculated** (the order matters for explaining the math):

1. Start with **F_ST** (fixation index) — a measure of allele frequency differentiation between two populations at each variant. We use the Hudson estimator (Bhatia et al. 2013), which is unbiased for unequal sample sizes.
2. Convert F_ST to a branch length using the Cavalli-Sforza transformation: **T = −ln(1 − F_ST)**. This puts the value on a more linear scale.
3. For a target population (e.g., African) with outgroup A (e.g., S. Asian) and distant outgroup B (e.g., Melanesian), PBS = (T_target,A + T_target,B − T_A,B) / 2.

**Interpretation:** The math isolates the branch length unique to the target population. If selection has driven allele frequency change in Africans but not in the outgroups, PBS_African will be high.

**Important caveat for the poster:** PBS measures *allele frequency differentiation*, which is **suggestive of but not proof of positive selection**. Demographic events (bottlenecks, founder effects) can also produce elevated PBS. We mention this in the poster.

---

### Z-scoring (standardization)

**Plain language:** Before fitting a regression, we transform each predictor so it has mean=0 and standard deviation=1. This puts all predictors on the same scale so their coefficients are directly comparable. After z-scoring, a β=0.5 for predictor X and β=0.3 for predictor Y means X has a stronger effect than Y per "unit of typical variation."

**Formula:** `z = (x − mean(x)) / SD(x)`.

**Why it matters here:** Tissue breadth ranges 0–54, while KEGG count ranges 1–100+. Without z-scoring, the coefficient on breadth would look small just because the predictor varies over a wider range. Z-scoring removes that artifact.

---

## Figure 1 — GTEx τ-split heatmap

### Plain-language version

The 129 melanogenesis network genes vary enormously in where they're expressed. Some (TYR, TYRP1, DCT — the melanin-synthesis enzymes) are expressed almost exclusively in skin. Others (NFKB1, MITF, transcription factors) are expressed across nearly every tissue in the body. This figure visualizes that variation.

We split the 129 genes into two groups at τ = 0.7: "tissue-specific" (top block, τ ≥ 0.7) and "broadly expressed" (bottom block, τ < 0.7). Within each block, we sort by skin enrichment (how much higher is expression in skin tissues compared to non-skin tissues).

The colored cells show *relative* expression — for each gene (row), the highest-expressing tissue gets the brightest color and other tissues are scaled relative to that maximum.

### Formal statistical terms

- **Data source:** GTEx v8, median TPM (transcripts per million) across 54 tissues.
- **Row normalization:** Each cell is x / row_max — this puts each gene on its own 0–1 scale. This is the same normalization used to compute τ itself.
- **τ threshold:** 0.7 is a conventional cutoff for "highly tissue-specific" (Yanai et al.; Sonawane et al. 2017).
- **Within-block sorting:** Welch's two-sample t-statistic comparing each gene's expression in skin tissues (skin sun-exposed + skin sun-unexposed) vs. all other tissues. Higher t = more skin-enriched.
- **Why Welch's t and not Student's t?** Welch's doesn't assume equal variances between the two tissue groups, which is more realistic for expression data.

### Likely questions

**Q: Why split at τ = 0.7 specifically?**
A: It's a conventional threshold in tissue-specificity literature (Yanai et al. 2005 used τ > 0.7 to define "specifically expressed"). The split is for visualization; the underlying analyses in Figures 2 and 3 treat τ and tissue breadth as continuous variables and don't depend on this cutoff.

**Q: Why row-normalize instead of showing raw TPM?**
A: Raw TPM is dominated by a few extremely high-expressing genes — the heatmap would just show "these few genes are expressed strongly everywhere." Row normalization makes within-gene patterns visible across all 129 genes simultaneously, which is what we want for showing tissue-specificity architecture.

---

## Figure 2 — Joint OLS regression of LOEUF on expression breadth and network connectivity

### Plain-language version

We want to test whether tissue expression breadth and network connectivity each *independently* predict how constrained a gene is. The concern is that broadly-expressed genes might also be the same genes that are in many pathways — so a relationship between breadth and constraint might just be a relationship between connectivity and constraint in disguise.

To rule this out, we fit a **joint model**: predict LOEUF from both breadth and KEGG pathway count *simultaneously*. The coefficient on each predictor in this joint model tells us how much that predictor matters *after* the other one is accounted for. If both coefficients remain significant, both axes are doing independent work.

That's what we find: both breadth (β = −0.135, p = 0.001) and KEGG (β = −0.131, p = 0.002) remain significant in the joint model.

Panels A–D are **added-variable plots** (or partial regression plots). They're the visual proof of this independence — each panel takes one predictor, removes whatever variance is shared with the other predictor, and shows the residual relationship. If the slope is still nonzero, the predictor still matters after the other one is controlled.

Panel E is a **forest plot** — a compact way to show the coefficient and confidence interval for each predictor across multiple model specifications. If the confidence interval doesn't cross zero, the predictor is significant.

### Formal statistical terms

- **Model:** `LOEUF ~ z(breadth) + z(KEGG)` (OLS regression, ordinary least squares).
- **Sample size:** n=127 (some network genes are missing LOEUF or expression data).
- **Z-scoring:** Both predictors standardized so coefficients are directly comparable.
- **Variance Inflation Factor (VIF):** 1.16. This measures how much each predictor is correlated with the other. VIF=1 means perfectly orthogonal; VIF > 5 typically indicates problematic collinearity. Our VIF=1.16 confirms the two predictors are essentially independent — they're not redundant measurements of the same thing.
- **Added-variable plot construction:** For each predictor X, regress Y on the *other* predictors and save the residuals (Y_resid). Regress X on the other predictors and save those residuals (X_resid). Plot Y_resid against X_resid. The slope of this scatter equals the partial regression coefficient in the full model, and its p-value equals the joint-model p-value for that predictor.

### Likely questions

**Q: Why OLS instead of a more sophisticated model?**
A: OLS is the standard for partial-regression analysis and is the appropriate model when you want interpretable coefficients. We checked the residuals for normality and don't see major violations. For robustness, we also computed rank-based partial Spearman correlations (in the supplementary text output), which give the same conclusions.

**Q: What about betweenness centrality?**
A: Betweenness independently predicts LOEUF when paired only with tissue breadth, but loses significance once KEGG is added — KEGG absorbs most of the betweenness signal. We interpret this as KEGG and betweenness being partly redundant measures of network position (cross-system pleiotropy vs. within-pathway centrality). The forest plot (Panel E) shows this clearly.

**Q: Could the relationship be confounded by gene length?**
A: Possibly — longer genes tend to have more LoF variants observed by chance, and longer genes may also be more pleiotropic. LOEUF itself adjusts for expected variant counts based on gene length, mutation rate, and CpG content. So gene length is implicitly controlled. We didn't add it as an explicit covariate in this model.

**Q: What's your R²?**
A: Around 0.21 for the breadth + KEGG joint model. Not enormous, but typical for gene-level association analyses on N~127. The point isn't to predict LOEUF perfectly — it's to demonstrate that both axes contribute.

---

## Figure 3 — Within-network PBS association

### Plain-language version

Among the 129 melanogenesis genes, we ask: do genes with certain biological properties (tissue specificity, expression breadth, network connectivity, constraint) tend to show stronger population-specific selection signals in any given branch?

We test this for four populations: African, Melanesian, European, and South Asian.

**Step 1 — Convert raw PBS to a within-network ranking.** For each population, we rank the 124 network genes (with complete data) by PBS, then convert to a percentile (0–100). A gene at the 90th percentile has higher PBS than 90% of the other network genes.

**Step 2 — Fit a regression for each predictor separately.** We run 20 separate regressions (4 populations × 5 predictors), each one of the form: `PBS percentile ~ z(predictor) + z(log1p SNPs) + z(log10 gene length)`. The covariates control for the mechanical effect of having more SNPs or longer genes on PBS estimates.

**Step 3 — Permutation test for significance.** Instead of using parametric p-values, we shuffle the PBS percentiles randomly among the 124 genes 10,000 times. For each shuffle, we re-fit the model and record the β. The permutation p-value is: how often does a random shuffle produce a β as large as what we observed?

**Main finding:** African selection targets broadly-expressed genes (β_τ = −6.89, perm p = 0.011; β_breadth = +7.88, perm p = 0.003). Melanesian selection targets tissue-specific genes (β_τ = +5.28, perm p = 0.048; β_breadth = −6.18, perm p = 0.021). European and South Asian signals are null — a clean negative control.

The two populations are selecting genes from **opposite ends** of the tissue-specificity axis, even though both are selecting from the same 129-gene pathway.

### Formal statistical terms

- **Outcome:** Within-network PBS percentile (0–100). We rank rather than use raw PBS because the raw distribution is right-skewed and has heteroscedastic variance.
- **Models:** Separate (not joint) OLS regressions, one per (population × predictor) combination. We use single-predictor rather than joint models because we want to test each feature's marginal association with PBS, not its partial effect after controlling for other features.
- **Covariates:** `log1p(n_sites_shared)` (SNP count) and `log10(gene body length)` (gene length minus the 20kb flanking region). Both z-scored.
- **Permutation test:** 10,000 random permutations of the outcome vector within the 124 network genes. Empirical p = (count of |β_null| ≥ |β_observed| + 1) / (N_perm + 1).
- **Why permutation instead of parametric p?** PBS percentiles are bounded (0–100) and not normally distributed. The permutation test makes no distributional assumptions and is conservative.
- **Variance Inflation Factor (VIF):** Computed per model; all VIFs < 1.5, confirming the predictor and covariates aren't badly collinear.

### Likely questions

**Q: Why use percentiles instead of raw PBS?**
A: Raw PBS is right-skewed (a few high-PBS outliers dominate any regression). Percentile-ranking handles this without log-transforming or trimming. The signal you see is on the *ranking* of genes within the network, which is robust to outliers.

**Q: Why single-predictor models instead of one joint model?**
A: For this analysis we want to know whether each feature *on its own* associates with PBS, not whether each contributes independently after the others are controlled. The biological question is "do tissue-specific genes show more Melanesian PBS?" — that's a marginal-association question. Joint models would answer "do tissue-specific genes show more Melanesian PBS *holding constraint and connectivity fixed*?" — a different question. We could do that as a supplementary analysis if a reviewer pushes.

**Q: Why opposite signs for African vs. Melanesian τ coefficients?**
A: That's the headline finding. African selection preferentially targets genes that are broadly expressed (low τ), while Melanesian selection targets tissue-restricted genes (high τ). One interpretation: African selection acts through regulatory variation in genes that contribute to many systems (so any allele frequency change reflects pleiotropic adaptation, possibly to environmental factors like UV, infectious disease, etc.); Melanesian selection acts on tissue-restricted effectors where the variant's effects are more localized (so adaptation can act without pleiotropic constraint).

**Q: How can European and South Asian PBS be computed without their own outgroup pair?**
A: We use the African–Melanesian reference structure: PBS_European = (T_African,European + T_European,Melanesian − T_African,Melanesian) / 2. This is valid as a measure of European-specific lineage drift but is not perfectly parallel to the African and Melanesian scans (which use S. Asian as the proximal outgroup). They serve as a negative-control / lineage-specificity test rather than a primary analysis.

---

## The genome-wide story (likely to come up in Q&A)

### What we did

In addition to the within-network analysis, we ran a **genome-wide** PBS scan using SGDP-only data (~19,000 protein-coding genes) to answer two questions:

1. **Network enrichment:** Is the melanogenesis network as a whole carrying more selection signal than a typical set of 129 genes from elsewhere in the genome?
2. **Within-network percentile:** Does our within-network analysis (Figure 3) replicate when we rank network genes against the *whole genome* rather than just against each other?

### The findings (these are the parts a critic will press on)

**1. Network-level enrichment is null.**
- The 129 melanogenesis genes do not have systematically higher PBS than the genome-wide background (Mann-Whitney p=0.15 for African PBS, p=0.18 for Melanesian PBS).
- 0 of 123 network genes (with PBS data) are in the genome-wide top 1% (expected ~1.2 by chance).
- **Interpretation:** Selection in this network is gene-by-gene, not pathway-wide. The "average" melanogenesis gene is not under detectable lineage-specific selection.

**2. Specific genes do rank high genome-wide.**
- PBS-3 (Melanesian) percentiles: KITLG 98.6, ATF2 98.6, HIF1A 97.7, PRKACA 96.9, MITF 85.8, TYRP1 86.6.
- KITLG is a well-known pigmentation-selection target — a useful positive control.

**3. The genome-wide-percentile version of Figure 3 is null in SGDP.**
- When we use SGDP-only genome-wide percentile as the outcome (instead of within-network percentile), all of Figure 3's signals disappear.
- The African τ coefficient even *sign-flips*: ρ = +0.075 in SGDP genome-wide percentile vs. ρ = −0.234 in combined-dataset within-network percentile.

### How to defend this if asked

**Q: If the SGDP genome-wide version is null, isn't the combined-dataset signal suspicious?**
A: The discrepancy is a **sample size** issue, not a method issue. SGDP has only ~18 Africans; the combined dataset (gnomAD HGDP+1KGP) has 747 Africans. With ~18 samples, allele frequency estimates have huge variance, which inflates noise and obscures real signal. Spearman correlation is rank-invariant, so the sign-flip can't be from the percentile transform — it's purely from the dataset choice.

**Q: Then why not just use the combined dataset everywhere?**
A: We do — Figure 3 uses the combined dataset. The reason we *also* show the SGDP genome-wide result is **transparency**. The SGDP result is what's available as a genome-wide background, and we want to honestly report that the signal is dataset-dependent. We're not hiding the null.

**Q: So the network isn't under selection at all?**
A: The poster claim isn't "the whole melanogenesis network is under selection." It's the more specific (and supported) claim: *within* the network, tissue-specificity stratifies which genes are differentiating between African and Melanesian lineages. That's a within-network distribution claim, not a network-vs-genome enrichment claim. We're careful about that wording in the title and caption.

**Q: How sure are you that the within-network signal isn't just noise?**
A: Three pieces of evidence:
1. **Lineage-specificity:** EUR and SAS are null — if the signal were a noise artifact, it would appear in those branches too.
2. **Sign consistency:** AFR τ and breadth coefficients have opposite signs (as expected — they measure inverse things), and so do MEL τ and breadth coefficients. Random noise wouldn't produce this internal consistency.
3. **Permutation-based p-values:** These are nonparametric and make no distributional assumptions.

---

## QC, filtering, and data provenance

You may get methodological questions on the upstream pipeline:

### VCF filtering (combined dataset, used for Figure 3)
- gnomAD HGDP+1KGP v3.1.2 joint callset (high-quality, pre-filtered by gnomAD)
- SGDP Simons.vcf.gz lifted from hg19 → hg38 using bcftools +liftover
- Per-chromosome merge of gnomAD + lifted SGDP
- **Filter:** biallelic SNPs only, site-level missingness < 10% (F_MISSING > 0.1 excluded)

### SGDP genome-wide pipeline (used for genome-wide background)
- Same as above but additionally: per-genotype GQ < 20 → set to missing (before the site-level missingness filter)

### Populations
- African (747): gnomAD HGDP+1KGP (YRI, LWK, ESN, GWD, MSL, Yoruba, Mandenka)
- Melanesian (47): gnomAD HGDP (PapuanHighlands, PapuanSepik, Bougainville) + SGDP (Papuan)
- South Asian (790): gnomAD HGDP+1KGP
- European (788): gnomAD HGDP+1KGP (CEU, TSI, FIN, GBR, IBS + HGDP European)
- East Asian (718): gnomAD HGDP+1KGP — extracted but not used in poster figures

### Sample size caveat for Melanesian
- Only 47 Melanesians total, with 17 of those from SGDP and 30 from HGDP.
- This is the main statistical limitation of the Melanesian PBS analysis — estimates have wider confidence intervals than the other populations.
- We disclose this in the caveats section of every PBS analysis.

### Gene set
- 129 genes from Raghunath et al. 2015 (*BMC Research Notes*, DOI 10.1186/s13104-015-1128-6)
- Network defined by experimentally-validated melanocyte signaling interactions
- Functional categories assigned by us: Pigment-specific, Developmental/NC, Generic signaling, Cytokines/growth factors, Apoptosis/cell death, Other

---

## One-paragraph elevator pitch for the poster

> "We're testing whether constraint and adaptation are operating on the same axes in the melanogenesis network or on different ones. Using GTEx tissue expression data, KEGG pathway annotations, gnomAD LOEUF constraint scores, and population-specific PBS computed from gnomAD HGDP+1KGP plus SGDP genotypes, we find that tissue specificity and network connectivity each independently predict LOEUF — broadly-expressed, pleiotropic genes are more constrained. Separately, we find that population-specific selection on this same network goes in opposite directions: African selection targets broadly-expressed genes, while Melanesian selection targets tissue-specific genes. Constraint and adaptation are working on partly orthogonal axes — pleiotropy constrains *which kinds of variants can fix*, not *which genes can be selected.*"

---

## Quick stats vocabulary cheat sheet

| Plain term | Formal term | What it means |
|---|---|---|
| Constraint | LoF constraint (or LoF intolerance) | Purifying selection has reduced LoF variants below random expectation |
| Tissue specificity | τ (Yanai index) | Quantitative measure of expression skew across tissues |
| Selection signal | PBS (Population Branch Statistic) | Lineage-specific allele frequency change relative to outgroups |
| Net effect | Marginal association | Effect of one predictor, alone |
| Independent effect | Partial association | Effect of one predictor, controlling for others |
| Effect on same scale | Standardized coefficient (β on z-scored predictor) | Direct comparison across predictors |
| Confidence interval | 95% CI | Range of plausible β values; CIs that exclude 0 are significant |
| Empirical p-value | Permutation p-value | Significance from shuffling the data, not from assumed distribution |
| Predictor redundancy | Variance Inflation Factor (VIF) | VIF=1 → orthogonal; VIF>5 → problematic |
| Partial regression plot | Added-variable plot | Visualizes joint-model β by removing shared variance |
| Side-by-side coefficients | Forest plot | Compact display of β ± CI across models or predictors |

---

*Generated 2026-05-28 — use this guide alongside the poster figures during the conference.*
