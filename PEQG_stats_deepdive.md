# PEQG 2026 — Under-the-Hood Stats & Math Guide

*A companion to POSTER_METHODS_GUIDE.md. That guide is plain-language Q&A. This one is the math layer: for each method, the equation, what it actually measures, why you chose it over alternatives, and the sharpest technical question a stats-savvy person could ask — with your answer ready.*

*Personal prep. Not for the repo.*

---

## How to use this

You will not be asked to derive anything at a poster. The point of knowing the math is **confidence**: when someone asks "why a permutation test?" or "isn't Hudson's FST biased?", you want to answer from understanding, not from a memorized line. Read each section for the *intuition first*, then the formula. If you can explain the intuition in one sentence, you own the method.

A reframe to hold onto: every method on this poster is a **standard, defensible, published choice**. None of them is exotic. The questions will mostly be "why this and not that," and for every one of those there's a clean answer below.

---

## Part 1 — The three statistical ideas that tie everything together

Before the individual methods, four framing distinctions. If you're solid on these, most "gotcha" questions dissolve.

### 1.1 Marginal vs. partial association (why joint for LOEUF, single for PBS)

- A **marginal** association is the effect of a predictor *on its own*: "do tissue-specific genes have higher Melanesian PBS?"
- A **partial** association is the effect *after holding other predictors fixed*: "controlling for connectivity and constraint, do tissue-specific genes have higher PBS?"

For **LOEUF (Figure 2)** you used a **joint** model on purpose — the whole question was *independence*: do breadth and connectivity each predict constraint after accounting for the other? That's a partial-effects question.

For **PBS (Figure 3)** you used **single-predictor** models on purpose — the question is whether each feature *on its own* tracks selection. That's a marginal-effects question. "Do tissue-specific genes show more Melanesian PBS?" is a marginal question by construction.

> **If pressed:** "Different questions need different models. For constraint I'm claiming the two pleiotropy axes are independent, which is a partial-effects claim, so I fit them jointly. For selection I'm asking whether each feature marginally tracks PBS, so single-predictor is the right design. I can run the joint PBS version as a supplement — it doesn't change the story."

### 1.2 Within-network vs. genome-wide (the distinction that protects you)

Figure 3 ranks genes **within the 129-gene network**, not against the genome. So the claim is about the *internal distribution* of selection signal across the network, **not** that the network is enriched for selection versus the rest of the genome.

This matters because your genome-wide enrichment test (SGDP) is **null**: the network as a whole is not above background (Mann-Whitney p = 0.15 African, 0.18 Melanesian; 0 of 123 genes in the genome-wide top 1%). That's fine — and not a contradiction — *because you never claimed network-level enrichment.* Keep these two claims in separate boxes in your head:

- **Within-network distribution** (your claim): which network genes carry more signal, ranked among themselves → tissue specificity stratifies this. ✅ supported.
- **Network-vs-genome enrichment** (NOT your claim): is the whole pathway a selection hotspot? → null, and you report it honestly.

### 1.3 PBS is a relative / outlier statistic

PBS isn't a yes/no selection test. It's a *continuous measure of lineage-specific allele-frequency change*, and it's meaningful mostly in the **tail**. The median PBS across the genome (~0.31 in SGDP) is just the demographic baseline — deep population divergence raises everyone. So "the median is high" is not evidence of selection; the **outliers** are what you interpret. This is why you rank and look at percentiles rather than reading absolute PBS values.

### 1.4 Allele-frequency differentiation ≠ proof of selection

PBS and FST measure *differentiation*. Drift, bottlenecks, and founder effects also produce differentiation. So the honest verb is **"consistent with selection,"** never "evidence of" or "proof of." Say it proactively; it makes you look careful, not weak.

---

## Part 2 — Constraint & gene-property metrics

### 2.1 LOEUF — Loss-of-function Observed/Expected Upper bound Fraction

**What it measures.** How strongly purifying selection has removed loss-of-function (LoF) variants from a gene in the human population. Low LOEUF = fewer LoF variants than expected = strong constraint. High LOEUF = LoF-tolerant.

**The math, conceptually.** gnomAD builds a sequence-context mutation model to predict the **expected** number of LoF variants for a transcript (given its length, base composition, CpG content, and trinucleotide mutation rates). It then counts the **observed** LoF variants in ~125k exomes. The point estimate is the ratio o/e. Because small genes have noisy ratios, gnomAD fits a Poisson model to the observed count given the expected and reports the **upper bound of the 90% confidence interval** on o/e — that upper bound is LOEUF. The "upper bound" is what makes it conservative: a short gene with 0 observed LoF doesn't get an artificially perfect o/e = 0; its CI upper bound stays sensible.

**Why this choice.** LOEUF (gnomAD v2.1.1) is the field-standard continuous constraint metric and is preferred over the older pLI because it's continuous (rankable) rather than a near-binary classifier.

> **Sharpest question:** *"LOEUF is computed on a European-default sample. Aren't your constraint results a European artifact?"*
> **Answer:** "I checked ancestry-stratified LOEUF (Han et al. 2025 — NFE/AFR/EAS/SAS). The rank order is concordant across ancestries, so the constraint architecture isn't a European artifact." (Also note: LoF constraint is dominated by strong purifying selection that's largely shared across human populations.)

> **Second question:** *"LOEUF only sees LoF variants. What about missense or regulatory selection?"*
> **Answer:** "Exactly right, and that's a key limitation I lean on — MC1R is the poster child: highest LOEUF in the set (1.97, totally LoF-tolerant) yet has documented African purifying selection on *missense* variants. LoF metrics are blind to that. It's also why African adaptation on broadly-expressed (constrained) genes is plausible: selection can act through regulatory/missense variation that LOEUF doesn't capture."

### 2.2 τ (tau) — tissue-specificity index (Yanai et al. 2005)

**What it measures.** How tissue-restricted a gene's expression is. τ = 0 → expressed evenly everywhere (housekeeping). τ = 1 → expressed in one tissue only (specialist).

**Formula.** For a gene with expression x_i across N tissues, normalize to the max: x̂_i = x_i / max(x). Then

  τ = Σ (1 − x̂_i) / (N − 1),  summed over all N tissues.

You computed it from log₂(TPM + 1) values across the 54 GTEx tissues.

**Intuition.** If a gene is equally expressed everywhere, every x̂_i ≈ 1, so each (1 − x̂_i) ≈ 0 and τ ≈ 0. If it's expressed in one tissue, all the other x̂_i ≈ 0, so the sum approaches (N − 1), and dividing by (N − 1) gives τ ≈ 1.

**Why this choice.** τ is the most widely used and best-behaved tissue-specificity index (it outperforms entropy- and "tissue-count"-based measures in benchmarking). It's independent of the gene's overall expression level — it only cares about the *shape* of the across-tissue profile.

> **Sharpest question:** *"τ depends on how you handle the expression transform and zeros. Robust?"*
> **Answer:** "I report tissue breadth alongside τ as a second, simpler measure of the same thing, and the joint-regression conclusions hold for both. So the signal isn't an artifact of one specific tissue-specificity definition."

### 2.3 Tissue breadth

**What it measures.** Number of GTEx tissues (of 54) where the gene's median TPM > 1. A plain count, 0–54.

**Relationship to τ.** Strongly (inversely) correlated — broad genes have low τ — but not identical. Breadth is a threshold count; τ is continuous and sensitive to *how much* expression varies, not just whether it crosses TPM > 1. Reporting both is the robustness move.

### 2.4 KEGG pathway count (and why log1p)

**What it measures.** How many distinct KEGG pathways a gene participates in — a proxy for **cross-system pleiotropy**. A gene in 1 pathway is a specialist; a gene in 30+ is a hub.

**Why log1p.** The raw distribution is heavily right-skewed (a few hubs with 50+, most with 1–5). `log1p(x) = log(1 + x)` compresses the long tail and handles zero counts gracefully, so a handful of hub genes don't dominate the regression.

### 2.5 Betweenness centrality (and why sqrt)

**What it measures.** For a gene in the melanogenesis network, the fraction of shortest paths between all other gene pairs that pass through it. High betweenness = a bottleneck / information-flow hub *within this specific network*.

**Formula (intuition).** betweenness(v) = Σ over node pairs (s,t) of σ_st(v) / σ_st, where σ_st is the number of shortest s–t paths and σ_st(v) how many of them go through v.

**KEGG vs. betweenness — the distinction.** KEGG count = *cross-system* pleiotropy (how many different pathways). Betweenness = *within-network* centrality (how central inside the melanogenesis graph). Related but distinct. In your joint model, KEGG absorbs the betweenness signal — betweenness predicts LOEUF on its own but loses significance once KEGG is in, because pleiotropic genes also tend to be central. **sqrt** transform is used because raw betweenness is heavily right-skewed.

---

## Part 3 — Population-genetic statistics

### 3.1 FST — Hudson estimator (Bhatia et al. 2013)

**What it measures.** Allele-frequency differentiation between two populations at a SNP (0 = identical frequencies, 1 = fixed for different alleles).

**Hudson per-SNP formula.** With allele frequencies p₁, p₂ and sample sizes n₁, n₂:

  Numerator: (p₁ − p₂)² − p₁(1−p₁)/(n₁−1) − p₂(1−p₂)/(n₂−1)
  Denominator: p₁(1−p₂) + p₂(1−p₁)

  FST = Σ Numerator / Σ Denominator  (summed over SNPs — a **ratio of sums**, not a mean of per-SNP ratios)

**Why Hudson over Weir-Cockerham.** Hudson's estimator is **unbiased with respect to unequal and small sample sizes** (Bhatia et al. 2013 showed Weir-Cockerham can be biased when sample sizes differ a lot). With n = 47 Melanesians vs 747 Africans, that bias-robustness matters a great deal — it's a deliberate, defensible choice you can name.

**Why ratio-of-sums.** Averaging per-SNP FST ratios is unstable (low-frequency SNPs blow up the ratio). Summing numerators and denominators separately, then dividing, is the standard, stable way to get a regional FST.

> **Sharpest question:** *"With 47 Melanesians, your FST estimates are noisy."*
> **Answer:** "Yes — small n inflates the *variance* of FST, which widens confidence intervals on Melanesian PBS, but the Hudson estimator is specifically chosen because it's *unbiased* under unequal/small samples, so it doesn't shift the point estimate. I disclose n = 47 as the main limitation."

### 3.2 The Cavalli-Sforza / Felsenstein transform: T = −ln(1 − FST)

**What it does.** Converts an FST value into an (approximate) **branch length** — a quantity that grows roughly linearly with divergence time under pure drift. FST itself saturates toward 1 and is non-linear in time; −ln(1 − FST) straightens it out so branch lengths can be added.

**Why you need it.** PBS adds and subtracts branch lengths between three populations. You can only do that arithmetic on a linearized (additive) scale, which is exactly what T provides.

### 3.3 PBS — Population Branch Statistic (Yi et al. 2010)

**What it measures.** For a target population A relative to two outgroups B and C, how much allele-frequency change is **specific to A's own branch** — a candidate signature of recent positive selection in A.

**Formula.** Compute the three pairwise FSTs, transform each to T, then:

  PBS_A = (T_AB + T_AC − T_BC) / 2

**Intuition (draw the triangle).** Three populations form a tree with three leaves. T_AB and T_AC are the two branches touching A (each is A's branch + an outgroup branch); T_BC is the branch between the two outgroups (no A in it). Adding T_AB + T_AC double-counts A's own branch and singly-counts each outgroup branch; subtracting T_BC removes the outgroup-only contribution; dividing by 2 isolates **A's branch alone**. A long A-branch = lineage-specific differentiation at that locus.

**Your scans (target / proximal outgroup / distant outgroup):**
- PBS-1: African / South Asian / Papuan
- PBS-3: Melanesian / South Asian / African

> **Sharpest question:** *"How do you get European and South Asian PBS without their own outgroup pair?"*
> **Answer:** "For EUR/SAS I reuse the African–Melanesian reference structure: PBS_EUR = (T_AFR,EUR + T_EUR,MEL − T_AFR,MEL)/2. That's a valid measure of European-specific lineage drift, but it's not perfectly parallel to the African/Melanesian scans, which use S. Asian as the proximal outgroup. So I treat EUR/SAS as a lineage-specificity **negative control**, not a primary scan — and they come out null, which is the point."

> **Second question:** *"PBS can't distinguish selection from demography."*
> **Answer:** "Correct — PBS measures lineage-specific allele-frequency change, which drift and bottlenecks can also produce. That's why I say 'consistent with selection,' and why the controls matter: EUR/SAS null, opposite AFR/MEL signs, and permutation-based significance. A pure demographic artifact wouldn't produce that structure."

---

## Part 4 — The regression & inference machinery

### 4.1 OLS joint regression (Figure 2) and z-scoring

**Model.** LOEUF ~ β₀ + β₁·z(breadth) + β₂·z(KEGG), ordinary least squares, n = 127.

**Z-scoring.** Each predictor transformed to mean 0, SD 1: z = (x − mean)/SD. This puts breadth (range 0–54) and KEGG (range 1–100+) on the **same scale**, so β₁ and β₂ are directly comparable as "change in LOEUF per 1 SD of predictor." Without it, breadth's coefficient would look artificially small just because it varies over a narrower numeric range.

**Your numbers.** β_breadth = −0.135 (p = 0.001), β_KEGG = −0.131 (p = 0.002). Both negative → more breadth / more pathways → lower LOEUF → more constrained. R² ≈ 0.21.

> **"R² is only 0.21."** "Right — and that's typical for gene-level association analysis at n ≈ 127. The claim isn't that I can predict LOEUF; it's that both axes contribute independently. The coefficients and their significance are the result, not the R²."

### 4.2 VIF — Variance Inflation Factor

**What it measures.** How much a predictor is collinear with the others. **VIF_j = 1 / (1 − R²_j)**, where R²_j comes from regressing predictor j on the remaining predictors.

**Your value.** VIF = 1.16 → R²_between-predictors ≈ 0.14. Rule of thumb: VIF = 1 is perfectly orthogonal, VIF > 5 (or > 10) signals problematic collinearity. **1.16 means breadth and KEGG are essentially independent** — they're not two names for the same variable. This is the statistical backbone of your "two independent axes of pleiotropy" claim.

### 4.3 Added-variable (partial regression) plots — Panels A–D

**What they show.** The relationship between Y and one predictor *after removing the variance both share with the other predictors*. Construction:

1. Regress Y (LOEUF) on the *other* predictor(s), keep residuals e_Y.
2. Regress the predictor of interest X on the *other* predictor(s), keep residuals e_X.
3. Scatter e_Y vs e_X.

The **slope of that scatter equals the full-model partial coefficient** for X, and its p-value equals the full-model p-value for X. So a nonzero slope is the *visual proof* that the predictor still matters once the other is accounted for. That's why these panels are your evidence of independence, not just decoration.

### 4.4 Forest plot — Panel E (and Figure 3)

**What it shows.** Each row is one model/predictor; the dot is the standardized β, the whiskers are the 95% CI. **If the CI crosses zero, the predictor is not significant.** A compact way to show, at a glance, that breadth and KEGG exclude zero while (e.g.) betweenness's CI crosses zero once KEGG is in the model.

### 4.5 Within-network PBS percentile regression (Figure 3)

**Outcome.** For each population, rank the 124 genes (complete data) by PBS and convert to a **percentile (0–100)**. Why rank: raw PBS is right-skewed with a few dominating outliers; ranking is robust without log-transforming or trimming, and isolates the *ordering* of genes, which is what you care about.

**Model (per population × predictor).** percentile ~ z(predictor) + z(log1p n_snps) + z(log10 gene length). The two covariates absorb the mechanical effect that genes with more SNPs or longer bodies can have differently-behaved PBS estimates.

**Your numbers.** AFR τ: β = −6.89 (perm p = 0.011); AFR breadth: β = +7.88 (perm p = 0.003); MEL τ: β = +5.28 (perm p = 0.048); MEL breadth: β = −6.18 (perm p = 0.021); EUR/SAS null; KEGG/betweenness/LOEUF ns everywhere.

**Internal consistency to point at:** AFR τ and breadth have **opposite signs** (as they must — they measure inverse things), and so do MEL τ and breadth. Random noise wouldn't manufacture that sign structure. The AFR and MEL τ signs are *also* opposite each other — that's the headline.

### 4.6 Permutation test — why, and the formula

**Procedure.** Shuffle the outcome (PBS percentiles) randomly across the 124 genes 10,000 times. Refit the model on each shuffle and record the null β. The empirical p-value is

  p = (#{ |β_null| ≥ |β_observed| } + 1) / (N_perm + 1)

(The "+1" in numerator and denominator prevents p = 0 and keeps the test valid.)

**Why permutation, not parametric.** PBS percentiles are bounded (0–100) and not normally distributed, so parametric t-test assumptions are shaky. Permutation makes **no distributional assumption** and is conservative. In practice your adjusted parametric p and permutation p agree closely, which is reassuring.

> **"Why permutation instead of parametric?"** "The outcome is a bounded percentile, not normal, so I didn't want to lean on normal-theory p-values. Permutation is assumption-free and conservative, and it agrees with the parametric p anyway."

### 4.7 Mann-Whitney U (network enrichment test)

**What it does.** A rank-based (nonparametric) test of whether one distribution is stochastically shifted from another — here, network-gene PBS vs. genome-wide background PBS. You used it **one-sided** (are network genes *higher*?). Result: p = 0.15 (African), 0.18 (Melanesian) → no enrichment. You also ran a **SNP-count-matched permutation null** to control for the gene-length/SNP-count confound, which agreed.

### 4.8 Spearman ρ — and why it kills the "percentile artifact" worry

**What it is.** Pearson correlation computed on **ranks** rather than raw values, so it captures monotonic association and is robust to skew and outliers.

**Why it matters for your SGDP story.** Spearman is **rank-invariant** — it only sees the ordering. So when AFR τ shows ρ = −0.234 in the combined data but ρ = +0.075 (sign-flip) in SGDP-only, that difference **cannot** be caused by the percentile transform (which preserves ranks). It can *only* be a **dataset** effect: SGDP has ~18 Africans vs 747 combined. This is your airtight rebuttal to "your percentile ranking created the signal."

### 4.9 Fisher's exact test (τ-split × selection contingency)

**What it does.** Tests association in a 2×2 table when counts are small (exact, not chi-square approximation). Your table: {tissue-specific vs broad} × {African-targeted vs Melanesian-targeted} among population-specifically-selected genes. Result: odds ratio = 0.10, p = 0.046 — the selection axis aligns with the tissue-specificity axis. The Mann-Whitney τ contrast (p = 0.032) and LOEUF non-difference (p = 0.97) are the supporting pair.

---

## Part 5 — The questions most likely to rattle you, pre-answered

**"Isn't the whole network just under selection?"**
No — network-level enrichment is null (MWU p = 0.15–0.18; 0/123 in genome-wide top 1%). My claim is the within-network *distribution*, not network-vs-genome enrichment. Individual effectors (KITLG, ATF2, HIF1A) do rank in the genome-wide top 2–3%, but those are specific genes, not the pathway.

**"If the SGDP genome-wide version is null, why trust the combined signal?"**
Sample size. SGDP has ~18 Africans; combined has 747. Spearman is rank-invariant, so the sign-flip isn't a transform artifact — it's the data. I report both; I'm not hiding the null.

**"Why ranks / percentiles instead of raw PBS?"**
Raw PBS is right-skewed; a few outliers dominate any regression. Ranking is robust without trimming or logging.

**"Why single-predictor for PBS but joint for LOEUF?"**
Different questions — marginal vs. partial. Constraint is an independence claim (joint); selection is a "does each feature track PBS on its own" claim (marginal). Joint PBS available as supplement.

**"τ might just be proxying constraint."**
In the earlier raw-PBS model, AFR τ survives LOEUF control (τ p = 0.022, LOEUF itself ns). And within selected genes, τ differs (p = 0.032) while LOEUF doesn't (p = 0.97). So τ isn't standing in for constraint.

**"n = 47 Melanesians is too small."**
It's the main limitation and I disclose it. It widens CIs on Melanesian PBS but, with the Hudson estimator, doesn't bias the point estimates. The Melanesian τ signal (perm p = 0.048) is real but the thinnest; I don't overstate it.

**"What about admixture in the Melanesian samples?"**
Real concern. HGDP samples were selected for minimal admixture; SGDP includes some admixed individuals. PBS partly accounts for shared variation via the outgroup structure, but for publication I'd want admixture-aware methods.

**"Background selection / mutation-rate variation as confounders?"**
Known omitted covariates. The within-network percentile partly controls for them because all 129 genes sit in similar gene-dense regions; I acknowledge it as a limitation.

**"Have you done haplotype tests (iHS)?"**
No — iHS needs phased whole-chromosome data and genome-wide normalization, incompatible with targeted extraction. Deferred to future work.

---

## Part 6 — When you don't know

You will get a question you can't answer. This is **expected and normal**, even for senior people. Good moves:

- "That's a great question — I haven't tested that directly. My instinct is [X], but I'd want to check before claiming it."
- "We haven't done that yet; it's on the list for the publication version."
- "I'm not sure — can you say more about what you'd expect to see?" (turns it into a conversation, buys time, and senior people *love* being asked to elaborate).
- "Let me write that down." (Then actually write it down. It signals you take the input seriously and ends the exchange gracefully.)

None of these is a failure. "I don't know" delivered calmly reads as scientific maturity. Bluffing is the only real mistake.

---

*Generated for PEQG 2026 prep. Re-verify any number against the latest CSVs before it goes on the poster.*
