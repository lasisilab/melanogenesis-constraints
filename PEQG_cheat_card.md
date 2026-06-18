# PEQG 2026 — Poster-Side Cheat Card

*One page. Glance, don't read. Keep at the poster.*

---

## The thesis in one breath

> Pleiotropy (broad expression + connectivity) predicts **constraint**. Tissue specificity predicts **where adaptation acts** — African selection on broadly-expressed genes, Melanesian on tissue-restricted ones. LoF constraint doesn't predict the adaptive target (p = 0.97). **Pleiotropy constrains the *type* of variation that can fix, not which genes can adapt.**

---

## The three figures

| Fig | Shows | Headline stat |
|---|---|---|
| **1** Tissue-specificity gradient | 129 genes, τ-split heatmap; specific enzymes (TYR/TYRP1/DCT) on top, broad regulators (MITF/KIT/EDNRB) below | descriptive — sets the stage |
| **2** Pleiotropy → constraint | Joint OLS, breadth + KEGG predict LOEUF, both independent | β_breadth = −0.135 (p=.001); β_KEGG = −0.131 (p=.002); **VIF = 1.16** |
| **3** τ → population-specific PBS | Within-network percentile, 4 branches | AFR τ −6.89 (p=.011); MEL τ +5.28 (p=.048); **opposite signs**; EUR/SAS null |

---

## Numbers that anchor you

- **Gene set:** 129 (Raghunath 2015); n = 127 in LOEUF regression; 124 regions in PBS.
- **Populations:** AFR 747 · MEL **47** (30 HGDP + 17 SGDP) · SAS 790 · EUR 788.
- **MC1R paradox:** LOEUF = 1.97 (most LoF-tolerant) yet African missense selection — LoF metrics miss it.
- **Supporting stats:** Fisher OR = 0.10, p = 0.046 · τ contrast p = 0.032 · LOEUF non-difference p = 0.97.
- **Genome-wide percentiles (Melanesian):** KITLG 98.6 · ATF2 98.6 · HIF1A 97.7 · PRKACA 96.9.

---

## One-line "why this method" answers

- **LOEUF** — field-standard continuous LoF-constraint; checked ancestry-stratified, concordant.
- **τ** — standard tissue-specificity index, level-independent; report breadth too for robustness.
- **Hudson FST** — *unbiased for unequal/small n* (matters at n=47 MEL).
- **PBS** — isolates the target population's own branch from 3 pairwise FSTs.
- **Z-scoring** — puts predictors on one scale so β's are comparable.
- **VIF 1.16** — breadth & KEGG essentially independent → two real axes of pleiotropy.
- **Added-variable plots** — slope = partial coefficient → visual proof of independence.
- **Percentile + permutation** — PBS is skewed & bounded; ranking is outlier-robust, permutation is assumption-free.
- **Spearman (rank-invariant)** — proves the SGDP sign-flip is a *dataset* effect, not the percentile transform.

---

## Disclose these proactively (it reads as careful, not weak)

1. **"Consistent with selection," never "evidence of."** PBS = differentiation; drift can mimic it.
2. **Within-network distribution, NOT network-vs-genome enrichment** (that test is null — and that's fine, never claimed it).
3. **Combined dataset for Fig 3**, not SGDP-only (SGDP genome-wide is null — sample size, n≈18 AFR).
4. **n = 47 Melanesians** — biggest limitation; widens CIs, doesn't bias point estimates.

---

## If your mind goes blank — fallback phrases

- "Great question — I haven't tested that directly; my instinct is ___, but I'd check first."
- "That's on the list for the publication version."
- "Say more about what you'd expect to see?" (turns it into a conversation)
- "Let me write that down." *(then do it)*

**"I don't know," said calmly, is scientific maturity — not a miss. Bluffing is the only real mistake.**

---

*You are the world expert on this specific analysis. Nobody in the room has spent more time on these genes, these populations, these tests than you.*
