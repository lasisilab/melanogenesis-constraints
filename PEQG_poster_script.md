# PEQG 2026 — Poster Walkthrough Script (Final)

*Verbatim, ~3 min spoken at the poster. Your revised flow, merged with two edits: a brief credit to prior regulatory-evolution work, and a sharper EUR/SAS line that won't be misheard as the network-enrichment claim. Bracketed cues mark figure moves and eye contact. Q&A backstop follows — reactive, not part of the walkthrough.*

---

## The walkthrough

**[Open — eye contact]**

Thanks for stopping by. This project started with a simple question: when we talk about evolutionary constraint and adaptation, are they acting on the same genes?

The usual assumption is yes. Genes that are highly constrained are often treated as genes that don't change, while genes that adapt are thought to sit at the more flexible end of the spectrum. I wanted to test that directly in a system where both processes can be measured — so I focused on the melanogenesis network, about 129 genes involved in pigmentation.

**[Gesture to Figure 1]**

What's useful about this pathway is that it spans a huge range of tissue specificity. On one end are broadly expressed regulators like MITF, KIT, and EDNRB, which function across many tissues. On the other are pigment-producing enzymes like TYR and TYRP1 that are largely restricted to melanocytes. That gradient in tissue specificity becomes the organizing axis for everything else in the poster.

**[Move to Figure 2]**

The first question was whether that axis predicts evolutionary constraint. And it does. Broadly expressed genes are significantly more constrained than tissue-restricted genes, and genes with greater network connectivity are more constrained as well — and both effects remain significant when modeled together. In other words, genes with more pleiotropic roles appear less tolerant of loss-of-function variation. That's exactly what we'd expect, so it serves as a useful baseline.

The more interesting question is whether the same factors predict where populations actually diverge.

**[Move to Figure 3]**

To test that, I used the Population Branch Statistic, or PBS, to identify genes showing elevated population-specific differentiation in African and Melanesian lineages. What emerged was a striking contrast: the two populations sit at opposite ends of the tissue-specificity axis. In Africans, differentiation is concentrated among broadly expressed genes. In Melanesians, it's concentrated among tissue-restricted genes. European and South Asian comparisons are non-significant, which helps emphasize that this pattern is lineage-specific rather than something generic to these genes across all lineages.

**[Gesture between Figures 2 and 3]**

The key point is that loss-of-function constraint doesn't separate these two sets of genes — the African- and Melanesian-associated targets are essentially indistinguishable in LOEUF. So within this pathway, constraint and adaptation are not tracking the same signal. Our interpretation, which echoes work in regulatory evolution, is that pleiotropy constrains the kinds of variants that can persist, but not necessarily which genes become targets of population-specific divergence.

**[Wrap — step back, eye contact]**

The contribution here isn't that tissue specificity matters for evolution — that's already well established. What's new is that two dark-skinned populations appear to draw on opposite ends of the same tissue-specificity spectrum within a single pigmentation pathway. Constraint helps define the evolutionary landscape, but it doesn't determine where adaptation occurs. These patterns are consistent with selection rather than direct proof of it, but together they suggest a more nuanced relationship between pleiotropy, constraint, and the evolution of pigmentation. I'd be happy to talk through any part of it.

---

## Delivery notes

- **Pace:** ~510 words ≈ 3 minutes at a calm pace; it may run closer to 3:15, so time it once. If it spills over, the Figure 1 sentence is the easiest trim.
- **Optional stickier closer:** swap *"Constraint helps define the evolutionary landscape, but it doesn't determine where adaptation occurs"* for *"Constraint sets the boundaries of the playing field; it doesn't pick the players."* — more memorable as a parting line if you want one.
- **Your safety line** if you freeze: *"Two dark-skinned populations are drawing on opposite ends of the tissue-specificity axis, within one pathway."* That sentence alone carries the contribution.
- **Numbers** are on the poster — point to the β's and p-values if asked rather than reciting them, which keeps the spoken version clean.

---

## Q&A backstop — reactive, not part of the script

*Genome-wide / SGDP, mechanism, and limitations live here. (Full set in PEQG_presentation_flow.md.)*

- **"Isn't this the known tissue-specificity/selection pattern?"** → "For Melanesians, yes — we recover the expected direction, which is reassuring. The African arm is the departure, and the lineage contrast is the new part."
- **"Does the African result hold genome-wide?"** → "Within-network and in the combined dataset, yes. The SGDP-only genome-wide scan is null for Africans, but that's an n ≈ 18 power issue — the genome-wide combined validation is the next step." *(Off-poster; raise only here.)*
- **"PBS isn't selection."** → "Right — it's differentiation, consistent with selection; demography can mimic it. That's why I lean on the controls: EUR/SAS null, opposite AFR/MEL signs, permutation p-values."
- **"How does adaptation act on constrained genes?"** → "Hypothesis: regulatory or missense variation that LoF metrics don't capture. We test it directly with missense Z and eQTL overlap — future work."
- **"Why no haplotype methods?"** → "Frequency-based only here; haplotype scans need phased whole-chromosome data my targeted extraction doesn't support. More demography-robust, and the natural next assay."
- **Sample size:** "n = 47 Melanesians is the main limitation — widens CIs, but the Hudson estimator is unbiased under unequal n, so it doesn't shift point estimates."

---

*Prep for PEQG 2026. Re-verify any stat against the latest CSVs before presenting.*
