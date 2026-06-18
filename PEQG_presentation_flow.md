# PEQG 2026 — Poster Walkthrough Flow & Talking Points

*The order that best showcases the strength and novelty of this work without over- or under-selling. Genome-wide / SGDP material stays OFF the poster and lives in the Q&A backstop. Personal prep — not for the repo.*

---

## The throughline: credit / claim / hold

One discipline runs through the whole talk. Get this right and nothing else can be over- or under-sold:

| Move | What it applies to | How it sounds |
|---|---|---|
| **Credit** the foundation | constraint↔pleiotropy; tissue-specificity↔selection in humans; regulatory-variation route | "It's well established that…" / "building on…" / "this echoes…" |
| **Claim** the contribution | the African↔Melanesian lineage contrast within this pathway | "What we find is…" / "we're not aware of a prior demonstration that…" |
| **Hold** as hypothesis | the regulatory/missense mechanism; the genome-wide framing | "our read is…" / "consistent with…" / "we test this directly in future work" |

Credit generously, claim narrowly, hold the mechanism lightly. That balance is what reads as a careful scientist.

---

## Walkthrough flow (~2.5–3 min)

### 1. Open with the tension, not a result *(~20 sec)*
- "Constraint and adaptation are usually discussed as if they act on the same genes — the genes you can't break are the genes that don't change."
- "I tested whether that's actually true, in a pathway where I can measure both: the melanogenesis network, 129 genes."
- *Why:* frames a question, not a foregone conclusion. Avoids sounding confirmatory.

### 2. The system — set the stage with Figure 1 *(~25 sec)*
- "These genes span a wide tissue-specificity gradient — broadly-expressed regulators like MITF and KIT at one end, skin-restricted enzymes like TYR and TYRP1 at the other."
- "That gradient is the axis everything else sits on."
- *Hold the line:* Figure 1 is context. Don't claim it shows constraint or selection — it's the scaffold.

### 3. Constraint strand (Figure 2) — the validated baseline *(~30 sec)*
- "First, constraint. It's well established that broadly-expressed, highly-connected genes are more constrained — and we confirm that here: breadth and connectivity each independently predict LoF constraint, and both hold in a joint model (VIF ≈ 1.16)."
- "I'm showing this *because it's expected* — it's the baseline. The real question is whether that same axis predicts where populations adapt."
- *Credit + frame:* explicitly call this the known pattern. This is your setup, not your headline — saying so is what makes you sound calibrated.

### 4. Adaptation strand (Figure 3) — the contribution *(~40 sec)*
- "Using PBS, I asked where African and Melanesian lineages show the most population-specific differentiation across these genes."
- "The striking result: the two populations sit at *opposite ends* of the tissue-specificity axis. African differentiation concentrates on broadly-expressed genes (β = −6.89, perm p = 0.011); Melanesian on tissue-restricted ones (β = +5.28, perm p = 0.048)."
- "European and South Asian signals are non-significant — a lineage-specificity control."
- *Language discipline:* "differentiation," "consistent with selection," never "proof." **The lineage contrast is the novel observation — let it land here.**

### 5. The decoupling — synthesis *(~25 sec)*
- "And LoF constraint doesn't distinguish those targets — the African- and Melanesian-associated genes don't differ in LOEUF (p = 0.97)."
- "So in this pathway, constraint and adaptation are tracking different things."
- "Our read — and this echoes work in regulatory evolution — is that pleiotropy constrains the *type* of variation that can fix, not *which* genes diverge between populations."
- *Hold as hypothesis:* "our read," "echoes work in regulatory evolution" — credits the idea and flags it as interpretation, not a measured mechanism.

### 6. Significance — what's actually new *(~20 sec)*
- "What's new isn't that tissue specificity tracks selection — that's known. It's that two dark-skinned populations draw on *opposite* ends of that axis, within one pathway."
- "Constraint sets the boundaries of the playing field; it doesn't pick the players."
- *Optional hook (only for a sharp room):* "And the African direction runs opposite to the usual genome-wide expectation that pleiotropy inhibits adaptation — which is the part we most want to explain."

---

## Anchor phrases (memorize a few)

- "Constraint sets the boundaries of the playing field; it doesn't pick the players."
- "Two dark-skinned populations draw on opposite ends of the tissue-specificity axis." ← **your safest single claim; survives any pushback**
- "Pleiotropy constrains the *type* of variation that can fix, not *which* genes diverge between populations." (deliver as interpretation)
- "We recover the expected direction in Melanesians, and a departure in Africans." (internal-consistency point)
- "Differentiation consistent with selection — not evidence of it."

---

## Q&A backstop — hold these, do NOT put them on the poster

*Deploy only if asked. Genome-wide / SGDP, mechanism, and limitations live here.*

- **"Is the network under selection genome-wide?"** → "My claim is the within-network *distribution*, not network-vs-genome enrichment. Individual effectors rank high genome-wide, but the pathway as a whole isn't a hotspot."

- **"Does the African result hold genome-wide?"** *(the one you keep off the poster)* → "Within-network and in the combined dataset, yes. The SGDP-only genome-wide scan is null for Africans, but that's an n ≈ 18 power issue — the genome-wide combined validation is the next step."

- **"Isn't this just the known tissue-specificity/selection pattern?"** → "For Melanesians, yes — we recover the expected direction, which is reassuring. The African arm is the departure, and the lineage contrast is the new part." *(Credit the adaptive-eQTL work as the prior finding you build on.)*

- **"PBS isn't selection."** → "Right — it's differentiation, consistent with selection; demography can mimic it. That's why I lean on the controls: EUR/SAS null, opposite AFR/MEL signs, permutation p-values."

- **"Why no haplotype-based methods?"** → "Frequency-based only here. Haplotype scans (iHS, XP-EHH) need phased, whole-chromosome data my targeted ±10 kb extraction doesn't support. They'd be more demography-robust and better at pinpointing recent sweeps — the natural next assay."

- **"How does adaptation act on constrained genes?"** → "Hypothesis: regulatory or missense variation that LoF metrics don't capture. We test it directly with gnomAD missense Z and eQTL overlap — future work."

- **"Why opposite signs for African vs. Melanesian τ?"** → "It's the headline. One interpretation: African selection acts on regulatory variation in pleiotropic genes; Melanesian on tissue-restricted effectors where variant effects are localized. Statistical pattern, not a mechanistic claim yet."

- **Sample size / admixture** → "n = 47 Melanesians is the main limitation — widens CIs, but the Hudson FST estimator is unbiased under unequal n, so it doesn't shift point estimates. SGDP includes some admixture; admixture-aware methods for the publication version."

- **EUR/SAS outgroup caveat** → "EUR/SAS PBS use the African–Melanesian reference structure, so they're valid lineage-specificity controls, not perfectly parallel scans."

---

## What to claim vs. what to credit (quick reference)

- **Established — credit it:** breadth & connectivity → constraint (Mähler 2017); tissue specificity tracks population-specific selection in humans, incl. with PBS (adaptive-eQTL work, PMC8756519); adaptation via regulatory variation (cis-regulatory literature; Hollis 2014).
- **Yours — claim it:** the African↔Melanesian *opposite-direction* contrast within the melanogenesis pathway; constraint-aware integration in a phenotype-anchored network.
- **Contested — engage it:** whether selection favors network hubs or periphery (Luisi 2015 vs Kim 2007) — your ρ ≈ 0 on network position adds data to an open question.

---

*Prep for PEQG 2026. Re-verify any stat against the latest CSVs before presenting.*
