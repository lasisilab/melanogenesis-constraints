# PEQG 2026 — Concepts From Scratch

*Every load-bearing idea on your poster, built from zero. Plain words and a small example first; the formula only if it helps. Read it in any order — each piece stands alone.*

*There is no shame in needing this. The people who look most fluent at conferences all built this scaffolding at some point; most just don't say so out loud. Once these click, the cheat-sheet phrases stop being things you recite and start being things you mean.*

---

## How the whole poster fits together (read this first)

Strip away the jargon and your poster is two measurements about each gene, asked of 129 genes:

1. **How constrained is it?** (Can evolution afford to break it? → LOEUF)
2. **Is it under recent population-specific selection?** (Did its DNA shift in one population? → PBS)

Plus a few **gene properties** you use to explain those two things: how broadly the gene is expressed (τ, breadth), and how connected it is (KEGG, betweenness).

Everything else — FST, branch lengths, regression, permutation tests — is just **machinery** for measuring those things and for honestly checking whether the patterns are real. That's it. Hold that map in your head and the parts stop feeling like a wall.

---

# PART 1 — Population genetics from zero

## Allele frequency — the atom everything is built from

At any given spot in the genome, people can carry different letters — say some chromosomes have a **G** there and some have an **A**. Each version is an **allele**.

The **allele frequency** is just the fraction of chromosomes in a population carrying one version. If you look at 100 chromosomes and 70 have G, then freq(G) = 0.70.

Different populations can have different frequencies at the same spot. **That difference is the raw material of everything in your population-genetics section.** When a population adapts, allele frequencies shift. When populations drift apart over time, frequencies drift apart too.

## FST — how different two populations are at a spot

**In plain words:** a number from 0 to 1 that measures how *differentiated* two populations' allele frequencies are. 0 = identical mixes. 1 = completely different (one population is all G, the other all A).

**Analogy:** two jars of colored marbles. FST asks how different the color *proportions* are. Same proportions → FST ≈ 0. One jar all red, the other all blue → FST = 1.

**Tiny example:** if freq(G) is 0.50 in both populations, there's no difference → FST ≈ 0. If it's 1.0 in one and 0.0 in the other, they're maximally different → FST = 1.

**"Hudson estimator"** is just *one specific recipe* for calculating FST. Its virtue: it doesn't get fooled when you've sampled far fewer people from one population than another. You sampled 747 Africans but only 47 Melanesians, so that fairness-under-unequal-samples is exactly why you picked it. You don't need the recipe — you need the one sentence: *"Hudson's FST is unbiased when sample sizes are unequal, which matters at n = 47."*

**Why it's on your poster:** FST is the building block of PBS. You can't get to selection signals without it.

## The branch-length transform: T = −ln(1 − FST)

Here's a small problem. FST is squished between 0 and 1, and it doesn't "add up" cleanly. If population A is some distance from B, and B some distance from C, the FST values don't combine by simple addition.

The transform **T = −ln(1 − FST)** stretches FST onto a new scale — a **branch length** — that *does* add up, and that grows steadily the longer two populations have been separated.

**Analogy:** imagine measuring distances along a curved road with a tape that bunches up near the end. Hard to add segments. The transform is like straightening the road into a flat ruler, so you can add and subtract segments normally.

**Why you need it:** the next step, PBS, *adds and subtracts* these distances between three populations. You can only do that arithmetic on a straight, additive ruler — which is what T gives you. (You'll never be asked to derive this. "It linearizes FST so branch lengths are additive" is the whole answer.)

## PBS — figuring out *which* population's DNA moved

**The setup:** pick a **target** population (say African) and two **outgroups** (South Asian and Papuan/Melanesian). Three populations form a little family tree with three tips.

**The intuition (picture a triangle):** between each pair of populations you have a distance (the branch length T). You want to know: *how much did the target population change on its own private branch* — the part of the tree that belongs only to it?

The trick is simple algebra. The two sides of the triangle that touch your target both contain the target's own branch. The third side — between the two outgroups — doesn't touch the target at all. So if you **add the two target-touching sides and subtract the far side** (then halve it), everything cancels except the target's own private branch.

  PBS_target = (T_target,out1 + T_target,out2 − T_out1,out2) / 2

**Analogy:** three siblings who split from one starting point. You want to know how much *one* sibling changed on their own. You compare the pairwise differences between all three and use the algebra to isolate that one sibling's personal change — separate from changes the other two share.

**What a high PBS means:** that population's allele frequencies moved a lot *on their own branch* at that gene — a candidate signature of recent selection **in that population specifically**. (Your scans: African target with S. Asian + Papuan outgroups; Melanesian target with S. Asian + African outgroups.)

## Selection vs. drift — why you always say "consistent with"

Allele frequencies can shift for two very different reasons:

- **Selection:** an allele helped (or hurt) survival/reproduction, so it spread (or shrank).
- **Drift:** pure random luck in which individuals happened to reproduce — much stronger in small populations, bottlenecks, and founder events.

**PBS sees the *shift*, but not the *reason*.** A high PBS is *consistent with* selection, but drift can produce the same footprint. That's the entire reason you say "consistent with selection," never "evidence of selection." Saying it yourself, before anyone challenges you, makes you look careful — which is exactly the impression you want.

---

# PART 2 — The gene-property metrics

## LOEUF — can the gene tolerate being broken?

**The idea:** some mutations completely **break** a gene (knock it out) — these are **loss-of-function (LoF)** variants. LOEUF asks: across huge human datasets, do we see **fewer** broken-gene mutations than we'd **expect** if breaking the gene were harmless?

- If we **expect** ~10 knockouts by chance but **observe** 0–1, that means people who broke this gene didn't make it into the dataset — natural selection removed them. The gene is **constrained** (low LOEUF).
- If observed ≈ expected, breaking it is apparently fine. The gene is **LoF-tolerant** (high LOEUF).

**Analogy:** you expect ~10 typos in a long document by random chance. You find zero in one section. Someone is clearly proofreading that section carefully — it's "constrained." A section riddled with the expected number of typos is unguarded.

LOEUF is roughly that observed/expected ratio (gnomAD reports a conservative upper bound of it). **Low = constrained, high = tolerant.** (Counterintuitive direction, so double-check yourself: *lower* LOEUF = *more* important to keep intact.)

**The crucial limitation (your MC1R point):** LOEUF only counts gene-*breaking* mutations. It's blind to subtler **missense** changes (one amino acid swapped) and **regulatory** changes (how much the gene is turned on). So a gene can look totally LoF-tolerant yet still be under selection through those other channels. That's the loophole that makes African selection on broadly-expressed, constrained genes make sense.

## τ (tau) and tissue breadth — how specialized is the gene's job?

Genes get "expressed" (switched on) in tissues. Some are on almost everywhere (housekeeping generalists); some only in one place (specialists).

- **Tissue breadth:** the simple count — *how many* of the 54 GTEx tissues is it on in (TPM > 1)?
- **τ (tau):** a polished 0-to-1 score for the same idea. **τ = 0** → expressed evenly everywhere. **τ = 1** → expressed in essentially one tissue.

**Analogy:** an employee who pitches in across every department (broad expression, low τ) versus one who only ever works in a single room (specialist, high τ).

**Why it's the star of your poster:** this is the axis everything sits on. Pigment enzymes (TYR, TYRP1, DCT) are skin specialists (high τ); master regulators (MITF, KIT, EDNRB) are generalists (low τ). You report both breadth and τ because if the same story shows up in two different ways of measuring specificity, it's not an artifact of one definition.

## KEGG count & betweenness — how connected / pleiotropic is the gene?

**Pleiotropy** = one gene affecting many things. Two ways to measure "many things":

- **KEGG pathway count:** how many different biological pathways the gene takes part in. Many pathways = a hub with many jobs across the cell.
- **Betweenness centrality:** *within the melanogenesis network specifically*, how often the gene sits on the shortest path between two other genes — i.e. how much of a bottleneck it is.

**Analogy:** KEGG count = how many company-wide committees you sit on (cross-department importance). Betweenness = how central you are inside your own department's org chart (local importance). Related, but not the same — which is why one can stand in for the other (KEGG absorbs betweenness in your model).

---

# PART 3 — Statistics from zero

## Regression and "predict" — drawing the best trend line

**Regression** finds the relationship that best links a **predictor** (x) to an **outcome** (y). "Breadth predicts LOEUF" just means: as breadth goes up, LOEUF tends to move in a consistent direction, and we can draw that trend.

The **coefficient (β)** is the **slope** of that trend: how much y changes for each unit of x. Its **sign** tells direction — negative β for breadth → LOEUF (your result) means *more breadth, lower LOEUF*.

**Important:** "predict" here means *statistical association*, not fortune-telling and not proof of cause. You're describing a trend, not claiming breadth *causes* constraint.

## Standard deviation — the "typical spread"

Quick foundation you'll need in a second: the **standard deviation (SD)** measures how spread out a set of numbers is around their average. Numbers clustered tightly → small SD. Numbers all over the place → large SD. "One SD" is shorthand for "a typical amount of variation."

## Z-scoring — putting different predictors on the same ruler

Your predictors live on wildly different scales: breadth runs 0–54, KEGG count runs 1–100+. If you compared their raw slopes, breadth's would look tiny just because its numbers span a narrower range — an unfair comparison.

**Z-scoring** fixes this. For each value: subtract the average, divide by the SD. Now every predictor is measured in the *same* unit: "how many SDs above or below average." z = 0 is average, z = +1 is one typical step above, and so on.

**Analogy:** comparing a score of 88 in an easy class to 88 in a brutal class is meaningless until you convert both to "SDs above the class average." Z-scoring is grading on a curve so the numbers are finally comparable.

After z-scoring, your β's *are* directly comparable: β_breadth = −0.135 and β_KEGG = −0.131 means the two effects are almost the same size **per typical amount of variation.**

## Marginal vs. partial — the difference behind "joint vs. single-predictor"

- **Marginal** effect: the effect of one thing *on its own*.
- **Partial** effect: the effect of one thing *after accounting for the others*.

**Analogy:** "Do taller people weigh more?" is marginal. "Do taller people weigh more *beyond what their age already explains*?" is partial — you've held age fixed.

**On your poster, this is a deliberate design choice:**
- **Figure 2 (LOEUF)** uses a **joint** model → partial effects. Because your claim is *independence*: does breadth predict constraint *even after* accounting for connectivity, and vice versa? That's a partial-effects question.
- **Figure 3 (PBS)** uses **single-predictor** models → marginal effects. Because the question is simply "does each feature, on its own, track selection?"

So if someone asks "why joint here but single there?", the honest answer is: *"different questions — one is about independence (partial), the other about each feature's own association (marginal)."*

## VIF — are my two predictors secretly the same thing?

If two predictors carry nearly identical information, a joint model can't tell them apart and the results get unstable. **VIF (variance inflation factor)** measures that redundancy. **VIF = 1** → the predictors are independent. **VIF > 5 (or 10)** → problematic overlap.

**Analogy:** predicting income from "years of schooling" and "height in cm" — unrelated, so low VIF, no problem. Predicting it from "height in cm" and "height in inches" — literally the same variable twice, so enormous VIF, and the model breaks.

**Yours is 1.16** → breadth and KEGG barely overlap → they're two genuinely different axes of pleiotropy. This single number is the backbone of your "two independent axes" claim.

## Added-variable plots — *seeing* a partial effect

These are the Panels A–D in Figure 2. An added-variable plot shows the relationship between the outcome and one predictor **after mathematically stripping out whatever the other predictors already explain.** If the points still form a slope after that stripping, the predictor still matters on its own.

**Analogy:** you want to know if coffee affects productivity, but coffee-drinkers also sleep less. An added-variable plot is like first removing the effect of sleep from both, then seeing whether coffee *still* has a slope. If it does, coffee earns its keep independently. That's the visual proof your two predictors aren't redundant.

## p-values — "could this just be a fluke?"

A **p-value** answers one specific question: *"If there were really no relationship at all, how often would random chance alone hand me a result at least this strong?"*

- Small p (say 0.01) → "this would almost never happen by luck, so I doubt it's luck."
- p = 0.05 is the conventional cutoff for "unlikely enough to take seriously."

**Analogy:** you flip a coin 10 times and get 9 heads. The p-value is the chance a *fair* coin would give you something that lopsided. It's small — so you start to suspect the coin isn't fair. The small p doesn't *prove* the coin is rigged; it just says "luck is a poor explanation."

**Common trap to avoid saying:** a p-value is **not** "the probability the result is true." It's strictly about how often chance alone would fake it.

## Parametric vs. nonparametric — and what a permutation test actually does

**Parametric** methods assume your data follow a known shape (often a bell curve) and use a formula to get the p-value. Fast, but if the data *aren't* that shape, the formula can mislead.

Your PBS percentiles are bounded between 0 and 100 and aren't bell-shaped, so you avoided leaning on those formulas. Instead you used a **permutation test**, which is **nonparametric** (assumes no particular shape). Here's the whole thing, from scratch:

1. You measured a real result — a slope β linking a gene feature to PBS.
2. You ask: *"If there were genuinely no link, what slopes would random chance produce?"*
3. You **simulate "no link" by shuffling** the PBS values randomly across the genes — this deliberately destroys any real gene-to-feature pairing while keeping all the numbers themselves identical. You recompute the slope on the shuffled data.
4. Do that **10,000 times** → a big pile of "luck-only" slopes.
5. If your *real* slope is more extreme than almost all 10,000 shuffled ones, luck is a poor explanation → small p-value.

**Analogy:** you're dealt a suspiciously good poker hand. To judge whether it's really suspicious, you reshuffle and re-deal thousands of times and see how often pure chance gives a hand that good. If almost never, your hand was genuinely unusual.

This is why, when challenged, you can say: *"It's assumption-free — I'm comparing my result against thousands of randomized versions of my own data, not against a bell curve I'd have to justify."* It's arguably the most intuitive test on the poster once you see it this way.

## Spearman, ranks, and rank-invariance — (the one from our chat)

**Ranks** = each gene's *position* when you line them up smallest to largest (1st, 2nd, 3rd…). **Percentiles** are just ranks rescaled to 0–100.

**Spearman correlation** measures association using *only the ranks* — the order — never the raw magnitudes. It asks: "when one variable climbs in rank, does the other tend to climb (or drop)?"

**Rank-invariant** means: any transform that keeps the *order* the same leaves Spearman completely unchanged — because Spearman never looked at anything but the order. (Race analogy: whether you record finishing *times* or finishing *places*, the order of who-beat-whom is identical, and a judge who only cares about order computes the same answer either way.)

**Why this is your power move:** it proves the percentile-ranking step can't have *manufactured* your Figure 3 signal (ranking only preserves order; Spearman only uses order), and it proves the SGDP sign-flip must come from the *data being different*, not from your math.

## Mann-Whitney U — "is one group generally higher?" (rank version)

A nonparametric test (ranks again, no bell-curve assumption) for whether one group's values tend to sit higher than another's. You used it for the enrichment check: *are the network genes' PBS values higher than the genome-wide background?* Answer: no (p = 0.15 African, 0.18 Melanesian) — no enrichment, which you report honestly.

## Fisher's exact test — "are these two categories linked?" in a small table

When you have a little 2×2 count table and want to know if the two categorizations are associated — and counts are small — Fisher's exact test gives an exact answer (no large-sample approximation needed).

**Yours:** {tissue-specific vs. broad} × {African-targeted vs. Melanesian-targeted} among the selected genes. Result p = 0.046 → the selection axis really does line up with the tissue-specificity axis.

## Percentiles, and within-network vs. genome-wide

A **percentile** is your position out of 100 within some group ("90th percentile" = higher than 90% of the group). The group you compare against is everything:

- **Within-network percentile:** rank a gene against the other ~124 network genes.
- **Genome-wide percentile:** rank it against all ~19,000 genes.

**Your claim lives at the within-network level** — how selection is *distributed inside the pathway*. That's why your null genome-wide *enrichment* result (the pathway as a whole isn't a hotspot) is **not** a contradiction: you never claimed the network beats the genome, only that *within* it, tissue specificity sorts which genes carry the signal. Keeping these two levels separate in your head defuses a whole category of tough questions.

---

# PART 4 — How it all chains together

Read top to bottom, this is the whole machine:

**Allele frequencies** differ between populations → you summarize each pairwise difference with **FST** (Hudson's recipe, fair for unequal samples) → you stretch FST into additive **branch lengths (T)** → you combine three of those with the **PBS** formula to isolate how much *one* population's DNA moved on its own branch = a per-gene, per-population **selection signal** (always "consistent with," because drift can mimic it).

Separately, each gene has **properties**: how broken-able it is (**LOEUF**), how specialized (**τ, breadth**), how connected (**KEGG, betweenness**).

Then **regression** is just the tool that links properties to outcomes — pleiotropy to constraint (Figure 2, **joint/partial**, checked with **VIF** and **added-variable plots**), and tissue specificity to selection (Figure 3, **single/marginal**, with **percentile** outcomes judged by **permutation tests**). **Spearman, Mann-Whitney, and Fisher** are the honest cross-checks that keep the conclusions robust to skew, sample size, and small counts.

That's the entire poster. Six measurements and the careful machinery for relating them.

---

*You don't need to hold all of this at once. Pick the three or four concepts that scare you most, get those solid, and the rest will fall into place — they all lean on the same handful of ideas.*
