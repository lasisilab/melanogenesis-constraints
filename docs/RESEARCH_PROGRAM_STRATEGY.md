# Research Program Strategy — Pigmentation Genomics (Lasisi Lab)

**Author:** Yemko Pryor · **Drafted:** 2026-06-23
**Sources:** the three planning-meeting transcripts (Genetic Research Planning, Scientific
Research Strategy, Vitamin D Study Analysis) + prior Cowork sessions ("First paper
analysis plan", "Research Pipeline Milestones"). Companion to
[`PROJECT_REORG_PLAN.md`](PROJECT_REORG_PLAN.md).

> This document does three things: (1) states the goals of each output, (2) lays out
> the overall strategy, sequencing, and shared infrastructure, and (3) connects the
> whole program to your larger Life Research goals. It corrects a couple of points in
> the earlier reorg plan now that the meeting details are in (see §2 note).

---

## 1. The program in one paragraph

You and Tina are running one data-and-method program that spins off **four writeups**
from two cohorts (the South African SAPPHIRE cohort; Melanesian WGS + public genomes)
and a shared methods/compute backbone. The intellectual spine is constraint →
selection → genotype-to-phenotype, studied specifically in African-descent and
Oceanian populations that are underrepresented in genomics. The near-term forcing
function is Tina's move (~3 weeks out), so the immediate job is to lock scope,
sequence, and data access before she's in packing mode.

---

## 2. The four outputs and their goals

These are ordered by maturity, not necessarily by submission order.

**0. Vitamin D / UV paper (SAPPHIRE) — Tina lead, Yemko co-author.**
South African longitudinal cohort (two groups: Xhosa and Cape "Coloured"/SAC).
Core result: serum vitamin D crashes in winter regardless of sun exposure (UV too
weak at that latitude to drive synthesis); no other longitudinal signal. *Decision
from the transcripts: this is split cleanly from the genetics and written up first,
ASAP, with Tina as first author.* Your involvement is supportive, not lead.

**1. Paper 1 — SAPPHIRE PRS (Yemko lead; possibly a technical note).**
*Goal:* using the SAPPHIRE genetic data (Axiom array + imputation), test whether
pigmentation polygenic scores from GWAS-catalog associations (a) predict individual
melanin index in rank order and (b) explain the group-level pigmentation difference
(Xhosa vs SAC) via pigmentation-increasing alleles. Methylation may be folded in as a
correspondence check against the measured phenotypes.
*Key constraints already worked out:* core vitamin D biology is off-array, so a vitamin
D PRS needs imputation; African pigmentation loci (MFSD12, DDB1) are off-array too.
The African-representation strategy is keep a powered anchor + frequency-filter to
African-segregating loci + supplement with African-ancestry GWAS (Crawford 2017,
Martin 2017 KhoeSan). Imputation reference panel work is scoped (Khoe-San Genome
Project, SABR, SAHGP, Zulu — each its own EGA access request).
*Decision: Melanesians are NOT in this paper* — pooling array data with WGS is hard to
defend to reviewers, and the SAPPHIRE PRS stands cleaner on its own.

**2. Paper 2 — Selection & introgression in Melanesians ("the big paper", Yemko lead).**
*Goal:* characterize the evolutionary history of pigmentation (and genome-wide) in
Melanesians using the new WGS (~47–82 individuals, final n TBD from François) plus
publicly available Melanesian and comparison genomes — more samples than anything
currently public. Move *away* from allele-frequency selection scans (PBS) toward
haplotype- and tree-based approaches (XP-EHH, ROH/homozygosity), plus archaic
introgression at pigmentation loci and ChromoPainter ancestry "painting".
*Prerequisites surfaced in the meetings:* you can't just use standard population
groupings — you need a demography/relatedness scaffold first (TreeMix → trees →
possibly ARGs) to decide who to compare, and Gideon (trees/phylogenetics) is the
collaborator for that. This is also where the "southern route / convergent vs IBD"
question about dark-skin retention lives.

**3. Paper 3 — Melanogenesis constraints & pleiotropy (the current repo, matured).**
*Reframed goal (from the strategy meeting):* show how **pleiotropy constrains genetic
variation**, and how reading constraint in conversation with pleiotropy gives a better
picture of the **evolutionary landscape around pigmentation phenotypes**. This is the
PEQG poster grown into its conceptual form: network position × multi-axis constraint
(LOEUF/PhyloP/GTEx) × population selection. It is the capstone that ties the program's
ideas together and maps most directly onto your central research goal (§5).

> **Correction to the earlier reorg plan.** That plan implied Paper 1 (PRS) would draw
> on "your populations + GWAS catalog" generically and that the 81 Melanesians might
> feed the PRS. The meetings resolve this: **Paper 1 PRS = South African SAPPHIRE
> cohort only (array + imputation), Melanesians excluded.** The 81 Melanesian WGS feed
> **Paper 2 (selection/introgression)**, not the PRS. The four-repo topology still
> holds; only the paper↔data mapping is refined.

---

## 3. Shared infrastructure (the backbone all four lean on)

Tina's "make a repo for every technical thing" idea (Genetic Research Planning meeting)
is the same instinct as the shared-core recommendation — just expressed as **modular
method repos** rather than one monolith. Reconcile them like this:

- **Modular technical repos** (each a reusable, teachable unit, pulled into papers):
  `gwas-catalog-pull`, `imputation-pipeline`, `selection-scans` (PBS + XP-EHH + ROH),
  `chromopainter-ancestry`, `demography-trees` (TreeMix/ARG). These double as teaching
  materials and as evidence of general skills — useful for Tina's grant and your own
  independence narrative.
- **A genome-data inventory + store** (one place, download-once): expand the
  "genome list" Tina started into a tracked inventory of what's available across 1KG,
  HGDP, SGDP, dbGaP, EGA — with access status and use terms. Raw genomes live on the
  cluster, never in git.
- **The existing `melanogenesis-constraints` callset pipeline** becomes the concrete
  home of `selection-scans` + the merged-callset build (per `PROJECT_REORG_PLAN.md`).

**Compute decision (open):** Great Lakes HPC vs buying a ~$5k Mac Studio. Lean:
benchmark first on HPC (flexible mem/GPU/CPU), record %-budget per run, then decide
whether a local box is worth it. Tina is fine writing the hardware into a grant either
way; she wants the specs/benchmarks before committing. PCA/PBS are cheap; **phasing**
is the long pole for the haplotype methods.

---

## 4. Strategy, sequencing & dependencies

The critical path is **data access → demography scaffold → everything else**, because
Paper 2's comparisons and Paper 1's imputation both stall without samples in hand.

1. **Unblock data (now).** Survey + request access to the Melanesian WGS and the
   African reference/imputation panels (EGA/dbGaP). These have long lead times.
2. **Build the demography scaffold (Paper 2 gating).** TreeMix/trees on available
   genomes to decide population groupings before running scans; review the literature
   on Oceanian/African population movement to motivate the groupings; take a single
   phylogeny to Gideon for feedback on whether to escalate to ARGs.
3. **Parallel track — Paper 1 PRS.** Largely independent of Paper 2; gated only on the
   African imputation panel. Can advance while the demography work runs.
4. **Paper 3 constraints** matures from the existing repo + re-baselined scans as the
   WGS and XP-EHH come online; it's the synthesis, so it lands later.
5. **Vitamin D paper** is on Tina's plate independently and first.

Sequencing principle: don't let the conceptually-final Paper 3 block the empirically-
foundational Paper 2; and keep Paper 1 moving on its own track so there's a near-term
output while the big selection paper's infrastructure is built.

---

## 5. Connecting this to your Life Research goals

Your All of Us research statement is the clearest articulation of your north star:

> *"By learning how genetic constraint and the epigenome prime and shape the variation
> and traits that natural selection acts on, I aim to ask and better answer questions
> about how stress and inequality intersect with genetic variation to shape human
> health"* — pursued rigorously and alongside underrepresented communities.

That sentence is a research program, and the four outputs are its load-bearing pieces.
Read as pillars:

- **Pillar A — Constraint: what variation is even possible.** What bounds the genetic
  landscape (LoF intolerance, pleiotropy, network position, conservation). → **Paper 3.**
- **Pillar B — Selection: what acts on that landscape, where and how.** Haplotype and
  tree-based history of pigmentation in understudied populations. → **Paper 2.**
- **Pillar C — Genotype→phenotype in real, underrepresented cohorts.** How well
  current (Euro-biased) genetic knowledge transports; where prediction and measurement
  break down. → **Paper 1 (PRS)** + the **vitamin D / erythema** work.
- **Pillar D — Epigenome, stress & inequality → expression.** The longer-horizon arc:
  how environments of stress and inequality intersect with genetic variation to shape
  health. → the methylation thread in Paper 1, and the future All of Us / grant work.
- **Throughline — equity & rigor as method.** Underrepresentation, measurement bias in
  darker skin, honest data over confirmatory stories. This is the "so what for human
  health" that bridges your basic-science questions to NIH-style framing, and it runs
  through all four.

How the pillars compose: **A** defines the board, **B** shows how the game has been
played on it in specific populations, **C** tests whether our tools read the board
correctly for people who've been left out, and **D** is where the program is heading —
constraint and the epigenome as the layer through which stress and inequality reach the
genome's expressed output. Paper 3 is the engine, Paper 2 is the empirical proving
ground, Paper 1 is the translational/equity wedge and your nearest-term first-author
output.

This framing also serves the career layer: the modular method repos + a coherent
constraint→selection→equity story are exactly what an independent-investigator
narrative (and Tina's K-01-adjacent training-grant pitch, with Gideon as a methods
mentor) are built from.

**The apex this all serves.** Pillars A–D and these four outputs are the *foundation*,
not the destination. The longer arc — the evolutionary and health histories of African
Americans (17th c → present), reconnecting living people to ancestries severed by the
slave trade, and making the biological impact of white supremacy concrete — lives in
[`LIFE_RESEARCH_FRAMEWORK.md`](LIFE_RESEARCH_FRAMEWORK.md). That document spells out how
each project here is a deliberate scaffold (P2's ancestry toolkit, P1's prediction/
representation methods, P3's landscape frame, the measurement and epigenome layers) for
those apex questions, so they stay on the map.

---

## 6. Immediate next steps (next ~2–3 weeks, before Tina's move)

From the meetings, lightly organized:

- **Yemko — "first thing": get the genomes.** Inventory what's available (1KG/HGDP/
  SGDP + dbGaP/EGA), confirm the Melanesian WGS count with François, and start the
  access requests (longest lead time first). Expand Tina's genome-list into the tracked
  inventory.
- **Tina — vitamin D writeup (ASAP, first author);** "deal with the GWAS catalog
  thing"; then switch to prioritizing teaching/TIG.
- **Both — milestones meeting next week:** set target dates per paper and a rough
  timeline while there's still runway before the move.
- **Reach out to Gideon** about trees/demography/ARGs as a collaboration; attend
  Nicole's lab meetings.
- **Compute:** pull HPC specs + benchmark a representative run (phasing especially) to
  inform the Mac-Studio-vs-cluster decision and any grant hardware line.
- **Architecture:** when ready to execute, stand up the modular method repos +
  genome-inventory repo and extract the pipeline per `PROJECT_REORG_PLAN.md`.

---

## 7. Open questions / decisions to confirm

1. **Paper submission order** — Paper 1 (PRS) as the near-term first-author output, or
   push Paper 2 once samples land?
2. **Final Melanesian n** (François: 47 vs 54 vs 82?) and the data-use terms before any
   sample IDs touch a public repo.
3. **Population-grouping strategy for Paper 2** — how far to go beyond standard labels;
   TreeMix-only vs full ARG (Gideon's call).
4. **Vitamin D PRS anchor trade-off** (Revez European power vs African representation) —
   already scoped in the "First paper analysis plan" session; confirm with Tina.
5. **Compute** — HPC vs local box, pending benchmarks.
6. **Where the "technical repos" live** — lab org, and how thin/standalone each is.

> To sharpen §5 against your actual Life Research notes, connect your Obsidian vault to
> this session (I currently only see the project folder) and I'll align the pillars and
> language to how you've already framed them there.
