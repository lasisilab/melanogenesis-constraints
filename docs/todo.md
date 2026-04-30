# TODO — Melanogenesis Constraints (PEQG 2026)

A lightweight task tracker for coordination between Yemko and Tina.
For full analysis details and phase breakdowns, see [plan.md](plan.md).

---

## Decisions — RESOLVED

- [x] **Population dataset:** Use gnomAD HGDP+1KGP joint callset + SGDP for Melanesian depth
- [x] **African populations:** Pool all 1KGP unambiguously continental African groups: YRI, LWK, ESN, GWD, MSL (~405 samples). Supplement with HGDP African (Yoruba, Mandenka, ~43). **Exclude** ASW (African Americans) and ACB (African Caribbeans in Barbados) — admixed.
- [x] **Melanesian populations:** HGDP Papuan + Bougainville (~37) + SGDP Papuan (~83). Use all.
- [x] **Outgroups (expanded):** Include East Asian (CHB, JPT, CHS, CDX, KHV from 1KGP + HGDP Han/Japanese) and South Asian (GIH, PJL, BEB, STU, ITU from 1KGP + HGDP South Asian) in addition to European (CEU) for PBS and FST comparisons.
- [ ] **Gene region definition:** Gene body only, or ± flanking region? If flanking, how much (e.g., ±10kb)?
- [ ] **iHS:** Include on poster, or mark as "in progress"?

## Decisions — PENDING

- [ ] **Gene region definition:** ±10kb recommended; confirm with Tina.
- [ ] **iHS:** Time permitting — `selscan` on phased 1KGP data; otherwise note as ongoing.

---

## Up next (ready to start)

- [ ] Start population data acquisition (Phase 2.1) — decisions resolved, ready to begin

---

## Done

- [x] LOEUF × functional category
- [x] LOEUF × disease class
- [x] Betweenness centrality × LOEUF
- [x] GWAS coding vs. regulatory constraint
- [x] LOEUF by clinical phenotype
- [x] Gene list overlap (Raghunath / Baxter / Bajpai)
- [x] PhyloP conservation scores — `data/phylop_scores.csv`, `data/network_constraint_phylop.csv`, `output/figure_phase1_phylop.png/pdf`
- [x] GTEx tissue expression breadth — `data/gtex_tissue_breadth.csv`, `data/network_constraint_gtex.csv`, `output/figure_phase1_gtex.png/pdf`
- [x] Quarto site with Phase 0 + Phase 1 results (`index.qmd`, `phase1.qmd`)
