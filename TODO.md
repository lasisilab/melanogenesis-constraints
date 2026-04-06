# TODO — Melanogenesis Constraints (PEQG 2026)

A lightweight task tracker for coordination between Yemko and Tina.
For full analysis details and phase breakdowns, see [PEQG_POSTER_PLAN.md](PEQG_POSTER_PLAN.md).

---

## Decisions needed

- [ ] **Population dataset:** Use gnomAD HGDP+1KGP joint callset, or download 1KGP and HGDP separately?
- [ ] **African populations:** YRI only, or pool multiple West African groups (YRI, LWK, ESN, GWD, MSL)?
- [ ] **Melanesian populations:** Bougainville, Papuan, or both?
- [ ] **PBS outgroup:** CEU, CHB, or both?
- [ ] **Gene region definition:** Gene body only, or ± flanking region? If flanking, how much (e.g., ±10kb)?
- [ ] **iHS:** Include on poster, or mark as "in progress"?

---

## Up next (ready to start)

- [ ] PhyloP conservation scores — download + merge to `network_constraint_data.csv` (Phase 1.1)
- [ ] GTEx tissue expression breadth — download + compute per gene (Phase 1.2)
- [ ] Start population data acquisition (Phase 2.1) — *blocked on decisions above*

---

## Done

- [x] LOEUF × functional category
- [x] LOEUF × disease class
- [x] Betweenness centrality × LOEUF
- [x] GWAS coding vs. regulatory constraint
- [x] LOEUF by clinical phenotype
- [x] Gene list overlap (Raghunath / Baxter / Bajpai)
