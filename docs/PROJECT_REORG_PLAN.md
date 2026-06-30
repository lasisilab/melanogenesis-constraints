# Project Reorganization Plan — Melanesia Pigmentation Genomics

**Status:** Draft for review (Yemko + Tina). Written 2026-06-23.
**Scope:** How to restructure the current `melanogenesis-constraints` repo into a
four-repo structure that supports three papers, incorporates 81 new Melanesian WGS
individuals, and adds haplotype-based scans (XP-EHH), ARGs, and Chromopainter.

> This plan is intentionally pre-otter: it fixes the *architecture* (repos, data
> handling, shared-code contract). Scientific scope per paper will be refined once
> the otter transcript is folded in. Open questions for that pass are flagged at the
> end.

---

## 1. Decisions locked

- **Topology:** one shared **core** repo + **three paper** repos (4 total).
- **Data:** raw + merged WGS/VCFs live on **Great Lakes only**. Git tracks code,
  small derived per-gene / per-locus tables, and figures. No VCFs in git.
- **This pass:** written plan only — no files moved or repos created yet.

---

## 2. Target repo structure

All under the `lasisilab` GitHub org.

```
lasisilab/
├── melanesia-popgen-core      # NEW — shared foundation (the durable asset)
├── melanogenesis-constraints  # EXISTING repo, repurposed → Paper 3
├── melanesia-prs              # NEW — Paper 1
└── melanesia-introgression    # NEW — Paper 2
```

Repo names are placeholders — adjust to lab convention.

### 2.1 `melanesia-popgen-core` — the shared foundation

Everything all three papers stand on. The single source of truth for the callset,
populations, and selection-scan machinery.

- **Callset pipeline** (migrated from current `analysis/cluster/00–17`): extract,
  SGDP liftover, the full `1KG ∪ HGDP ∪ SGDP` merge, sample/population lists,
  per-pop VCF construction. This is where the **81 Melanesian WGS** individuals get
  wired in.
- **Selection-scan engine:** PBS (existing) **plus the new XP-EHH** and any
  haplotype-based scans. Phasing step (Eagle/SHAPEIT) lives here too, since XP-EHH,
  ARGs, and Chromopainter all need phased haplotypes.
- **ARG / Chromopainter infrastructure:** ARG inference (e.g. tsinfer/tsdate, Relate,
  ARG-Needle) and Chromopainter/fineSTRUCTURE runs — shared because Papers 2 and 3
  both consume them.
- **Shared annotations:** gene regions/coordinates, the melanogenesis network gene
  set, and the **GWAS-catalog pigmentation loci** (the through-line across all three
  papers).
- **Shared Python package** (`core/` installable via `pip install -e .`): plotting
  palette, gene-list utilities, population label parsing, I/O helpers — so papers
  don't re-implement or drift.
- **Versioned derived-table releases** (see §4): the small CSVs papers actually pull.

### 2.2 `melanogenesis-constraints` (existing) → Paper 3: constraints landscape

The repo name already matches this paper, and the Quarto/Pages site + PEQG poster are
its natural content. Keep the repo, its URL, history, and collaborators; **extract the
pipeline out to core** and leave the analysis behind.

- Stays: `phase0_/phase1_/phase2_*` constraint + network analysis, the Quarto site
  (`pages/`, `_quarto.yml`, GitHub Pages workflow), poster artifacts, constraint data
  (LOEUF/PhyloP/GTEx).
- Leaves (moves to core): `analysis/cluster/*`, the merged-callset machinery, scan
  engine, shared annotations.
- Consumes core's derived tables (π, PBS, XP-EHH, ARG/Chromopainter summaries per
  gene) by pinning a core release version.

### 2.3 `melanesia-prs` — Paper 1: PRS

- PRS construction from your population genotypes + GWAS-catalog pigmentation
  associations; portability / transferability analysis across populations.
- Lightest coupling to core: needs the genotypes (callset) and the curated GWAS loci,
  not the haplotype scans. Pulls core's genotype manifest + GWAS-loci table.

### 2.4 `melanesia-introgression` — Paper 2: archaic introgression

- Archaic introgression in Melanesians, focused on pigmentation loci.
- Heaviest coupling to core's haplotype layer: phased data, ARGs, Chromopainter,
  archaic reference panels (Altai/Vindija Neanderthal, Denisova). Introgression
  callers (e.g. Sprime, IBDmix, ArchaicSeeker, or ARG-based) live here; the phased
  callset + ARGs they consume come from core.

---

## 3. How the current repo maps over (migration)

The existing `melanogenesis-constraints` becomes **Paper 3**, and **core is carved out
of it.** Concretely:

| Current location | Destination |
|---|---|
| `analysis/cluster/00–17`, `SGDP_PBS_CHECKLIST.md`, `push_to_greatlakes.sh` | **core** (pipeline) |
| `analysis/cluster/08_compute_pbs*`, `07_compute_pi*`, genome-wide scans | **core** (scan engine) |
| `data/gene_regions.bed`, network gene set, `gwas_*` tables | **core** (shared annotations) |
| `analysis/phase0/1/2_*` constraint + network scripts | stay in **constraints (Paper 3)** |
| `pages/`, `_quarto.yml`, `.github/workflows/publish.yml`, poster files | stay in **Paper 3** |
| `data/network_constraint_*`, LOEUF/PhyloP/GTEx tables | stay in **Paper 3** |
| `output/pbs_per_gene.csv`, `pi_per_gene.csv` | produced by **core**, vendored into Paper 3 |
| `SUBSTRUCTURE_PLAN.md` / `SUBSTRUCTURE_INSTRUCTIONS.md` | **core** (it's pipeline work) |

**History.** Two options for core:

- *Simple (recommended to start):* create `melanesia-popgen-core` fresh, copy the
  pipeline in with a `MIGRATED_FROM.md` noting the source commit SHA in
  `melanogenesis-constraints`. Fast, clean, loses per-file history in core.
- *History-preserving:* use `git filter-repo` to extract `analysis/cluster/` (and the
  shared-annotation paths) into the new core repo with their commits intact. More
  work; do this only if commit-level provenance of the pipeline matters.

**Repo size.** Current repo is ~410 MB (200 MB `data/`, 105 MB `.git`). When carving
core out, this is the moment to (a) move VCF-scale data out of git entirely (§4) and
(b) consider a one-time history prune if large blobs were ever committed
(`git filter-repo --strip-blobs-bigger-than 10M`) — flag to Tina first since it
rewrites history.

---

## 4. Data handling & the core ↔ paper contract

**Raw / merged WGS stays on Great Lakes.** Each repo's `.gitignore` excludes
`*.vcf.gz`, `*.bcf`, `*.bgz`, BAM/CRAM, and large intermediates. Git holds the code
that produces them and the small tables that summarize them.

**The contract between core and the papers** is *versioned derived tables*, not a live
filesystem dependency:

1. Core runs the pipeline on the cluster and emits small per-gene / per-locus tables
   (π, PBS, XP-EHH, ARG/Chromopainter summaries, GWAS-loci annotations).
2. Core tags a release (`core-v0.1`, …) with those tables + a manifest
   (sample counts per population, callset build params, cluster paths, git SHA).
3. Each paper **vendors a snapshot** of the tables it needs into its own `data/`,
   with a `PROVENANCE.md` recording the core release/SHA it pinned. Reproducible,
   and each paper is self-contained at submission time.
4. Shared *code* (palette, utils) is consumed via `pip install` of core's package
   (pin a version), not copy-paste.

This keeps papers lightweight and independently reviewable while guaranteeing they all
trace back to one callset build.

---

## 5. Where the new methods land

All pipeline-level and therefore **in core**, consumed by whichever papers need them:

- **81 Melanesian WGS** → fold into the merge step (extends the SGDP/HGDP Melanesian
  panel already described in `SUBSTRUCTURE_PLAN.md`). Re-baselines downstream scans;
  archive prior outputs and record per-population *n*.
- **Phasing** (Eagle/SHAPEIT) → new core step; prerequisite for everything haplotype-based.
- **XP-EHH** (replacing PBS as the headline scan) → core scan engine, e.g. `selscan`
  + normalization. PBS can stay as a complementary frequency-based scan.
- **ARGs** (tsinfer/tsdate, Relate, or ARG-Needle) → core; outputs feed Paper 2
  (introgression dating/local ancestry) and Paper 3 (selection landscape).
- **Chromopainter / fineSTRUCTURE** → core; ancestry "painting" of pigmentation /
  GWAS-catalog loci, shared by Papers 2 and 3.

Migration of PBS→XP-EHH is additive: keep PBS columns, append XP-EHH, re-baseline
together when the 81 WGS samples enter.

---

## 6. Suggested sequence (when you're ready to execute)

1. Create `melanesia-popgen-core`; move the pipeline in; wire `.gitignore` for VCFs.
2. Add the 81 WGS Melanesians to the merge; phase; re-baseline PBS; add XP-EHH.
3. Cut `core-v0.1` with refreshed derived tables + manifest.
4. Slim `melanogenesis-constraints` to Paper 3 (drop pipeline; vendor core-v0.1
   tables; keep site + poster).
5. Scaffold `melanesia-prs` and `melanesia-introgression` from a common paper-repo
   template (data/, analysis/, figures/, PROVENANCE.md, pinned core dependency).
6. Layer ARG + Chromopainter into core; release `core-v0.2`; papers bump their pin.

---

## 7. Open questions (resolve with the otter transcript)

1. Exact author/lead split and **timeline ordering** of the three papers — which is
   the next submission target?
2. Final repo names + whether the lab prefers a `melanesia-*` prefix or per-paper
   short names.
3. For Paper 2: which **archaic references** and which **introgression caller(s)**
   (Sprime / IBDmix / ARG-based)?
4. For Paper 1: PRS source weights — GWAS-catalog summary stats only, or also
   externally fit weights? Which target populations for the transferability test?
5. ARG tool choice (Relate vs tsinfer/tsdate vs ARG-Needle) — driven by sample size
   and whether you need branch-length dating.
6. Do we want history-preserving extraction of the pipeline into core (§3), and a
   one-time large-blob prune of the existing repo?
7. Confirm the 81 WGS samples' consent/data-use terms before any are referenced in a
   public repo (even in sample lists).
