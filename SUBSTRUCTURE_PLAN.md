# Population Substructure — Implementation Plan

Derived from `SUBSTRUCTURE_INSTRUCTIONS.md`, grounded in the actual state of the
cluster pipeline. This plan supersedes the instruction file where they conflict;
conflicts are called out under **Corrections** below.

## Locked decisions

- **Standard dataset going forward = `1KG ∪ HGDP ∪ SGDP` for *every* population**,
  not just Melanesian. The full SGDP callset is lifted to hg38 and merged in.
- **East African source:** SGDP East Africans (Somali, Masai, Dinka, Luo,
  BantuKenya) **+ 1KG LWK**, as a blend — a natural consequence of the full merge,
  not a separate cross-dataset bolt-on.
- **PBS-1…4 will be re-baselined.** Adding SGDP to African/SAS/EUR/EAS changes the
  FST pairs feeding the approved scans, so their *values* move (better-powered FST,
  same column names). The current `output/pbs_per_gene.csv` is archived as
  `pbs_per_gene_hgdp1kg_only.csv` before the rebuild, and the re-baseline is flagged
  to Tina Lasisi (scans were marked "approved").
- **Phasing:** Phase 0 (data-foundation rebuild) → Phase 1 (substructure branches).
  The earlier "gnomAD-only first / defer East African" two-phase split is
  **superseded**: with full SGDP in the merge, East African is just another branch.

---

## Corrections to the original instructions

1. **SGDP is Melanesian-only end-to-end today, not just at the merge.**
   `02_extract_vcfs.sh` extracts only `sgdp_melanesian.chr{N}`; `04_liftover_sgdp.sh`
   lifts only that subset (`04_liftover_sgdp.sh:71-72`); `05` merges only it
   (`05_merge_filter_vcfs.sh:56`). The full merge therefore requires changes to
   **02, 04, and 05**, not a one-line addition. (This is Phase 0.)

2. **`08_compute_pbs.py` reads per-population VCFs, not sample lists.** It reads
   `vcf/final/{pop}.vcf.gz` (`08_compute_pbs.py:185`), produced by the `05` subset
   loop (`05_merge_filter_vcfs.sh:108-141`) → `06_concat_vcfs.sh`. Every new branch
   needs its own per-pop VCF built through that chain before `08` can score it.

3. **East African was internally inconsistent in the instructions** (§1 said gnomAD
   LWK, §3 said SGDP). Resolved by the full merge: `african_east = LWK ∪ SGDP East`.

4. **EUR/SAS are never PBS targets** — `pbs_per_gene.csv` has only `pbs1/2_african`
   and `pbs3/4_melanesian`; EUR/SAS appear only as outgroups. The instructions'
   `["AFR","MEL","EUR","SAS"]` branch-loop reference does not map onto PBS target
   columns; verify the real branch→column mapping before editing the regression
   (Phase 1, Step 1.4).

5. **Classification caveat (flag, don't silently accept):** instruction §2 lumps
   `Mbuti, Ju_hoan_North, Khomani_San` into `AFRICAN_WEST`. Those are Central African
   rainforest hunter-gatherers and Southern African Khoe-San — deeply diverged from
   West African farmers, not "West African." Confirm the intended grouping before
   building `african_west` (see Open items).

---

## Phase 0 — Data-foundation rebuild (full SGDP merge)

Goal: produce a study VCF where every population is `1KG ∪ HGDP ∪ SGDP`, re-baseline
the existing PBS, and archive the old output. Nothing substructure-specific yet.

### Step 0.0 — Archive + sanity checks (before touching data)

- Copy `output/pbs_per_gene.csv` → `output/pbs_per_gene_hgdp1kg_only.csv` (the
  pre-SGDP baseline for comparison).
- `bcftools query -l` the SGDP source (Simons) VCF and confirm the full sample set /
  population labels parse via `parse_population_from_sample_id`
  (`13_make_sgdp_only_sample_lists.py:83-91`).

### Step 0.1 — `02_extract_vcfs.sh`: extract the full SGDP panel

- Extend extraction from Melanesian-only to all SGDP samples in the 5-pop panel
  (the union of the `AFRICAN/MELANESIAN/EAST_ASIAN/SOUTH_ASIAN/EUROPEAN` sets in
  script 13). Output `sgdp.chr{N}.vcf.gz` (rename from `sgdp_melanesian.*`).

### Step 0.2 — `04_liftover_sgdp.sh`: lift the full SGDP callset

- Generalize input/output from `sgdp_melanesian.chr{N}` to `sgdp.chr{N}`. Liftover
  logic (`+liftover` → sort → norm) is unchanged; only file names + sample count
  grow. Expect more variants and a larger `rejected` set — log both.

### Step 0.3 — `13_make_sgdp_only_sample_lists.py`: produce SGDP per-pop lists

- This already enumerates SGDP samples per population from the VCF header. Run it
  against the full lifted SGDP VCF to emit `samples_sgdp_{african,melanesian,
  eastasian,southasian,european}.txt`. Also add the African sub-splits here now
  (West/East/South per instruction §2) so Phase 1 can consume them.
- **Print n per group; flag n < 15.** SGDP has few individuals per population.

### Step 0.4 — `05_merge_filter_vcfs.sh`: merge full SGDP, generalize per-pop subset

- Merge gnomAD + full lifted SGDP (`sgdp.chr{N}.hg38.vcf.gz`).
- **Generalize the per-pop subset loop**: today only Melanesian cats HGDP+SGDP
  (`05:121-125`). Make *every* population's sample file `cat(gnomAD_list, sgdp_list)`
  so african/SAS/EUR/EAS all gain their SGDP samples.
- `samples_all_keep.txt` (`03_make_sample_lists.py:149-151`) must include the full
  SGDP set, not just the 17 Melanesians — otherwise `05`'s `--samples-file ALL_KEEP`
  silently drops them. Update `03` (or have `05` union in the SGDP lists).

### Step 0.5 — `06_concat_vcfs.sh` + re-run PBS-1…4

- Concat per-chrom → `vcf/final/{pop}.vcf.gz` for the five base pops.
- Re-run `08_compute_pbs.py` unchanged to regenerate PBS-1…4 on the new dataset.
- Diff against the archived baseline; record the magnitude of the re-baseline for
  the writeup / Tina.

### Phase 0 done-when

`vcf/final/{pop}.vcf.gz` each contain 1KG∪HGDP∪SGDP samples; PBS-1…4 regenerated;
old baseline archived; per-pop n's recorded.

---

## Phase 1 — Substructure branches

Now branches off the SGDP-inclusive dataset. East African is no longer special.

### Step 1.0 — Verification (gating)

- Confirm gnomAD metadata actually carries Highland/Lowland Papuan labels
  (`PapuanHighlands`, `PapuanSepik`) — unverified; `03:80` prints unique pops.
- Per-group n for `african_west`, `african_east`, `melanesian_highlands`,
  `melanesian_lowlands` (each now 1KG/HGDP + SGDP where available). Flag n < 15.
- Note Highland/Lowland is effectively HGDP-only: the 17 SGDP `S_Papuan-*` have no
  highland/lowland tag, so they cannot be split and stay only in combined Melanesian.

### Step 1.1 — `03_make_sample_lists.py`: gnomAD split lists

- Keep combined `samples_african.txt` / `samples_melanesian_hgdp.txt` (backward
  compat). Add `AFRICAN_WEST_*`/`AFRICAN_EAST_*` and Highland/Lowland sets per §1;
  write the four new gnomAD lists.

### Step 1.2 — Build per-pop VCFs for the four new branches

- Add `african_west`, `african_east`, `melanesian_highlands`, `melanesian_lowlands`
  to the `05` subset loop (or a re-runnable `05b`), each with
  `cat(gnomAD_split_list, sgdp_split_list)`:
  - `african_west`  = gnomAD West ∪ SGDP West
  - `african_east`  = gnomAD LWK ∪ SGDP East
  - `melanesian_highlands` / `melanesian_lowlands` = HGDP-only subsets
- Extend `06_concat_vcfs.sh` to concat these four into `vcf/final/`.

### Step 1.3 — `08_compute_pbs.py`: add PBS-5…8

- Append to `PBS_SCANS` (`08:64-69`), leaving PBS-1…4 definitions intact:
  - `pbs5_african_west`  = (african_west, southasian, melanesian)
  - `pbs6_african_east`  = (african_east, southasian, melanesian)
  - `pbs7_mel_highlands` = (melanesian_highlands, southasian, african_west)
  - `pbs8_mel_lowlands`  = (melanesian_lowlands, southasian, african_west)
- `FST_PAIRS`, `pops_needed`, the per-gene loop, and column ordering all derive from
  `PBS_SCANS` — no other logic changes. New `fst_*`/`pbs*` columns append; verify
  `08:267-274` keeps existing columns first.

### Step 1.4 — Downstream within-network regression

- First read `analysis/phase2_pbs_within_network.py` /
  `phase2_pbs_within_network_combined.py` to find the real branch→column mapping
  (Correction #4). Add `AFR_WEST`, `AFR_EAST`, `MEL_HIGHLANDS`, `MEL_LOWLANDS` with
  the same recipe (Spearman ρ + adjusted linear model + permutation over τ,
  betweenness, KEGG count, LOEUF, tissue breadth). Append new sections to
  `output/phase2_pbs_within_network_combined.txt`; preserve existing ones.

---

## Methodological notes for the writeup

- **Re-baseline disclosure.** PBS-1…4 values change vs the published/approved run
  because all populations now include SGDP. Report old-vs-new alongside the archived
  `pbs_per_gene_hgdp1kg_only.csv`.
- **Residual cross-dataset effects.** The full merge makes every FST pair a
  gnomAD+SGDP blend on both sides (symmetric), which is far better than an SGDP-only
  branch vs gnomAD outgroups — but SGDP samples per population are few, so branches
  leaning heavily on SGDP (african_east especially) are still n-limited. Report n.
- **Liftover loss.** Track variants dropped during the full-SGDP liftover
  (`04` rejected file) — larger now than the Melanesian-only run.

## What NOT to touch

- PBS-1…4 *scan definitions* (their values re-baseline, but the code/columns stay).
- Existing columns in `pbs_per_gene.csv` (append only) and existing sections in
  `phase2_pbs_within_network_combined.txt` (append only).
- Figure-generation scripts (separate pass after results confirmed).
- Genome-wide PBS scripts (11, 15) and the SGDP-only genome-wide pipeline (14–16) —
  out of scope here; revisit after per-gene results look sensible.

## Open items to verify on cluster (blocking)

1. Full SGDP sample set + population labels parse cleanly (Step 0.0).
2. gnomAD metadata carries Highland/Lowland Papuan labels (Step 1.0).
3. Per-group n for all SGDP pops and all four Phase-1 branches; flag n < 15.
4. Intended `african_west` grouping — are Mbuti / Ju|'hoan / Khomani San really meant
   to be "West African"? (Correction #5.)
5. Branch→PBS-column mapping in the phase2 regression scripts (Correction #4).
6. Magnitude of the PBS-1…4 re-baseline vs the archived HGDP+1KG-only run.
