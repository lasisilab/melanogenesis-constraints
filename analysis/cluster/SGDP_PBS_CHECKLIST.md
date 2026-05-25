# SGDP-only genome-wide PBS — runbook & test plan

Keep this open next to your terminal while you submit. Each step has a
green-light condition you can grep for before moving on.

`BASE` = `/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints`

---

## Resource footprint (after the "don't hog cores" pass)

| Script | cpus | mem | time | array concurrency | peak cores in use |
|---|---|---|---|---|---|
| 14 (filter + liftover + split) | 2 | 16 G | 12 h | `1-22%2` | **4** |
| 15 (PBS compute) | 1 | 8 G | 8 h | `1-22%3` | **3** |

So at most ~4 cores at a time across all your jobs. Comfortable for a shared
account.

---

## Pre-flight checks (5 min, one-time)

```bash
ssh greatlakes
cd /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
module load Bioinformatics bcftools/1.21
```

- [ ] **bcftools has +liftover**
  ```bash
  bcftools +liftover --version 2>&1 | head -5
  # expect: bcftools version >= 1.17 and "About: Liftover ..."
  ```

- [ ] **Reference files exist** (already there from previous Melanesian liftover)
  ```bash
  ls -la data/reference/hg19.fa data/reference/hg38.fa data/hg19ToHg38.over.chain.gz
  # all three should print sizes (hg19.fa ~ 3 GB, hg38.fa ~ 3 GB, chain ~ 1 MB)
  ```

- [ ] **SGDP VCF naming**
  ```bash
  ls -la /nfs/turbo/lsa-tlasisi1/tlasisi/reference-genomes/sgdp/
  ```
  Look for one of:
  - per-chrom: `sgdp.chr{N}.vcf.gz`, `Simons.chr{N}.vcf.gz`, `chr{N}.vcf.gz`
  - single-multichrom: `Simons.vcf.gz`, `cteam_extended.v4.PS2_phase.public.vcf.gz`

  If the file naming doesn't match any of the patterns hard-coded in
  [14_filter_sgdp_genomewide.sh](14_filter_sgdp_genomewide.sh) lines 69–86,
  add yours there before submitting.

- [ ] **Logs directory exists**
  ```bash
  mkdir -p logs
  ```

---

## Step 1 — Build sample lists  (1 min, run on login node)

- [ ] **Discover the VCF if you don't know which one to point at**
  ```bash
  python analysis/cluster/13_make_sgdp_only_sample_lists.py
  # script prints candidate VCFs and exits — pick one
  ```

- [ ] **Run for real, passing --vcf**
  ```bash
  python analysis/cluster/13_make_sgdp_only_sample_lists.py \
      --vcf /nfs/turbo/lsa-tlasisi1/tlasisi/reference-genomes/sgdp/<picked-file>
  ```

- [ ] **Green-light check**: counts should look roughly like
  ```
  african    : n ≈  18–25
  melanesian : n =   17                ← same 17 you used before
  eastasian  : n ≈  35–45
  southasian : n ≈  35–45
  european   : n ≈  25–35
  excluded   : n ≈  60–80 (NorthAfr, Siberian, NativeAm, etc.)
  ```
  ```bash
  wc -l data/sgdp_only/samples_sgdp_*.txt
  ```

- [ ] **If African comes back as 0**: the SGDP sample-ID parser didn't recognize
  your VCF's naming. Inspect the audit TSV:
  ```bash
  head -50 data/sgdp_only/sgdp_population_assignments.tsv
  awk -F'\t' '$3=="EXCLUDED"' data/sgdp_only/sgdp_population_assignments.tsv | cut -f2 | sort -u
  ```
  Add missing population names to [13_make_sgdp_only_sample_lists.py](13_make_sgdp_only_sample_lists.py) lines 34–65 and rerun.

---

## Step 2 — Smoke test on chr22  (15–25 min)

**Do this before the full submission.** chr22 is the smallest autosome — a
clean run end-to-end here means the path/format assumptions are right.

- [ ] **Submit chr22 only**
  ```bash
  sbatch --array=22 analysis/cluster/14_filter_sgdp_genomewide.sh
  squeue -u $USER
  ```

- [ ] **While it runs, watch the live log**
  ```bash
  tail -f logs/sgdp_gw_filter_*_22.out
  ```
  Stages you should see (in order):
  ```
  [...] Using per-chrom SGDP VCF: ...        OR  Extracting chr22 from ...
  [...] chr22: applying quality filters in hg19 space...
  [...] chr22: NNNNNN variants after hg19 filtering.
  [...] chr22: lifting hg19 → hg38...
  [...] chr22: NNNNNN variants lifted; NN rejected.
    Writing african... african chr22: NN variants, NN samples
    Writing melanesian... melanesian chr22: NN variants, 17 samples   ← ✓
    Writing eastasian...
    Writing southasian...
    Writing european...
  [...] chr22 complete.
  ```

- [ ] **Green-light**: all 5 per-pop VCFs present
  ```bash
  ls -la vcf/sgdp_genomewide/*.chr22.vcf.gz
  # expect 5 files (african, eastasian, european, melanesian, southasian)
  bcftools query -l vcf/sgdp_genomewide/melanesian.chr22.vcf.gz | wc -l   # should be 17
  ```

- [ ] **PBS for chr22**
  ```bash
  sbatch --array=22 analysis/cluster/15_compute_pbs_sgdp_genomewide.sh
  tail -f logs/sgdp_gw_pbs_*_22.out
  ```

- [ ] **Sanity-check the PBS output**
  ```bash
  head output/pbs_sgdp_genomewide_chr22.csv
  awk -F, 'NR>1 && $13>0' output/pbs_sgdp_genomewide_chr22.csv | wc -l
  # expect: ~50–200 genes with non-zero pbs1_african out of ~500 chr22 genes
  ```
  PBS values should be small positive floats, mostly < 0.1, occasional 0.2–0.5
  outliers. If everything is 0 → an FST estimator went wrong; if everything is
  NaN → the VCFs are empty for those gene regions (could be liftover dropouts).

**🟢 If chr22 looks right, proceed to step 3.**

---

## Step 3 — Full submission  (~12–18 h wall-clock)

- [ ] **Submit step 14 for all 22 chrs**
  ```bash
  sbatch analysis/cluster/14_filter_sgdp_genomewide.sh
  squeue -u $USER
  ```
  Concurrency `%2` means only 2 jobs at a time → ~11 batches.
  Each chrom 30 min–3 h depending on size (chr1 / chr2 slowest).

- [ ] **Watch for liftover dropouts** (occasional rejected sites are fine;
  >5% rejected means a reference-version mismatch)
  ```bash
  grep "rejected" logs/sgdp_gw_filter_*.out | head -22
  ```

- [ ] **After step 14 finishes for all 22**, submit step 15
  ```bash
  # confirm all 22 per-pop VCFs exist
  for c in $(seq 1 22); do
      for p in african melanesian eastasian southasian european; do
          [[ ! -f "vcf/sgdp_genomewide/${p}.chr${c}.vcf.gz" ]] \
              && echo "MISSING: ${p}.chr${c}"
      done
  done
  # ...should print nothing

  sbatch analysis/cluster/15_compute_pbs_sgdp_genomewide.sh
  ```

- [ ] **Green-light**: 22 PBS CSVs land in `output/`
  ```bash
  ls output/pbs_sgdp_genomewide_chr*.csv | wc -l   # should be 22
  ```

---

## Step 4 — Compare HGDP vs SGDP  (run locally, ~1 min)

- [ ] **Pull SGDP results back to your laptop**
  ```bash
  scp greatlakes:/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/output/pbs_sgdp_genomewide_chr*.csv \
      ~/GitHub/melanogenesis-constraints/output/
  ```

- [ ] **Run the comparison**
  ```bash
  cd ~/GitHub/melanogenesis-constraints
  python analysis/cluster/16_compare_pbs_hgdp_vs_sgdp.py --base .
  ```

- [ ] **Look at the deliverables**
  - `output/figure_pbs_hgdp_vs_sgdp_scatter.{png,pdf}` — 2×2 scatter
  - `output/pbs_hgdp_vs_sgdp_summary.txt` — Spearman ρ + top-percentile Jaccard
  - `output/pbs_hgdp_vs_sgdp_per_gene.csv` — merged table for downstream work

**🎯 Headline question for the poster**: does **ρ_net** on PBS-3 (Melanesian
target) come back ≥ 0.5? If yes → strong replication claim. If not → dataset
effects, worth reporting honestly as a caveat.

---

## Troubleshooting cheatsheet

| Symptom | Likely cause | Fix |
|---|---|---|
| `13` reports 0 African samples | SGDP pop name missing from dict | Inspect audit TSV; add to `AFRICAN` set in [13_make_sgdp_only_sample_lists.py](13_make_sgdp_only_sample_lists.py#L34) |
| `14` aborts with "could not find SGDP VCF" | File naming not matched | Add your pattern to the loop at [14_filter_sgdp_genomewide.sh:69-77](14_filter_sgdp_genomewide.sh#L69-L77) |
| `14` >50% rejected in liftover | hg19.fa is wrong build (b37 vs hg19) or chain file mismatch | `samtools faidx data/reference/hg19.fa` should show "chr1..chr22"; if it shows "1..22" the FASTA is GRCh37 and won't match the UCSC chain |
| `14` OOM at liftover sort | Memory too tight for chr1/chr2 sort | Bump `--mem=24G` on chr1+chr2 only: `sbatch --array=1,2 --mem=24G analysis/cluster/14_filter_sgdp_genomewide.sh` |
| `15` very slow (>10 h on small chrom) | bcftools subprocess per gene region | Expected — ~19k genes × ~3 s each. Not a bug. |
| All PBS values NaN | Per-pop VCFs empty after liftover | Re-check step 14 — likely a reference mismatch |
| `16` reports `ρ_net = NaN` | Network gene names don't match between datasets | Check `merged.is_network.sum()` — should be ~100+ |

---

## Order of operations summary

```
13 (login node, ~1 min)
   ↓
[smoke test: sbatch --array=22 14, --array=22 15]
   ↓
14 (sbatch, all 22, ~12 h)
   ↓
15 (sbatch, all 22, ~4 h)
   ↓
[scp results back]
   ↓
16 (local, ~1 min)
```
