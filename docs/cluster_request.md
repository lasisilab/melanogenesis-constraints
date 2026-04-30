# Cluster Resource Request: Melanogenesis Constraint Population Analysis
**Prepared by:** Yemko Pryor  
**Date:** 2026-04-07  
**Location:** `/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/`

---

## What this analysis does

We are computing nucleotide diversity (π) and the Population Branch Statistic (PBS) across **128 melanogenesis network genes** in African and Melanesian populations, using East and South Asian populations as outgroups. This is the population genomics component of the PEQG 2026 poster.

The key design choice is **targeted extraction**: instead of downloading full genome VCFs, we use remote tabix-indexed access to stream only our 128 gene regions (gene body ± 10 kb flanking = **13.4 Mb of sequence**) directly from cloud storage. This keeps storage requirements very low.

> **Scope note:** iHS (integrated haplotype score) requires phased whole-chromosome data and a genome-wide normalization baseline, making it incompatible with targeted extraction. iHS is **out of scope for PEQG 2026** and deferred to future work. The PEQG poster deliverables are π and PBS only.

---

## Datasets

| Dataset | Source | Access | Samples in study | Notes |
|---------|--------|--------|-----------------|-------|
| gnomAD HGDP+1KGP v3.1.2 | Google Cloud (public) | HTTP stream | ~1,744 | hg38; no download needed |
| SGDP (Papuan populations) | Reich Lab (public) | HTTP stream | ~83 | hg19 → must liftover to hg38 |

### Population breakdown

| Role | Source | Populations | N |
|------|--------|-------------|---|
| **Target: African** | 1KGP | YRI, LWK, ESN, GWD, MSL | ~405 |
| **Target: African** | HGDP | Yoruba, Mandenka | ~43 |
| *Excluded* | 1KGP | ASW, ACB (admixed) | — |
| **Target: Melanesian** | HGDP | PapuanHighlands, PapuanSepik, Bougainville | ~30 |
| **Target: Melanesian** | SGDP | Papuan populations | ~15 |
| **Outgroup: South Asian** | 1KGP + HGDP | GIH, PJL, BEB, STU, ITU, Balochi, etc. | ~790 |
| **Outgroup: European** | 1KGP | CEU, TSI, FIN, GBR, IBS | ~? |
| **Outgroup: European** | HGDP | French, Sardinian, Basque, Orcadian, etc. | ~? |
| **Outgroup: East Asian** | 1KGP + HGDP | CHB, JPT, CHS, CDX, KHV, Han, Japanese | ~718 |
| **Total study samples** | | | **~2,285+** |

*Note: East Asian samples are retained in the dataset but are not used in the four primary PBS comparisons.*

### PBS scan design (approved by T. Lasisi)

| Scan | Target | Outgroup | Distant outgroup |
|------|--------|----------|-----------------|
| PBS-1 | African | South Asian | Papuan |
| PBS-2 | African | European | Papuan |
| PBS-3 | Papuan | South Asian | African |
| PBS-4 | Papuan | European | African |

---

## Storage estimate

### Analysis footprint (what we are storing)

| File category | Description | Estimated size |
|---------------|-------------|---------------|
| Gene regions BED | 128 genes ± 10kb, input to bcftools | < 1 MB |
| Sample list files | Per-population .txt files | < 1 MB |
| gnomAD metadata | HGDP+1KGP sample table | ~5 MB |
| SGDP metadata | Sample table | < 1 MB |
| **Raw extracted VCFs — gnomAD** | 128 gene regions × 4,400 samples × chr1–22 | **~300–500 MB** |
| **Raw extracted VCFs — SGDP** | 128 gene regions × 83 Melanesian samples × chr1–22 | **~20–50 MB** |
| Lifted SGDP VCFs | After hg19→hg38 CrossMap liftover | ~20–50 MB |
| Chain file (hg19→hg38) | UCSC download, needed once | ~100 MB |
| **Merged + filtered study VCFs** | ~1,744 samples, biallelic SNPs only | **~300–500 MB** |
| Per-population subset VCFs | African, Melanesian, E.Asian, S.Asian | ~200–400 MB |
| Analysis outputs (π, PBS) | Per-gene summary tables and figures | ~50–100 MB |
| SLURM logs | stdout/stderr for all jobs | ~10 MB |
| **Total** | | **~1.0–1.7 GB** |

### A note on how these sizes are calculated

The VCF storage estimates above cover **all ~1,744 individuals combined** — VCFs are multi-sample files, so all populations are stored together in one file per chromosome. The numbers are small because we are extracting only **13.4 Mb of target sequence** (128 genes ± 10 kb) out of ~3,000 Mb total genome (~0.4%). Even with 1,744 samples, a compressed multi-sample VCF for a tiny region is modest in size.

For comparison: a single whole-genome sequence is ~30–50 GB compressed. At 0.4% of the genome, the proportional expected size for one person would be ~120–200 MB — and across 1,744 people in a shared VCF format (which stores genotypes efficiently), ~300–500 MB total is reasonable.

### Reference genome FASTA (hg38)

The hg38 reference FASTA is **not currently on the cluster** and will need to be downloaded. It is required by CrossMap (liftover of SGDP) and bcftools norm (allele renormalization).

| File | Compressed | Uncompressed |
|------|-----------|-------------|
| hg38.fa.gz (UCSC) | ~3 GB download | ~30 GB on disk |

> CrossMap and bcftools both require the uncompressed (or bgzip-compressed) FASTA to be indexed. Plan for **~30 GB** for the reference alone.

### Total storage summary

| Category | Size |
|----------|------|
| All extracted + filtered VCFs (all ~1,744 samples combined) | ~1–2 GB |
| hg38 reference FASTA (uncompressed + index) | ~30 GB |
| Chain file, metadata, logs, results | ~0.5 GB |
| **Total** | **~32–33 GB** |

---

## Compute estimate

The pipeline runs as SLURM array jobs (chr1–22 in parallel). All jobs run on a single node.

### Per-job resource requirements

| Script | CPUs | Memory | Wall time / job | Concurrency |
|--------|------|--------|----------------|-------------|
| `02_extract_vcfs.sh` (stream + extract) | 1 | 8 GB | ~1–4 hr | 4 at a time |
| `04_liftover_sgdp.sh` (CrossMap) | 1 | 8 GB | ~15–30 min | 4 at a time |
| `05_merge_filter_vcfs.sh` (merge + filter) | 1 | 16 GB | ~15–30 min | 4 at a time |
| `06_concat_vcfs.sh` (concatenate per population) | 1 | 8 GB | ~15 min | single job |

### Total compute

| Stage | Jobs | CPU-hours | Wall time (sequential batches) |
|-------|------|-----------|-------------------------------|
| Extract gnomAD + SGDP | 22 array jobs (%4) | ~22–88 | ~3–9 hr |
| Liftover SGDP | 22 array jobs (%4) | ~6–11 | ~30–60 min |
| Merge + filter | 22 array jobs (%4) | ~6–11 | ~1–2 hr |
| Concatenate per population | 1 job | ~1 | ~15 min |
| π / PBS analysis | 1–2 jobs | ~5–10 | ~2–4 hr |
| **Total** | | **~40–125 CPU-hours** | **~7–16 hr wall time** |

> **Note:** The extraction step (`02_extract_vcfs.sh`) dominates wall time because it streams from remote HTTP. Actual CPU utilization is low during network I/O — the jobs are mostly waiting on data transfer, not computing. A faster network partition (if available) would reduce this substantially.

---

## Step-by-step execution plan

```
Step 1  bash 00_setup_dirs.sh          # create directory tree (~1 min)
Step 2  python 03_make_sample_lists.py  # download metadata, write sample files (~5 min)
Step 3  sbatch 02_extract_vcfs.sh       # array: extract gene regions from gnomAD + SGDP
Step 4  sbatch 04_liftover_sgdp.sh      # array: liftover SGDP hg19 → hg38
Step 5  sbatch 05_merge_filter_vcfs.sh  # array: merge datasets, filter to study samples
Step 6  sbatch 06_concat_vcfs.sh        # concatenate per-chromosome VCFs into per-population files
Step 7  [downstream] compute π and PBS per gene per population
```

Steps 3–6 depend on each other and must run in order. Submit as a dependency chain:

```bash
JOB1=$(sbatch --parsable analysis/cluster/02_extract_vcfs.sh)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 analysis/cluster/04_liftover_sgdp.sh)
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 analysis/cluster/05_merge_filter_vcfs.sh)
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 analysis/cluster/06_concat_vcfs.sh)
echo "Chain: $JOB1 → $JOB2 → $JOB3 → $JOB4"
```

---

## Key decisions still open

| Decision | Recommendation | Status |
|----------|---------------|--------|
| Gene region flanking | ±10 kb (currently set) | Confirm with Tina |
| iHS | Out of scope for PEQG 2026 — deferred to future work | **Decided** |
| hg38 reference FASTA | Downloaded and indexed at `data/reference/hg38.fa` | **Done** |

---

## What we are NOT downloading

To be explicit: we are **not** downloading full-chromosome VCFs from gnomAD (which would be ~100–200 GB per chromosome, ~2–4 TB total). The targeted extraction approach streams only the bytes covering our 128 gene regions via tabix-indexed HTTP access. The storage estimates above already reflect this.

---

*Scripts: `analysis/cluster/` in the melanogenesis-constraints repo (peqg-poster branch)*
