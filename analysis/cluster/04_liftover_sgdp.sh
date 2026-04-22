#!/usr/bin/env bash
#SBATCH --job-name=sgdp_liftover
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/liftover_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/liftover_%A_%a.err
#SBATCH --array=1-22%6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 04_liftover_sgdp.sh
#
# Lifts SGDP Melanesian VCFs from hg19 → hg38 using bcftools +liftover,
# then re-normalizes alleles against the hg38 reference.
#
# Why this is needed:
#   SGDP (Reich Lab) was called on GRCh37/hg19.
#   gnomAD HGDP+1KGP was called on GRCh38/hg38.
#   They cannot be merged until they share the same reference coordinates.
#
# Prereqs:
#   module load Bioinformatics bcftools/1.21 htslib
#   hg38 reference FASTA at FASTA_HG38
#   hg19 reference FASTA at FASTA_HG19 (needed by bcftools +liftover)
#   Chain file downloaded to $BASE/data/ (done in this script if missing)
#
# Submit after 02_extract_vcfs.sh completes:
#   sbatch 04_liftover_sgdp.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
RAW="${BASE}/vcf/raw"
LIFTED="${BASE}/vcf/lifted"
DATA="${BASE}/data"
CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p "${LIFTED}"

# ─── Reference files ──────────────────────────────────────────────────────────
FASTA_HG38="${BASE}/data/reference/hg38.fa"

# hg19 reference FASTA — needed by bcftools +liftover as the source reference.
# Download if not present (smaller than hg38, ~3 GB compressed / ~28 GB uncompressed).
FASTA_HG19="${BASE}/data/reference/hg19.fa"
if [[ ! -f "${FASTA_HG19}" ]]; then
    echo "[$(date)] Downloading hg19 reference FASTA (~3 GB)..."
    wget -q -O "${FASTA_HG19}.gz" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    gunzip "${FASTA_HG19}.gz"
    samtools faidx "${FASTA_HG19}"
    echo "[$(date)] hg19 FASTA ready → ${FASTA_HG19}"
fi

# Chain file: hg19 → hg38
CHAIN="${DATA}/hg19ToHg38.over.chain.gz"
if [[ ! -f "${CHAIN}" ]]; then
    echo "[$(date)] Downloading hg19→hg38 chain file..."
    wget -q -O "${CHAIN}" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    echo "[$(date)] Chain file saved → ${CHAIN}"
fi

# ─── Input: SGDP Melanesian (hg19), extracted in 02_extract_vcfs.sh ──────────
IN_VCF="${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz"

if [[ ! -f "${IN_VCF}" ]]; then
    echo "No SGDP VCF for chr${CHR} — skipping."
    exit 0
fi

# ─── Step 1: bcftools +liftover hg19 → hg38 ──────────────────────────────────
# bcftools +liftover is available as a built-in plugin (bcftools >= 1.17).
# It handles strand flips, REF/ALT swaps, and reports rejected variants.
LIFTED_VCF="${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_raw.vcf.gz"
REJECTED="${LIFTED}/sgdp_melanesian.chr${CHR}.rejected.vcf.gz"

echo "[$(date)] Lifting chr${CHR} hg19 → hg38 with bcftools +liftover..."
bcftools +liftover \
    --output-type z \
    --threads 2 \
    --output "${LIFTED_VCF}" \
    "${IN_VCF}" \
    -- \
    --src-fasta-ref "${FASTA_HG19}" \
    --fasta-ref "${FASTA_HG38}" \
    --chain "${CHAIN}" \
    --reject "${REJECTED}" \
    --reject-type z

tabix -p vcf "${LIFTED_VCF}"

if [[ -f "${REJECTED}" ]]; then
    N_DROPPED=$(bcftools view -H "${REJECTED}" 2>/dev/null | wc -l || echo 0)
    echo "  Variants dropped (unmapped): ${N_DROPPED}"
fi

# ─── Step 2: Sort (liftover may scramble order) ───────────────────────────────
echo "[$(date)] Sorting..."
SORTED="${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz"
bcftools sort \
    --output-type z \
    --threads 2 \
    --output "${SORTED}" \
    "${LIFTED_VCF}"

tabix -p vcf "${SORTED}"

# ─── Step 3: Re-normalize alleles against hg38 reference ─────────────────────
echo "[$(date)] Normalizing against hg38..."
bcftools norm \
    --fasta-ref "${FASTA_HG38}" \
    --check-ref w \
    --multiallelics -any \
    --threads 2 \
    --output-type z \
    --output "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz" \
    "${SORTED}"

tabix -p vcf "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz"

# Clean up intermediates
rm -f "${LIFTED_VCF}" "${LIFTED_VCF}.tbi" \
      "${SORTED}" "${SORTED}.tbi"

N_LIFTED=$(bcftools view -H "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz" | wc -l)
echo "[$(date)] chr${CHR} done. Variants after liftover + normalization: ${N_LIFTED}"
