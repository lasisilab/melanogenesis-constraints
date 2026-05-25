#!/usr/bin/env bash
#SBATCH --job-name=sgdp_gw_filter
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/sgdp_gw_filter_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/sgdp_gw_filter_%A_%a.err
#SBATCH --array=1-22%2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 14_filter_sgdp_genomewide.sh
#
# Genome-wide filtering and per-population splitting of the SGDP VCFs for the
# SGDP-only PBS pipeline. Mirrors 10_download_filter_genomewide.sh but reads
# from a local SGDP VCF (hg19) and lifts it to hg38 per chromosome before
# splitting.
#
# Per chromosome:
#   1. Extract the SGDP VCF chunk for this chromosome (hg19)
#   2. Apply per-sample GQ < 20 → missing, drop sites with >10% missing,
#      keep biallelic SNPs
#   3. Lift hg19 → hg38 with bcftools +liftover (reuses chain file already
#      downloaded by 04_liftover_sgdp.sh)
#   4. Write per-population subsets using lists from 13_make_sgdp_only_sample_lists.py
#
# Outputs (in ${BASE}/vcf/sgdp_genomewide/):
#   african.chr${CHR}.vcf.gz       — per-pop SGDP VCFs in hg38
#   melanesian.chr${CHR}.vcf.gz
#   eastasian.chr${CHR}.vcf.gz
#   southasian.chr${CHR}.vcf.gz
#   european.chr${CHR}.vcf.gz
#
# Prerequisites:
#   - SGDP VCF(s) in ${SGDP_DIR}, hg19
#   - data/sgdp_only/samples_sgdp_*.txt  (from 13_make_sgdp_only_sample_lists.py)
#   - data/reference/hg19.fa, data/reference/hg38.fa, data/hg19ToHg38.over.chain.gz
#     (already in place from 04_liftover_sgdp.sh)
#
# Submit after 13_make_sgdp_only_sample_lists.py has run:
#   sbatch 14_filter_sgdp_genomewide.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib samtools

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
SGDP_DIR=/nfs/turbo/lsa-tlasisi1/tlasisi/reference-genomes/sgdp
DATA="${BASE}/data"
SGDP_OUT="${BASE}/vcf/sgdp_genomewide"
LOGS="${BASE}/logs"
CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p "${SGDP_OUT}" "${LOGS}"

# ── Sample lists ──────────────────────────────────────────────────────────────
SGDP_DATA="${DATA}/sgdp_only"
AFRICAN="${SGDP_DATA}/samples_sgdp_african.txt"
MELANESIAN="${SGDP_DATA}/samples_sgdp_melanesian.txt"
EASTASIAN="${SGDP_DATA}/samples_sgdp_eastasian.txt"
SOUTHASIAN="${SGDP_DATA}/samples_sgdp_southasian.txt"
EUROPEAN="${SGDP_DATA}/samples_sgdp_european.txt"
ALL_KEEP="${SGDP_DATA}/samples_sgdp_all_keep.txt"

for f in "${AFRICAN}" "${MELANESIAN}" "${EASTASIAN}" "${SOUTHASIAN}" \
         "${EUROPEAN}" "${ALL_KEEP}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: sample list missing: ${f}" >&2
        echo "Run 13_make_sgdp_only_sample_lists.py first." >&2
        exit 1
    fi
done

# ── Reference files for liftover ──────────────────────────────────────────────
FASTA_HG19="${BASE}/data/reference/hg19.fa"
FASTA_HG38="${BASE}/data/reference/hg38.fa"
CHAIN="${DATA}/hg19ToHg38.over.chain.gz"

for f in "${FASTA_HG19}" "${FASTA_HG38}" "${CHAIN}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: reference file missing: ${f}" >&2
        echo "These should already exist from 04_liftover_sgdp.sh." >&2
        exit 1
    fi
done

# ── Locate SGDP VCF(s) ────────────────────────────────────────────────────────
# Try several layouts in priority order. The first match wins.
#   1) per-chrom files: sgdp.chr${CHR}.vcf.gz  /  Simons.chr${CHR}.vcf.gz
#   2) single multi-chrom VCF (extract this chrom via bcftools)
SGDP_CHR_VCF=""
for pat in "sgdp.chr${CHR}.vcf.gz"      "Simons.chr${CHR}.vcf.gz" \
           "sgdp_chr${CHR}.vcf.gz"      "Simons_chr${CHR}.vcf.gz" \
           "chr${CHR}.vcf.gz"; do
    if [[ -f "${SGDP_DIR}/${pat}" ]]; then
        SGDP_CHR_VCF="${SGDP_DIR}/${pat}"
        echo "[$(date)] Using per-chrom SGDP VCF: ${SGDP_CHR_VCF}"
        break
    fi
done

# Whether we extracted this chrom from a multi-chrom VCF (so we know to delete)
EXTRACTED_FROM_MULTI=0

if [[ -z "${SGDP_CHR_VCF}" ]]; then
    # Look for a single multi-chrom VCF
    SGDP_FULL_VCF=""
    for pat in Simons.vcf.gz cteam_extended.v4.PS2_phase.public.vcf.gz \
               sgdp.vcf.gz SGDP.vcf.gz; do
        if [[ -f "${SGDP_DIR}/${pat}" ]]; then
            SGDP_FULL_VCF="${SGDP_DIR}/${pat}"
            break
        fi
    done
    if [[ -z "${SGDP_FULL_VCF}" ]]; then
        echo "ERROR: could not find SGDP VCF for chr${CHR} in ${SGDP_DIR}" >&2
        echo "       Tried per-chrom patterns and common multi-chrom names." >&2
        echo "       ls ${SGDP_DIR}:" >&2
        ls -la "${SGDP_DIR}" >&2
        exit 1
    fi
    echo "[$(date)] Extracting chr${CHR} from ${SGDP_FULL_VCF}..."

    # SGDP uses hg19 chromosome naming — try both "1" and "chr1"
    CHR_TAG="${CHR}"
    if bcftools view -h "${SGDP_FULL_VCF}" 2>/dev/null | grep -q "##contig=<ID=chr${CHR},"; then
        CHR_TAG="chr${CHR}"
    fi

    SGDP_CHR_VCF="${SGDP_OUT}/sgdp_hg19.chr${CHR}.vcf.gz"
    bcftools view -r "${CHR_TAG}" \
        --samples-file "${ALL_KEEP}" \
        --force-samples \
        --min-ac 1 \
        --types snps \
        -m2 -M2 \
        --output-type z \
        --threads 1 \
        --output "${SGDP_CHR_VCF}" \
        "${SGDP_FULL_VCF}"
    tabix -p vcf "${SGDP_CHR_VCF}"
    EXTRACTED_FROM_MULTI=1
    echo "[$(date)] Extracted ${SGDP_CHR_VCF}"
fi

# ── Step 1: Filter (hg19): keep biallelic SNPs, GQ<20→missing, drop high-miss ──
FILTERED_HG19="${SGDP_OUT}/filtered_hg19.chr${CHR}.vcf.gz"
echo "[$(date)] chr${CHR}: applying quality filters in hg19 space..."

bcftools view \
    --samples-file "${ALL_KEEP}" \
    --force-samples \
    --min-ac 1 \
    --types snps \
    -m2 -M2 \
    --output-type u \
    --threads 1 \
    "${SGDP_CHR_VCF}" \
| bcftools filter \
    --set-GTs . \
    --exclude 'FORMAT/GQ < 20' \
    --output-type u \
    --threads 1 \
| bcftools filter \
    --exclude 'F_MISSING > 0.1' \
    --output-type u \
    --threads 1 \
| bcftools +fill-tags \
    --output-type z \
    --threads 1 \
    --output "${FILTERED_HG19}" \
    -- -t AC,AN,AF

tabix -p vcf "${FILTERED_HG19}"

N_VARS_HG19=$(bcftools view -H "${FILTERED_HG19}" | wc -l)
echo "[$(date)] chr${CHR}: ${N_VARS_HG19} variants after hg19 filtering."

# ── Step 2: Liftover hg19 → hg38 ──────────────────────────────────────────────
FILTERED_HG38="${SGDP_OUT}/filtered.chr${CHR}.vcf.gz"
REJECTED="${SGDP_OUT}/rejected.chr${CHR}.vcf.gz"

echo "[$(date)] chr${CHR}: lifting hg19 → hg38..."
bcftools +liftover \
    --output-type z \
    --threads 1 \
    --output "${FILTERED_HG38}" \
    "${FILTERED_HG19}" \
    -- \
    --src-fasta-ref "${FASTA_HG19}" \
    --fasta-ref "${FASTA_HG38}" \
    --chain "${CHAIN}" \
    --reject "${REJECTED}" --write-reject

# Sort post-liftover (chain remaps can disorder records) and re-index
bcftools sort \
    --output-type z \
    --output "${FILTERED_HG38}.sorted" \
    "${FILTERED_HG38}"
mv "${FILTERED_HG38}.sorted" "${FILTERED_HG38}"
tabix -f -p vcf "${FILTERED_HG38}"

N_VARS_HG38=$(bcftools view -H "${FILTERED_HG38}" | wc -l)
N_REJ=$(bcftools view -H "${REJECTED}" 2>/dev/null | wc -l || echo 0)
echo "[$(date)] chr${CHR}: ${N_VARS_HG38} variants lifted; ${N_REJ} rejected."

# ── Step 3: Write per-population subsets ─────────────────────────────────────
declare -A POPS=(
    ["african"]="${AFRICAN}"
    ["melanesian"]="${MELANESIAN}"
    ["eastasian"]="${EASTASIAN}"
    ["southasian"]="${SOUTHASIAN}"
    ["european"]="${EUROPEAN}"
)

for POP in "${!POPS[@]}"; do
    OUT="${SGDP_OUT}/${POP}.chr${CHR}.vcf.gz"
    if [[ -f "${OUT}" ]]; then
        echo "  ${POP} chr${CHR} already exists — skipping."
        continue
    fi
    echo "  Writing ${POP}..."
    bcftools view \
        --samples-file "${POPS[$POP]}" \
        --force-samples \
        --min-ac 1 \
        --output-type z \
        --threads 1 \
        --output "${OUT}" \
        "${FILTERED_HG38}"
    tabix -p vcf "${OUT}"
    N=$(bcftools view -H "${OUT}" | wc -l)
    NS=$(bcftools query -l "${OUT}" | wc -l)
    echo "  ${POP} chr${CHR}: ${N} variants, ${NS} samples"
done

# ── Step 4: Cleanup temp files ────────────────────────────────────────────────
rm -f "${FILTERED_HG19}" "${FILTERED_HG19}.tbi"
rm -f "${FILTERED_HG38}" "${FILTERED_HG38}.tbi"
[[ "${EXTRACTED_FROM_MULTI}" -eq 1 ]] && rm -f "${SGDP_CHR_VCF}" "${SGDP_CHR_VCF}.tbi"

echo "[$(date)] chr${CHR} complete."
