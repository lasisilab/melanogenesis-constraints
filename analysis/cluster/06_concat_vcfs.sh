#!/usr/bin/env bash
#SBATCH --job-name=melano_concat
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/concat_%j.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/concat_%j.err

# 06_concat_vcfs.sh
#
# Concatenates per-chromosome filtered VCFs into one file per population.
# Run after 05_merge_filter_vcfs.sh completes.
#
# Inputs (from 05):
#   vcf/filtered/{population}.chr{1-22}.vcf.gz
#
# Outputs:
#   vcf/final/african.vcf.gz
#   vcf/final/melanesian.vcf.gz
#   vcf/final/eastasian.vcf.gz
#   vcf/final/southasian.vcf.gz
#
# Submit after 05 completes:
#   sbatch 06_concat_vcfs.sh
# Or via dependency chain:
#   sbatch --dependency=afterok:<05_job_id> 06_concat_vcfs.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
FILTERED="${BASE}/vcf/filtered"
FINAL="${BASE}/vcf/final"

mkdir -p "${FINAL}"

POPULATIONS=(african melanesian eastasian southasian european)

for POP in "${POPULATIONS[@]}"; do
    echo "[$(date)] Concatenating chromosomes for ${POP}..."

    # Build list of per-chromosome VCFs in order (chr1–22)
    CHR_VCFS=()
    for CHR in $(seq 1 22); do
        VCF="${FILTERED}/${POP}.chr${CHR}.vcf.gz"
        if [[ -f "${VCF}" ]]; then
            CHR_VCFS+=("${VCF}")
        else
            echo "  WARNING: Missing ${VCF} — skipping chr${CHR} for ${POP}"
        fi
    done

    if [[ ${#CHR_VCFS[@]} -eq 0 ]]; then
        echo "  ERROR: No VCFs found for ${POP} — skipping."
        continue
    fi

    echo "  Found ${#CHR_VCFS[@]} chromosome VCFs"

    OUT="${FINAL}/${POP}.vcf.gz"
    bcftools concat \
        --allow-overlaps \
        --output-type z \
        --output "${OUT}" \
        "${CHR_VCFS[@]}"

    tabix -p vcf "${OUT}"

    N_VARS=$(bcftools view -H "${OUT}" | wc -l)
    N_SAMPS=$(bcftools query -l "${OUT}" | wc -l)
    echo "  [$(date)] ${POP}: ${N_VARS} variants, ${N_SAMPS} samples → ${OUT}"
done

echo "[$(date)] All populations concatenated."
echo ""
echo "Output files:"
for POP in "${POPULATIONS[@]}"; do
    OUT="${FINAL}/${POP}.vcf.gz"
    if [[ -f "${OUT}" ]]; then
        N_VARS=$(bcftools view -H "${OUT}" | wc -l)
        N_SAMPS=$(bcftools query -l "${OUT}" | wc -l)
        echo "  ${POP}: ${N_VARS} variants, ${N_SAMPS} samples"
    fi
done
