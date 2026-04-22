#!/usr/bin/env bash
#SBATCH --job-name=melano_vcf_extract
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/extract_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/extract_%A_%a.err
#SBATCH --array=1-22%6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 02_extract_vcfs.sh
#
# Streams gene regions (±10kb) from the gnomAD HGDP+1KGP v3.1.2 joint callset
# and the SGDP (Melanesian populations) via remote tabix queries.
# No full-chromosome download needed — bcftools streams only the requested regions.
#
# Prereqs on cluster:
#   module load bcftools htslib
#   gene_regions.bed and sample list files in $BASE/data/
#
# Submit:
#   sbatch 02_extract_vcfs.sh

set -euo pipefail

module load bcftools htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
REGIONS="${BASE}/data/gene_regions.bed"
RAW="${BASE}/vcf/raw"
CHR=${SLURM_ARRAY_TASK_ID}

# ─── gnomAD HGDP+1KGP v3.1.2 ─────────────────────────────────────────────────
# Streams regions from Google Cloud (HTTP accessible without gsutil)
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${CHR}.vcf.bgz"

echo "[$(date)] Extracting chr${CHR} from gnomAD HGDP+1KGP..."
bcftools view \
    --regions-file "${REGIONS}" \
    --min-ac 1 \
    --types snps \
    --output-type z \
    --threads 1 \
    --output "${RAW}/hgdp_1kgp.chr${CHR}.vcf.gz" \
    "${GNOMAD_URL}"

tabix -p vcf "${RAW}/hgdp_1kgp.chr${CHR}.vcf.gz"
echo "[$(date)] gnomAD chr${CHR} done."

# ─── SGDP (Melanesian populations only) ───────────────────────────────────────
# Local copy of Simons.vcf.gz (all chromosomes, hg19, chr naming: "1" not "chr1").
# Indexed with tabix before this job runs.
# Regions extracted using gene_regions_hg19.bed (chr prefix stripped).
SGDP_VCF="${BASE}/data/sgdp/Simons.vcf.gz"
SGDP_SAMPLES="${BASE}/data/samples_melanesian_sgdp.txt"
REGIONS_HG19="${BASE}/data/gene_regions_hg19.bed"

if [[ -f "${SGDP_SAMPLES}" && -f "${SGDP_VCF}" && -f "${REGIONS_HG19}" ]]; then
    echo "[$(date)] Extracting chr${CHR} from local SGDP (Melanesian)..."
    bcftools view \
        --regions-file "${REGIONS_HG19}" \
        --samples-file "${SGDP_SAMPLES}" \
        --min-ac 1 \
        --types snps \
        --output-type z \
        --threads 1 \
        --output "${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz" \
        "${SGDP_VCF}"

    tabix -p vcf "${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz"
    echo "[$(date)] SGDP chr${CHR} done."
else
    [[ ! -f "${SGDP_VCF}" ]]      && echo "WARNING: ${SGDP_VCF} not found — skipping SGDP."
    [[ ! -f "${SGDP_SAMPLES}" ]]  && echo "WARNING: ${SGDP_SAMPLES} not found — run 03_make_sample_lists.py."
    [[ ! -f "${REGIONS_HG19}" ]]  && echo "WARNING: ${REGIONS_HG19} not found — run: sed 's/^chr//' gene_regions.bed > gene_regions_hg19.bed"
fi

echo "[$(date)] chr${CHR} extraction complete."
