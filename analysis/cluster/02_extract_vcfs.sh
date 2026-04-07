#!/usr/bin/env bash
#SBATCH --job-name=melano_vcf_extract
#SBATCH --account=tlasisi1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints/logs/extract_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints/logs/extract_%A_%a.err
#SBATCH --array=1-22

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

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints
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
    --threads 4 \
    --output "${RAW}/hgdp_1kgp.chr${CHR}.vcf.gz" \
    "${GNOMAD_URL}"

tabix -p vcf "${RAW}/hgdp_1kgp.chr${CHR}.vcf.gz"
echo "[$(date)] gnomAD chr${CHR} done."

# ─── SGDP (Melanesian populations only) ───────────────────────────────────────
# Simons Genome Diversity Project — per-chromosome VCFs from Reich Lab
# NOTE: SGDP is hg19 (GRCh37). You will need to liftover to hg38 before merging
#       with gnomAD, OR liftover the gene_regions.bed to hg19 for this step.
#       Liftover step is handled in 04_liftover_merge.sh (to be written).
#
# SGDP VCFs: https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/
SGDP_URL="https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/cteam_extended.v4.PS2_phase.public.chr${CHR}.vcf.gz"
SGDP_SAMPLES="${BASE}/data/samples_sgdp_melanesian.txt"

if [[ -f "${SGDP_SAMPLES}" ]]; then
    echo "[$(date)] Extracting chr${CHR} from SGDP (Melanesian)..."
    bcftools view \
        --regions-file "${REGIONS}" \
        --samples-file "${SGDP_SAMPLES}" \
        --min-ac 1 \
        --types snps \
        --output-type z \
        --threads 4 \
        --output "${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz" \
        "${SGDP_URL}"

    tabix -p vcf "${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz"
    echo "[$(date)] SGDP chr${CHR} done."
else
    echo "WARNING: ${SGDP_SAMPLES} not found — skipping SGDP. Run 03_make_sample_lists.py first."
fi

echo "[$(date)] chr${CHR} extraction complete."
