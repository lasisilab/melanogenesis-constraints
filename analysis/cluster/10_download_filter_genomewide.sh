#!/usr/bin/env bash
#SBATCH --job-name=melano_gw_download
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/gw_download_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/gw_download_%A_%a.err
#SBATCH --array=1-22%4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 10_download_filter_genomewide.sh
#
# Streams one full chromosome from gnomAD HGDP+1KGP v3.1.2, applies quality
# filters, then writes per-population whole-chromosome VCFs.  These are the
# inputs for 11_compute_pbs_genomewide.sh.
#
# Changes vs. the targeted pipeline (02/05):
#   - No --regions-file: whole chromosome retained
#   - FILTER="PASS" enforced (gnomAD VQSR/RF passed sites only)
#   - Per-sample GQ < 20 set to missing, then sites >10% missing dropped
#   - Melanesian = HGDP only (PapuanHighlands, PapuanSepik, Bougainville)
#     SGDP samples excluded to avoid pre-calling asymmetry at genome-wide scale
#
# Storage estimate per chromosome:
#   Temp filtered all-sample VCF:  ~8–15 GB (deleted after pop subsets written)
#   Per-pop VCFs (5 pops):         ~1–3 GB each → ~5–15 GB per chromosome
#   Total persistent (22 chrs):    ~110–330 GB
#
# Prerequisites:
#   module load Bioinformatics bcftools/1.21 htslib
#   data/samples_*.txt written by 03_make_sample_lists.py
#
# Submit after 09_fetch_gene_annotation.py has run:
#   sbatch 10_download_filter_genomewide.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
DATA="${BASE}/data"
GW_DIR="${BASE}/vcf/genomewide"
LOGS="${BASE}/logs"
CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p "${GW_DIR}" "${LOGS}"

# ── Sample lists ──────────────────────────────────────────────────────────────
AFRICAN="${DATA}/samples_african.txt"
MELANESIAN="${DATA}/samples_melanesian_hgdp.txt"   # HGDP only — no SGDP
EASTASIAN="${DATA}/samples_eastasian.txt"
SOUTHASIAN="${DATA}/samples_southasian.txt"
EUROPEAN="${DATA}/samples_european.txt"

# Union of all study samples (for initial filter step)
ALL_KEEP="${DATA}/samples_all_keep_hgdp_only.txt"
if [[ ! -f "${ALL_KEEP}" ]]; then
    cat "${AFRICAN}" "${MELANESIAN}" "${EASTASIAN}" "${SOUTHASIAN}" "${EUROPEAN}" \
        | sort -u > "${ALL_KEEP}"
    echo "  Created ${ALL_KEEP}"
fi

# ── gnomAD source URL ─────────────────────────────────────────────────────────
GNOMAD_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.hgdp_tgp.chr${CHR}.vcf.bgz"

# ── Temp file (deleted at end) ────────────────────────────────────────────────
FILTERED="${GW_DIR}/filtered.chr${CHR}.vcf.gz"

echo "[$(date)] chr${CHR}: streaming gnomAD + applying filters..."

# Step 1: Stream → keep PASS sites → keep study samples → biallelic SNPs
#         → set GQ<20 to missing → drop >10% missing sites
bcftools view \
    --apply-filters PASS \
    --samples-file "${ALL_KEEP}" \
    --force-samples \
    --min-ac 1 \
    --types snps \
    -m2 -M2 \
    --output-type u \
    --threads 2 \
    "${GNOMAD_URL}" \
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
    --threads 2 \
    --output "${FILTERED}" \
    -- -t AC,AN,AF

tabix -p vcf "${FILTERED}"

N_VARS=$(bcftools view -H "${FILTERED}" | wc -l)
N_SAMP=$(bcftools query -l "${FILTERED}" | wc -l)
echo "[$(date)] chr${CHR}: ${N_VARS} variants, ${N_SAMP} samples after filtering."

# ── Step 2: Write per-population subsets ─────────────────────────────────────
declare -A POPS=(
    ["african"]="${AFRICAN}"
    ["melanesian"]="${MELANESIAN}"
    ["eastasian"]="${EASTASIAN}"
    ["southasian"]="${SOUTHASIAN}"
    ["european"]="${EUROPEAN}"
)

for POP in "${!POPS[@]}"; do
    OUT="${GW_DIR}/${POP}.chr${CHR}.vcf.gz"
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
        --threads 2 \
        --output "${OUT}" \
        "${FILTERED}"

    tabix -p vcf "${OUT}"
    N=$(bcftools view -H "${OUT}" | wc -l)
    NS=$(bcftools query -l "${OUT}" | wc -l)
    echo "  ${POP} chr${CHR}: ${N} variants, ${NS} samples"
done

# ── Step 3: Delete temp filtered VCF to save space ───────────────────────────
rm -f "${FILTERED}" "${FILTERED}.tbi"
echo "[$(date)] chr${CHR}: temp filtered VCF deleted."

echo "[$(date)] chr${CHR} complete."
