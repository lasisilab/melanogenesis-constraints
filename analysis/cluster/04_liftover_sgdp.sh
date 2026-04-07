#!/usr/bin/env bash
#SBATCH --job-name=sgdp_liftover
#SBATCH --account=tlasisi1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints/logs/liftover_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints/logs/liftover_%A_%a.err
#SBATCH --array=1-22

# 04_liftover_sgdp.sh
#
# Lifts SGDP Melanesian VCFs from hg19 → hg38 using CrossMap,
# then re-normalizes alleles against the hg38 reference.
#
# Why this is needed:
#   SGDP (Reich Lab) was called on GRCh37/hg19.
#   gnomAD HGDP+1KGP was called on GRCh38/hg38.
#   They cannot be merged until they share the same reference coordinates.
#
# What happens during liftover:
#   - Most SNPs (~98%) transfer cleanly.
#   - A small fraction are dropped: variants in regions that moved,
#     split into multiple intervals, or whose strand flipped ambiguously.
#   - After liftover, bcftools norm re-aligns alleles to hg38 to correct
#     any strand issues.
#
# Prereqs:
#   module load crossmap bcftools htslib
#   hg38 reference FASTA (see FASTA_HG38 below — set to your cluster's copy)
#   Chain file downloaded to $BASE/data/ (done in this script if missing)
#
# Submit after 02_extract_vcfs.sh completes:
#   sbatch 04_liftover_sgdp.sh

set -euo pipefail

module load crossmap bcftools htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints
RAW="${BASE}/vcf/raw"
LIFTED="${BASE}/vcf/lifted"
DATA="${BASE}/data"
CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p "${LIFTED}"

# ─── Reference files ──────────────────────────────────────────────────────────
# hg38 reference FASTA — update this path to your cluster's copy.
# Common locations on UM Great Lakes:
#   /nfs/turbo/lsa-tlasisi1/resources/hg38/hg38.fa
#   /scratch/reference/hg38/GRCh38_no_alt.fa
# If you don't have one, download from UCSC:
#   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
FASTA_HG38="/path/to/hg38.fa"   # <-- UPDATE THIS

# Chain file: hg19 → hg38
CHAIN="${DATA}/hg19ToHg38.over.chain.gz"
if [[ ! -f "${CHAIN}" ]]; then
    echo "Downloading hg19→hg38 chain file..."
    wget -q -O "${CHAIN}" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    echo "Chain file saved → ${CHAIN}"
fi

# ─── Input: SGDP Melanesian (hg19), extracted in 02_extract_vcfs.sh ──────────
IN_VCF="${RAW}/sgdp_melanesian.chr${CHR}.vcf.gz"

if [[ ! -f "${IN_VCF}" ]]; then
    echo "No SGDP VCF for chr${CHR} — skipping."
    exit 0
fi

# ─── Step 1: CrossMap liftover hg19 → hg38 ───────────────────────────────────
LIFTED_RAW="${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_raw.vcf"
REJECTED="${LIFTED}/sgdp_melanesian.chr${CHR}.rejected.vcf"

echo "[$(date)] Lifting chr${CHR} hg19 → hg38..."
CrossMap.py vcf \
    "${CHAIN}" \
    "${IN_VCF}" \
    "${FASTA_HG38}" \
    "${LIFTED_RAW}" \
    --chromid s           # output chrom format: "chr1" not "1"

# CrossMap writes a plain VCF and a .unmap file for failed variants
UNMAPPED="${LIFTED_RAW}.unmap"
if [[ -f "${UNMAPPED}" ]]; then
    N_DROPPED=$(grep -vc "^#" "${UNMAPPED}" 2>/dev/null || echo 0)
    echo "  Variants dropped (unmapped): ${N_DROPPED}"
fi

# ─── Step 2: Sort, compress, index ───────────────────────────────────────────
echo "[$(date)] Sorting and compressing..."
bcftools sort \
    --output-type z \
    --output "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz" \
    "${LIFTED_RAW}"

tabix -p vcf "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz"

# ─── Step 3: Re-normalize alleles against hg38 reference ─────────────────────
# This corrects any strand inconsistencies introduced by liftover.
# -m -any splits multiallelic sites; -f left-aligns indels (good practice for SNPs too).
echo "[$(date)] Normalizing against hg38..."
bcftools norm \
    --fasta-ref "${FASTA_HG38}" \
    --check-ref w \          # warn but don't fail on ref mismatches; review warnings
    --multiallelics -any \
    --output-type z \
    --output "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz" \
    "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz"

tabix -p vcf "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz"

# Clean up intermediate files
rm -f "${LIFTED_RAW}" \
      "${LIFTED_RAW}.unmap" \
      "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz" \
      "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38_sorted.vcf.gz.tbi"

N_LIFTED=$(bcftools view -H "${LIFTED}/sgdp_melanesian.chr${CHR}.hg38.vcf.gz" | wc -l)
echo "[$(date)] chr${CHR} done. Variants after liftover + normalization: ${N_LIFTED}"
