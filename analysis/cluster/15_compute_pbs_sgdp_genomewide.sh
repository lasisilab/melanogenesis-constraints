#!/usr/bin/env bash
#SBATCH --job-name=sgdp_gw_pbs
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=8:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/sgdp_gw_pbs_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/sgdp_gw_pbs_%A_%a.err
#SBATCH --array=1-22%3
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 15_compute_pbs_sgdp_genomewide.sh
#
# Submits 15_compute_pbs_sgdp_genomewide.py for one chromosome (SLURM array).
# Computes per-gene PBS for all protein-coding genes using the SGDP-only
# per-population VCFs from step 14.
#
# Submit after 14_filter_sgdp_genomewide.sh has completed all 22 chromosomes:
#   sbatch 15_compute_pbs_sgdp_genomewide.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib
module load python3.10-anaconda

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
CHR=${SLURM_ARRAY_TASK_ID}

echo "[$(date)] Starting SGDP-only genome-wide PBS for chr${CHR}..."

python "${BASE}/analysis/cluster/15_compute_pbs_sgdp_genomewide.py" \
    --chr "${CHR}" \
    --base "${BASE}"

echo "[$(date)] SGDP chr${CHR} PBS complete."
