#!/usr/bin/env bash
#SBATCH --job-name=melano_gw_pbs
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/gw_pbs_%A_%a.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/gw_pbs_%A_%a.err
#SBATCH --array=1-22%6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu

# 11_compute_pbs_genomewide.sh
#
# Runs 11_compute_pbs_genomewide.py for one chromosome (SLURM array).
# Computes per-gene PBS for all ~19,000 protein-coding genes using the
# whole-chromosome per-population VCFs from step 10.
#
# Submit after 10_download_filter_genomewide.sh completes all 22 chromosomes:
#   sbatch 11_compute_pbs_genomewide.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib
module load python3.10-anaconda

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
CHR=${SLURM_ARRAY_TASK_ID}

echo "[$(date)] Starting genome-wide PBS for chr${CHR}..."

python "${BASE}/analysis/cluster/11_compute_pbs_genomewide.py" \
    --chr "${CHR}" \
    --base "${BASE}"

echo "[$(date)] chr${CHR} PBS complete."
