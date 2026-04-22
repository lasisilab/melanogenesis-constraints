#!/usr/bin/env bash
#SBATCH --job-name=melano_pbs
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/pbs_%j.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/pbs_%j.err

# 08_compute_pbs.sh
#
# Computes per-gene PBS (Population Branch Statistic) for 4 approved scans
# using the Hudson FST estimator (Bhatia et al. 2013 Genome Research).
#
# PBS formula:
#   T_XY = -ln(1 - FST_XY)         (transform FST to branch length)
#   PBS_X = (T_XY + T_XZ - T_YZ) / 2
#
# Approved scans (Tina Lasisi):
#   PBS-1: Target=African,    Outgroup=SouthAsian, Distant=Papuan
#   PBS-2: Target=African,    Outgroup=European,   Distant=Papuan
#   PBS-3: Target=Papuan,     Outgroup=SouthAsian, Distant=African
#   PBS-4: Target=Papuan,     Outgroup=European,   Distant=African
#
# Inputs:
#   data/gene_regions.bed             — chr, start, end, gene (128 regions)
#   vcf/final/{pop}.vcf.gz            — one per population (5 files)
#
# Output:
#   output/pbs_per_gene.csv
#
# Can run concurrently with 07 (both read from vcf/final/, neither writes to it).
# Submit after 06 completes:
#   sbatch 08_compute_pbs.sh
# Or via dependency:
#   sbatch --dependency=afterok:<06_job_id> 08_compute_pbs.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints

echo "[$(date)] Starting PBS computation"
echo "  Base: ${BASE}"

python3 "${BASE}/analysis/cluster/08_compute_pbs.py" --base "${BASE}"

echo "[$(date)] Done. Output: ${BASE}/output/pbs_per_gene.csv"
