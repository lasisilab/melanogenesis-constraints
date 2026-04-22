#!/usr/bin/env bash
#SBATCH --job-name=melano_pi
#SBATCH --account=tlasisi0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ypryor@umich.edu
#SBATCH --output=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/pi_%j.out
#SBATCH --error=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/logs/pi_%j.err

# 07_compute_pi.sh
#
# Computes per-gene nucleotide diversity (π) for all 5 populations using
# the concatenated per-population VCFs from 06_concat_vcfs.sh.
#
# π = (1 / L) × Σ_i  2 · p_i · (1 − p_i)
#   where p_i = ALT allele frequency at site i,
#         L   = number of base pairs in the gene region (body ± 10 kb).
#
# Inputs:
#   data/gene_regions.bed             — chr, start, end, gene (128 regions)
#   vcf/final/{pop}.vcf.gz            — one per population (5 files)
#
# Output:
#   output/pi_per_gene.csv
#
# Submit after 06 completes:
#   sbatch 07_compute_pi.sh
# Or via dependency:
#   sbatch --dependency=afterok:<06_job_id> 07_compute_pi.sh

set -euo pipefail

module load Bioinformatics bcftools/1.21 htslib

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints

echo "[$(date)] Starting π computation"
echo "  Base: ${BASE}"

python3 "${BASE}/analysis/cluster/07_compute_pi.py" --base "${BASE}"

echo "[$(date)] Done. Output: ${BASE}/output/pi_per_gene.csv"
