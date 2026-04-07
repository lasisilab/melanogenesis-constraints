#!/usr/bin/env bash
# 00_setup_dirs.sh
# Run once on the cluster to create the working directory tree.
# Usage: bash 00_setup_dirs.sh

BASE=/nfs/turbo/lsa-tlasisi1/tlasisi/melanosome-constraints

mkdir -p "${BASE}"/{data,logs,vcf/{raw,filtered},results/{pi,pbs,ihs}}

# data/      — gene region BED file, sample lists, merged VCFs
# vcf/raw/   — per-chromosome region extracts (one per chrom × dataset)
# vcf/filtered/ — merged, biallelic-only, population-filtered VCFs
# results/   — π, PBS, iHS outputs
# logs/      — SLURM stdout/stderr

echo "Directory structure created under ${BASE}:"
find "${BASE}" -type d | sort
