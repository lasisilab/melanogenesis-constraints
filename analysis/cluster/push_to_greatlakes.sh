#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
scp \
  analysis/cluster/13_make_sgdp_only_sample_lists.py \
  analysis/cluster/14_filter_sgdp_genomewide.sh \
  analysis/cluster/15_compute_pbs_sgdp_genomewide.py \
  analysis/cluster/15_compute_pbs_sgdp_genomewide.sh \
  analysis/cluster/16_compare_pbs_hgdp_vs_sgdp.py \
  analysis/cluster/SGDP_PBS_CHECKLIST.md \
  ypryor@greatlakes.arc-ts.umich.edu:/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints/analysis/cluster/
