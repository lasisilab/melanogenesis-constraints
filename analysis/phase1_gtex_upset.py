"""
phase1_gtex_upset.py

Phase 1.2 supplementary: UpSet plot of GTEx tissue expression patterns for
the Raghunath melanogenesis network, with LOEUF distribution per intersection.

Reads GTEx v8 median TPM, marks each gene as expressed (TPM > 1) in each of
the 54 tissues individually, then plots tissue intersections with a boxplot of
LOEUF on top.

Output:
  data/gtex_tissue_membership.csv          — gene × tissue boolean matrix
  output/figure_phase1_gtex_upset.png/pdf  — top-N intersections (≥2 genes)
  output/figure_phase1_gtex_upset_full.png/pdf — every unique intersection
"""

import os
import gzip
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from upsetplot import UpSet

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GTEX_GCT    = os.path.join(PROJECT_DIR, 'data', 'GTEx_v8_gene_median_tpm.gct.gz')
NETWORK_CSV = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
OUT_CSV     = os.path.join(PROJECT_DIR, 'data', 'gtex_tissue_membership.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

TPM_THRESHOLD   = 1.0
MIN_SUBSET_SIZE = 2     # for the top-N plot
MAX_INTERSECTIONS = 25  # cap top-N plot

# Color palette — block the two stacked panels visually
COLOR_LOEUF = '#4878CF'    # blue: LOEUF boxplot panel
COLOR_COUNT = '#D94040'    # red:  intersection count panel (gene-n bars)


def style_panels(axes, n_intersections):
    """Apply panel-level color blocking and labels."""
    loeuf_ax = axes.get('extra1') or axes.get('LOEUF')
    if loeuf_ax is not None:
        loeuf_ax.set_ylabel('LOEUF\n(higher = less constrained)', fontsize=11)
        loeuf_ax.axhline(1.0, color='#999999', linestyle='--', lw=0.7, alpha=0.5)
        loeuf_ax.set_facecolor('#F4F8FE')   # pale blue background
        loeuf_ax.tick_params(labelsize=9)
        for spine in loeuf_ax.spines.values():
            spine.set_color(COLOR_LOEUF)

    inter_ax = axes.get('intersections')
    if inter_ax is not None:
        inter_ax.set_ylabel('Genes per\nintersection (n)', fontsize=11)
        inter_ax.set_facecolor('#FDEFEF')   # pale red background
        inter_ax.tick_params(labelsize=9)
        for bar in inter_ax.patches:
            bar.set_facecolor(COLOR_COUNT)
            bar.set_edgecolor('#7A1F1F')
            bar.set_linewidth(0.4)
        for spine in inter_ax.spines.values():
            spine.set_color(COLOR_COUNT)

    totals_ax = axes.get('totals')
    if totals_ax is not None:
        totals_ax.set_xlabel('Genes expressed (n)', fontsize=10)
        totals_ax.tick_params(labelsize=8)
        for bar in totals_ax.patches:
            bar.set_facecolor('#888888')


# ── Step 1: Load GTEx and build per-gene tissue membership ────────────────
print("Parsing GTEx GCT...")
with gzip.open(GTEX_GCT, 'rt') as f:
    f.readline(); f.readline()
    gtex_df = pd.read_csv(f, sep='\t')

TISSUE_COLS = [c for c in gtex_df.columns if c not in ('Name', 'Description')]
print(f"  {len(TISSUE_COLS)} tissues")
gtex_df['gene'] = gtex_df['Description'].str.upper()

expressed_matrix = (gtex_df[TISSUE_COLS].astype(float) > TPM_THRESHOLD)
membership = pd.concat([gtex_df[['gene']], expressed_matrix], axis=1)
membership = membership.groupby('gene')[TISSUE_COLS].any().reset_index()

# IL8 → CXCL8 alias (IL8 was renamed; recover it)
if 'CXCL8' in membership['gene'].values and 'IL8' not in membership['gene'].values:
    cxcl8 = membership[membership['gene'] == 'CXCL8'].copy()
    cxcl8['gene'] = 'IL8'
    membership = pd.concat([membership, cxcl8], ignore_index=True)
    print("  Aliased IL8 → CXCL8")

# ── Step 2: Merge with network LOEUF ──────────────────────────────────────
network = pd.read_csv(NETWORK_CSV)
network['gene'] = network['gene'].str.upper()
df = network.merge(membership, on='gene', how='inner')
df = df.dropna(subset=['LOEUF'])
print(f"  {len(df)} network genes with LOEUF + GTEx tissue data "
      f"({len(network) - len(df)} missing)")

df.to_csv(OUT_CSV, index=False)
print(f"Saved → {OUT_CSV}")

upset_df = df.set_index(TISSUE_COLS)[['LOEUF', 'gene', 'functional_category']]
n_tissues = upset_df.index.to_frame().sum(axis=1)
upset_df = upset_df[n_tissues > 0]


# ===========================================================================
# Plot 1 — Top-N intersections with ≥ MIN_SUBSET_SIZE genes (clean view)
# ===========================================================================
counts = (upset_df.index.to_frame().reset_index(drop=True)
          .groupby(TISSUE_COLS).size().sort_values(ascending=False))
keep = counts[counts >= MIN_SUBSET_SIZE].head(MAX_INTERSECTIONS).index
top_df = upset_df[upset_df.index.isin(keep)]
n_keep = len(keep)
print(f"\nTop-N plot: {n_keep} intersections (≥{MIN_SUBSET_SIZE} genes), "
      f"{len(top_df)} genes; "
      f"{(counts == 1).sum()} singletons hidden")

fig = plt.figure(figsize=(14, 18))
upset = UpSet(top_df,
              subset_size='count', show_counts=True,
              sort_by='cardinality', sort_categories_by='cardinality',
              min_subset_size=MIN_SUBSET_SIZE,
              element_size=None,
              intersection_plot_elements=3)   # smaller count panel
upset.add_catplot(value='LOEUF', kind='box', color=COLOR_LOEUF,
                  elements=8)                 # taller LOEUF panel for breathing room
axes = upset.plot(fig=fig)
style_panels(axes, n_keep)

# Legend explaining the two color panels
fig.legend(handles=[
    Patch(facecolor=COLOR_LOEUF, edgecolor='gray', label='LOEUF (per intersection)'),
    Patch(facecolor=COLOR_COUNT, edgecolor='gray', label='Genes per intersection (n)'),
    Patch(facecolor='#888888',   edgecolor='gray', label='Genes expressed per tissue (totals)'),
], loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.95,
   edgecolor='gray', fontsize=10)

fig.suptitle(
    f'Tissue expression patterns and LOEUF — top {n_keep} intersections '
    f'(≥{MIN_SUBSET_SIZE} genes), {len(TISSUE_COLS)} GTEx tissues, '
    f'{len(top_df)} of {len(df)} genes shown',
    fontsize=13, fontweight='bold', y=0.995)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'figure_phase1_gtex_upset.{ext}'),
                dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  Saved figure_phase1_gtex_upset.png/pdf")


# ===========================================================================
# Plot 2 — Full UpSet: every unique intersection (no min size filter)
# ===========================================================================
n_unique = upset_df.index.unique().size
print(f"\nFull plot: {n_unique} unique intersections, {len(upset_df)} genes")

# With ~70+ intersections, scale figure width with intersection count.
fig_full = plt.figure(figsize=(max(20, n_unique * 0.42), 22))
upset_full = UpSet(upset_df,
                   subset_size='count', show_counts=True,
                   sort_by='cardinality', sort_categories_by='cardinality',
                   min_subset_size=1,
                   element_size=None,
                   intersection_plot_elements=3)
upset_full.add_catplot(value='LOEUF', kind='box', color=COLOR_LOEUF, elements=8)
axes_full = upset_full.plot(fig=fig_full)
style_panels(axes_full, n_unique)

fig_full.legend(handles=[
    Patch(facecolor=COLOR_LOEUF, edgecolor='gray', label='LOEUF (per intersection)'),
    Patch(facecolor=COLOR_COUNT, edgecolor='gray', label='Genes per intersection (n)'),
    Patch(facecolor='#888888',   edgecolor='gray', label='Genes expressed per tissue (totals)'),
], loc='upper right', bbox_to_anchor=(0.99, 0.99), framealpha=0.95,
   edgecolor='gray', fontsize=10)

fig_full.suptitle(
    f'Tissue expression patterns and LOEUF — all {n_unique} unique intersections, '
    f'{len(TISSUE_COLS)} GTEx tissues, {len(upset_df)} genes',
    fontsize=14, fontweight='bold', y=0.995)

for ext in ('png', 'pdf'):
    fig_full.savefig(os.path.join(OUT_DIR, f'figure_phase1_gtex_upset_full.{ext}'),
                     dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig_full)
print("  Saved figure_phase1_gtex_upset_full.png/pdf")

print("\nDone!")
