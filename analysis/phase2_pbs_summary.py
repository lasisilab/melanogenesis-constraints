"""
phase2_pbs_summary.py

Cross-population PBS summary figure for the PEQG 2026 poster.

Summarises the four PBS scans into two population-averaged signals:
  mean_afr = (PBS-1 + PBS-2) / 2   African selection, averaged over outgroup choice
  mean_mel = (PBS-3 + PBS-4) / 2   Melanesian selection, averaged over outgroup choice

Panel A — 2-D scatter:
  X = mean_afr,  Y = mean_mel
  Nodes sized by rank-sum score (signal breadth across all 4 scans).
  Threshold lines at 75th percentile of each axis divide the plot into
  four quadrants (African-specific / Melanesian-specific / Both / Neither).
  Diagonal y = x marks equal signal in both populations.

Panel B — Diverging bar chart (top 20 genes by rank-sum):
  African bar extends LEFT (mean_afr), Melanesian bar extends RIGHT (mean_mel).
  Sorted by mean_afr + mean_mel.  Gene labels coloured by functional category.

Output: output/figure_phase2_pbs_summary.png / .pdf
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from scipy import stats
from adjustText import adjust_text

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_gtex.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared constants ───────────────────────────────────────────────────────
CATEGORY_ORDER = [
    'Pigment-specific',
    'Developmental/NC',
    'Generic signaling',
    'Cytokines/growth factors',
    'Apoptosis/cell death',
    'Other',
]
CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',
    'Developmental/NC':          '#E8907E',
    'Generic signaling':         '#F5C242',
    'Cytokines/growth factors':  '#4878CF',
    'Apoptosis/cell death':      '#6BAD6B',
    'Other':                     '#B0B0B0',
}

# ── Load & merge ───────────────────────────────────────────────────────────
print("Loading data...")
master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()

pbs = pd.read_csv(PBS_CSV)
pbs['gene'] = pbs['gene'].str.upper()
for col in ['pbs1_african', 'pbs2_african', 'pbs3_melanesian', 'pbs4_melanesian']:
    pbs[col] = pbs[col].clip(lower=0)

df = pbs.merge(master[['gene', 'functional_category', 'LOEUF']], on='gene', how='left')

# ── Derived summary scores ─────────────────────────────────────────────────
df['mean_afr'] = (df['pbs1_african']    + df['pbs2_african'])    / 2
df['mean_mel'] = (df['pbs3_melanesian'] + df['pbs4_melanesian']) / 2

# Rank-sum: rank within each scan (ascending → low PBS = low rank)
for col in ['pbs1_african', 'pbs2_african', 'pbs3_melanesian', 'pbs4_melanesian']:
    df[f'rank_{col}'] = df[col].rank(ascending=True, method='average')
df['rank_sum'] = (df['rank_pbs1_african'] + df['rank_pbs2_african'] +
                  df['rank_pbs3_melanesian'] + df['rank_pbs4_melanesian'])
n_genes = len(df)
df['rank_pct'] = df['rank_sum'] / (4 * n_genes)   # 0–1, higher = more broadly selected

print(f"  {n_genes} genes in PBS dataset")

# ── Threshold lines (75th percentile of each averaged axis) ───────────────
thr_afr = df['mean_afr'].quantile(0.75)
thr_mel = df['mean_mel'].quantile(0.75)
print(f"  Threshold lines: mean_afr ≥ {thr_afr:.3f}  /  mean_mel ≥ {thr_mel:.3f}")

# ── Quadrant classification ────────────────────────────────────────────────
def quadrant(row):
    hi_a = row['mean_afr'] >= thr_afr
    hi_m = row['mean_mel'] >= thr_mel
    if hi_a and hi_m:   return 'both'
    if hi_a:             return 'african'
    if hi_m:             return 'melanesian'
    return 'neither'
df['quadrant'] = df.apply(quadrant, axis=1)

print("\n=== Quadrant counts ===")
print(df['quadrant'].value_counts())
for q in ['african', 'melanesian', 'both']:
    sub = df[df['quadrant'] == q]
    cats = sub['functional_category'].value_counts()
    print(f"\n{q.capitalize()} quadrant — top genes:")
    print(sub.nlargest(5, 'rank_sum')[['gene', 'functional_category',
                                        'mean_afr', 'mean_mel']].to_string(index=False))

print("\n=== Top 20 by rank-sum ===")
top20 = df.nlargest(20, 'rank_sum')[['gene', 'functional_category',
                                      'mean_afr', 'mean_mel', 'rank_pct']]
print(top20.to_string(index=False))

# ── Node size: rank_sum percentile → marker area ──────────────────────────
S_MIN, S_MAX = 40, 700
df['s_scatter'] = S_MIN + (S_MAX - S_MIN) * df['rank_pct'] ** 1.2

# ── Gene label set ────────────────────────────────────────────────────────
ALWAYS_LABEL = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R', 'OCA2', 'MLANA',
    'MITF', 'SOX10', 'PAX3', 'TFAP2A', 'KIT', 'KITLG', 'EDNRB',
}
top_rank = set(df.nlargest(12, 'rank_sum')['gene'])
label_genes = ALWAYS_LABEL | top_rank


# ══════════════════════════════════════════════════════════════════════════
# Figure
# ══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(20, 9))
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.35, 1],
                       left=0.06, right=0.97, top=0.88, bottom=0.12,
                       wspace=0.32)
ax_scatter = fig.add_subplot(gs[0])
ax_bar     = fig.add_subplot(gs[1])

# ══════════════════════════════════════════════════════════════════════════
# Panel A — 2-D scatter
# ══════════════════════════════════════════════════════════════════════════
# Threshold reference lines
ax_scatter.axvline(thr_afr, color='#CCCCCC', lw=1.2, linestyle='--', zorder=1)
ax_scatter.axhline(thr_mel, color='#CCCCCC', lw=1.2, linestyle='--', zorder=1)

# Diagonal y = x (equal signal in both populations)
xy_max = max(df['mean_afr'].max(), df['mean_mel'].max()) * 1.05
ax_scatter.plot([0, xy_max], [0, xy_max],
                color='#AAAAAA', lw=1.0, linestyle=':', zorder=1,
                label='Equal African / Melanesian signal')

# Scatter by functional category
for cat in CATEGORY_ORDER:
    sub = df[df['functional_category'] == cat]
    if len(sub) == 0:
        continue
    ax_scatter.scatter(sub['mean_afr'], sub['mean_mel'],
                       s=sub['s_scatter'],
                       c=CATEGORY_COLORS[cat],
                       alpha=0.82,
                       edgecolors='white', linewidths=0.5,
                       zorder=3)

# Gene labels
texts = []
top_rank_set = set(df.nlargest(12, 'rank_sum')['gene'])
for _, row in df.iterrows():
    if row['gene'] not in label_genes:
        continue
    is_top = row['gene'] in top_rank_set
    texts.append(
        ax_scatter.text(row['mean_afr'], row['mean_mel'],
                        row['gene'],
                        fontsize=8.5 if is_top else 7.5,
                        fontweight='bold' if is_top else 'normal',
                        style='italic',
                        color='#111111',
                        zorder=5)
    )
adjust_text(texts, ax=ax_scatter,
            arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.5),
            force_points=(2.0, 2.0), force_text=(2.0, 2.0),
            expand_points=(2.5, 2.5), expand_text=(2.0, 2.0),
            iter_lim=1000)

# Quadrant labels
quad_style = dict(fontsize=9, color='#555555', style='italic',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='#F8F8F8',
                            edgecolor='#DDDDDD', alpha=0.88))
x_max_plot = df['mean_afr'].max() * 1.05
y_max_plot = df['mean_mel'].max() * 1.05

ax_scatter.text(thr_afr * 0.35, y_max_plot * 0.94,
                'Melanesian-specific', ha='center', **quad_style)
ax_scatter.text(x_max_plot * 0.72, y_max_plot * 0.94,
                'Broadly selected\n(both populations)', ha='center', **quad_style)
ax_scatter.text(x_max_plot * 0.72, thr_mel * 0.25,
                'African-specific', ha='center', **quad_style)
ax_scatter.text(thr_afr * 0.35, thr_mel * 0.25,
                'Neither', ha='center', **quad_style)

# Size legend (rank-sum percentile scale)
size_handles = []
for pct, label in [(0.25, '25th'), (0.50, '50th'), (0.75, '75th'), (0.95, '95th')]:
    s = S_MIN + (S_MAX - S_MIN) * pct ** 1.2
    h = ax_scatter.scatter([], [], s=s, c='#888888', alpha=0.75,
                           edgecolors='white', linewidths=0.5,
                           label=f'{label} percentile')
    size_handles.append(h)
size_leg = ax_scatter.legend(handles=size_handles,
                              title='Node size =\nrank-sum percentile',
                              title_fontsize=8.5, fontsize=8,
                              loc='lower right', framealpha=0.92,
                              edgecolor='gray', borderpad=0.8,
                              handletextpad=1.0, labelspacing=0.7)
ax_scatter.add_artist(size_leg)

ax_scatter.set_xlim(-0.01, x_max_plot)
ax_scatter.set_ylim(-0.01, y_max_plot)
ax_scatter.set_xlabel('Mean African PBS\n[(PBS-1 + PBS-2) / 2]', fontsize=13)
ax_scatter.set_ylabel('Mean Melanesian PBS\n[(PBS-3 + PBS-4) / 2]', fontsize=13)
ax_scatter.set_title('African vs. Melanesian population-specific selection',
                     fontsize=12, fontweight='bold', loc='left', pad=10)
ax_scatter.tick_params(labelsize=10)
ax_scatter.text(-0.08, 1.06, 'A', transform=ax_scatter.transAxes,
                fontsize=22, fontweight='bold', va='top')


# ══════════════════════════════════════════════════════════════════════════
# Panel B — Diverging bar chart (top 20 by rank-sum)
# ══════════════════════════════════════════════════════════════════════════
top20 = (df.nlargest(20, 'rank_sum')
           .sort_values('rank_sum', ascending=True)
           .reset_index(drop=True))

y_pos = np.arange(len(top20))
BAR_HEIGHT = 0.38

# African bars (extend left = negative x)
ax_bar.barh(y_pos + BAR_HEIGHT / 2,
            -top20['mean_afr'],
            height=BAR_HEIGHT,
            color='#C0392B', alpha=0.85, label='Mean African PBS')

# Melanesian bars (extend right = positive x)
ax_bar.barh(y_pos - BAR_HEIGHT / 2,
            top20['mean_mel'],
            height=BAR_HEIGHT,
            color='#2471A3', alpha=0.85, label='Mean Melanesian PBS')

# Centre line
ax_bar.axvline(0, color='#333333', lw=1.2, zorder=4)

# Gene labels (centred on x=0), coloured by functional category
for i, row in top20.iterrows():
    cat = row['functional_category'] if pd.notna(row['functional_category']) else 'Other'
    color = CATEGORY_COLORS.get(cat, '#B0B0B0')
    ax_bar.text(0, i, row['gene'],
                ha='center', va='center',
                fontsize=8.5, fontstyle='italic', fontweight='bold',
                color=color, zorder=5,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          edgecolor='none', alpha=0.85))

# x-axis: show positive values on both sides
x_extent = top20[['mean_afr', 'mean_mel']].max().max() * 1.15
ax_bar.set_xlim(-x_extent, x_extent)
tick_vals = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
tick_vals = tick_vals[tick_vals <= x_extent]
ax_bar.set_xticks(np.concatenate([-tick_vals[1:][::-1], tick_vals]))
ax_bar.set_xticklabels(
    [f'{v:.1f}' for v in tick_vals[1:][::-1]] + [f'{v:.1f}' for v in tick_vals],
    fontsize=9)

ax_bar.set_yticks([])
ax_bar.set_xlabel('PBS value', fontsize=12)
ax_bar.set_title('Top 20 genes by cross-population rank-sum',
                 fontsize=12, fontweight='bold', loc='left', pad=10)

# Population bar legend
bar_handles = [
    mpatches.Patch(facecolor='#C0392B', alpha=0.85, label='Mean African PBS\n[(PBS-1+PBS-2)/2]'),
    mpatches.Patch(facecolor='#2471A3', alpha=0.85, label='Mean Melanesian PBS\n[(PBS-3+PBS-4)/2]'),
]
ax_bar.legend(handles=bar_handles, fontsize=8.5, loc='lower left',
              framealpha=0.9, edgecolor='gray')

# Population direction labels above the bars
ax_bar.text(-x_extent * 0.55, len(top20) - 0.2, '← African',
            fontsize=9, color='#C0392B', fontstyle='italic', ha='center')
ax_bar.text(x_extent * 0.55, len(top20) - 0.2, 'Melanesian →',
            fontsize=9, color='#2471A3', fontstyle='italic', ha='center')

ax_bar.tick_params(labelsize=9)
ax_bar.text(-0.10, 1.06, 'B', transform=ax_bar.transAxes,
            fontsize=22, fontweight='bold', va='top')

# ── Shared category legend (bottom) ───────────────────────────────────────
cat_patches = [
    mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85,
                   edgecolor='white', linewidth=0.5)
    for c in CATEGORY_ORDER
]
fig.legend(handles=cat_patches,
           fontsize=9.5, ncol=6,
           loc='lower center', bbox_to_anchor=(0.5, 0.002),
           framealpha=0.92, edgecolor='gray',
           title='Functional category (node color / gene label color)',
           title_fontsize=9.5)

fig.suptitle(
    'Cross-population PBS summary — African vs. Melanesian selection in the melanogenesis network\n'
    'Node size = rank-sum percentile across all 4 scans  |  Threshold lines = 75th percentile',
    fontsize=12, fontweight='bold', y=0.975)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_summary.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {path}")
plt.close(fig)
print("\nDone!")
