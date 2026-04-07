"""
phase1_phylop_analysis.py

Phase 1.1: PhyloP evolutionary conservation scores for the Raghunath
melanogenesis network.

Merges mean PhyloP 100-way vertebrate conservation scores (data/phylop_scores.csv)
with the LOEUF network data (LOEUF_by_functional_category.xlsx) and produces:

  1. data/network_constraint_phylop.csv  — merged dataset (LOEUF + PhyloP)
  2. output/figure_phase1_phylop.png/pdf — two-panel figure:
       Panel A: PhyloP vs. LOEUF scatter (colored by functional category)
       Panel B: PhyloP by functional category (boxplot)

Expected relationship:
  Higher PhyloP → more conserved → lower LOEUF (more LoF-intolerant)
  i.e., Spearman ρ(PhyloP, LOEUF) should be NEGATIVE.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from adjustText import adjust_text
import scikit_posthocs as sp

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOEUF_FILE   = os.path.join(PROJECT_DIR, 'LOEUF_by_functional_category.xlsx')
PHYLOP_FILE  = os.path.join(PROJECT_DIR, 'data', 'phylop_scores.csv')
OUT_CSV      = os.path.join(PROJECT_DIR, 'data', 'network_constraint_phylop.csv')
OUT_DIR      = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load and merge ─────────────────────────────────────────────────────────
print("Loading LOEUF network data...")
loeuf = pd.read_excel(LOEUF_FILE, sheet_name='All Genes by Category')
loeuf.columns = ['gene', 'functional_category', 'disease_class',
                  'LOEUF', 'pLI', 'betweenness_centrality']
loeuf['gene'] = loeuf['gene'].str.upper()

print("Loading PhyloP scores...")
phylop = pd.read_csv(PHYLOP_FILE)[['gene', 'mean_phylop_100way']]
phylop['gene'] = phylop['gene'].str.upper()

df = loeuf.merge(phylop, on='gene', how='left')
n_missing = df['mean_phylop_100way'].isna().sum()
print(f"Merged: {len(df)} genes, {n_missing} missing PhyloP scores")
if n_missing > 0:
    print("  Missing:", df.loc[df['mean_phylop_100way'].isna(), 'gene'].tolist())

df.to_csv(OUT_CSV, index=False)
print(f"Saved merged data → {OUT_CSV}")

# ── Color scheme (consistent with generate_figures.py) ────────────────────
CATEGORY_ORDER = [
    'Pigment-specific',
    'Developmental/NC',
    'Generic signaling',
    'Cytokines/growth factors',
    'Apoptosis/cell death',
    'Other',
]

CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',   # Red
    'Developmental/NC':          '#E8907E',   # Salmon
    'Generic signaling':         '#F5C242',   # Gold
    'Cytokines/growth factors':  '#4878CF',   # Blue
    'Apoptosis/cell death':      '#6BAD6B',   # Green
    'Other':                     '#B0B0B0',   # Gray
}

CATEGORY_SHORT = {
    'Pigment-specific':          'Pigment-\nspecific',
    'Developmental/NC':          'Developmental\n/NC',
    'Generic signaling':         'Generic\nsignaling',
    'Cytokines/growth factors':  'Cytokines/\nGF',
    'Apoptosis/cell death':      'Apoptosis/\ncell death',
    'Other':                     'Other',
}

# ── Statistics ────────────────────────────────────────────────────────────
df_clean = df.dropna(subset=['mean_phylop_100way', 'LOEUF'])
rho, pval = stats.spearmanr(df_clean['mean_phylop_100way'], df_clean['LOEUF'])
print(f"\nSpearman ρ(PhyloP, LOEUF) = {rho:.3f}, p = {pval:.3e}")

groups = [df_clean.loc[df_clean['functional_category'] == c, 'mean_phylop_100way'].values
          for c in CATEGORY_ORDER if c in df_clean['functional_category'].values]
kw_stat, kw_p = stats.kruskal(*[g for g in groups if len(g) > 0])
print(f"Kruskal-Wallis PhyloP across categories: H = {kw_stat:.2f}, p = {kw_p:.3e}")

dunn = sp.posthoc_dunn(
    df_clean, val_col='mean_phylop_100way', group_col='functional_category',
    p_adjust='bonferroni')
print("\nDunn's posthoc (Bonferroni):")
cats_present = [c for c in CATEGORY_ORDER if c in df_clean['functional_category'].values]
for i in range(len(cats_present)):
    for j in range(i + 1, len(cats_present)):
        p = dunn.loc[cats_present[i], cats_present[j]]
        if p < 0.05:
            print(f"  {cats_present[i]} vs {cats_present[j]}: p = {p:.4f}")

# ── Figure ─────────────────────────────────────────────────────────────────
LABEL_GENES = {'TYR', 'TYRP1', 'DCT', 'OCA2', 'MC1R', 'SOX10', 'MITF',
               'PAX3', 'TFAP2A', 'FOS', 'DUSP1', 'JUN', 'AKT1', 'TP53'}

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.subplots_adjust(wspace=0.32, left=0.07, right=0.98, top=0.88, bottom=0.38)

# ── Panel A: PhyloP vs. LOEUF scatter ─────────────────────────────────────
ax = axes[0]
for cat in CATEGORY_ORDER:
    sub = df_clean[df_clean['functional_category'] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(sub['mean_phylop_100way'], sub['LOEUF'],
               c=CATEGORY_COLORS[cat], s=65, alpha=0.85,
               edgecolors='white', linewidths=0.4,
               label=f'{cat} (n={len(sub)})', zorder=3)

# Trend line
slope, intercept, *_ = stats.linregress(df_clean['mean_phylop_100way'], df_clean['LOEUF'])
x_line = np.linspace(df_clean['mean_phylop_100way'].min(),
                     df_clean['mean_phylop_100way'].max(), 100)
ax.plot(x_line, slope * x_line + intercept, '--',
        color='#555555', alpha=0.6, lw=1, zorder=1)

# Gene labels
texts = []
for _, row in df_clean.iterrows():
    if row['gene'] in LABEL_GENES:
        texts.append(ax.text(row['mean_phylop_100way'], row['LOEUF'],
                             row['gene'], fontsize=13, fontweight='bold',
                             color='#333333', style='italic'))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7),
            force_points=(2.5, 2.5), force_text=(2.5, 2.5),
            expand_points=(3.0, 3.0), expand_text=(2.0, 2.0),
            iter_lim=1000)

ax.text(0.98, -0.20, f'Spearman ρ = {rho:.3f},  p = {pval:.2e}',
        transform=ax.transAxes, fontsize=13, ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.9))

ax.set_xlabel('Mean PhyloP 100-way\n(higher = more conserved)', fontsize=16)
ax.set_ylabel('LOEUF (lower = more constrained)', fontsize=16)
ax.set_title('Evolutionary conservation vs.\nLoF intolerance', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=11, loc='upper center', bbox_to_anchor=(0.5, -0.32),
          ncol=2, framealpha=0.9, edgecolor='gray',
          handletextpad=0.4, borderpad=0.5, markerscale=1.2)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')


# ── Panel B: PhyloP boxplot by functional category ────────────────────────
ax = axes[1]
bp_data, colors, ns = [], [], []
for cat in CATEGORY_ORDER:
    vals = df_clean.loc[df_clean['functional_category'] == cat, 'mean_phylop_100way'].values
    bp_data.append(vals)
    colors.append(CATEGORY_COLORS[cat])
    ns.append(len(vals))

positions = list(range(len(CATEGORY_ORDER)))
bp = ax.boxplot(bp_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=False, zorder=2,
                medianprops=dict(color='black', linewidth=1.5),
                whiskerprops=dict(color='gray'), capprops=dict(color='gray'))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)
    patch.set_edgecolor('gray')

rng = np.random.default_rng(42)
for i, (vals, cat) in enumerate(zip(bp_data, CATEGORY_ORDER)):
    if len(vals) == 0:
        continue
    jitter = rng.uniform(-0.15, 0.15, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals,
               c=CATEGORY_COLORS[cat], s=25, alpha=0.8,
               edgecolors='white', linewidths=0.3, zorder=3)
    ax.text(i, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else -0.9,
            f'n={ns[i]}', ha='center', fontsize=11, color='#555555',
            fontstyle='italic', zorder=4)


def draw_bracket(ax, x1, x2, y, p_text, lw=1.0, fontsize=14):
    h = 0.03
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], color='black', lw=lw, zorder=5)
    ax.text((x1 + x2) / 2, y + 0.01, p_text,
            ha='center', va='bottom', fontsize=fontsize, zorder=5)


sig_pairs = []
for i in range(len(CATEGORY_ORDER)):
    for j in range(i + 1, len(CATEGORY_ORDER)):
        c1, c2 = CATEGORY_ORDER[i], CATEGORY_ORDER[j]
        if c1 in dunn.index and c2 in dunn.columns:
            p = dunn.loc[c1, c2]
            if p < 0.05:
                sig_pairs.append((i, j, p))
sig_pairs.sort(key=lambda x: x[1] - x[0])

y_vals = [v for vals in bp_data for v in vals]
y_max = max(y_vals) if y_vals else 1.5
bracket_y_start = y_max + 0.08
bracket_spacing = 0.15
for k, (i, j, p) in enumerate(sig_pairs):
    y_bracket = bracket_y_start + k * bracket_spacing
    p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
    draw_bracket(ax, i, j, y_bracket, p_text)

top_y = bracket_y_start + max(len(sig_pairs), 1) * bracket_spacing + 0.05
ax.set_ylim(ax.get_ylim()[0], top_y + 0.2)
ax.text(0.98, 0.98, f'Kruskal-Wallis p = {kw_p:.2e}',
        transform=ax.transAxes, fontsize=13, ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.9))

ax.axhline(0, color='#aaaaaa', lw=0.8, ls=':', zorder=0)
ax.set_xticks(positions)
ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                   fontsize=13, rotation=30, ha='right')
ax.set_ylabel('Mean PhyloP 100-way\n(higher = more conserved)', fontsize=16)
ax.set_title('PhyloP conservation by\nfunctional category', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

# Save
for ext in ('png', 'pdf'):
    out_path = os.path.join(OUT_DIR, f'figure_phase1_phylop.{ext}')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {out_path}")
plt.close(fig)

# ── Summary table ─────────────────────────────────────────────────────────
print("\nPhyloP summary by functional category:")
summary = (df_clean.groupby('functional_category')['mean_phylop_100way']
           .agg(['count', 'median', 'mean'])
           .round(3)
           .rename(columns={'count': 'N', 'median': 'Median PhyloP', 'mean': 'Mean PhyloP'}))
summary = summary.loc[[c for c in CATEGORY_ORDER if c in summary.index]]
print(summary.to_string())

print("\nDone!")
