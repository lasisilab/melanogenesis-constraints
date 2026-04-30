"""
phase2_network_constraint_loeuf.py

LoF intolerance vs. population-specific selection — two-panel figure
for the PEQG 2026 poster.

Variant on phase2_network_selection.py using LOEUF instead of betweenness
centrality to test whether constrained genes show different PBS patterns.

Layout (1 × 2):
  Panel A: LoF intolerance (Y, higher = less constrained) × Cross-system connectivity (X)
           Node size ∝ PBS-1 (African selection, S. Asian outgroup)
  Panel B: Same axes
           Node size ∝ PBS-3 (Melanesian selection, S. Asian outgroup)

Node color = PEQG poster functional category (6 categories).
Labeled genes: key pigment genes + top PBS genes per panel.

Output: output/figure_phase2_network_constraint_loeuf.png / .pdf
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from adjustText import adjust_text

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_gtex.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data',   'kegg_pathway_counts.csv')
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

kegg = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()

pbs = pd.read_csv(PBS_CSV)[['gene', 'pbs1_african', 'pbs3_melanesian']]
pbs['gene'] = pbs['gene'].str.upper()
pbs['pbs1_african']    = pbs['pbs1_african'].clip(lower=0)
pbs['pbs3_melanesian'] = pbs['pbs3_melanesian'].clip(lower=0)

df = (master[['gene', 'functional_category', 'LOEUF']]
      .merge(kegg, on='gene', how='left')
      .merge(pbs,  on='gene', how='left'))

df = df.dropna(subset=['LOEUF', 'kegg_pathway_count'])
print(f"  {len(df)} genes with LOEUF + KEGG data")

# ── Square-root transform for plotting (spread skewed distributions) ──────
# KEGG is skewed; sqrt spreads the mass more evenly.
# LOEUF is not transformed.
df['x_plot'] = np.sqrt(df['kegg_pathway_count'].clip(lower=0))
df['y_plot'] = df['LOEUF']

# ── Node size scaling ──────────────────────────────────────────────────────
# Map PBS → marker area (s parameter in scatter).
# Use global max across both PBS columns for consistent scale between panels.
S_MIN, S_MAX = 20, 1400
pbs_global_max = max(df['pbs1_african'].max(), df['pbs3_melanesian'].max())

def pbs_to_size(pbs_col):
    vals = df[pbs_col].fillna(0)
    if pbs_global_max == 0:
        return np.full(len(vals), S_MIN)
    return S_MIN + (S_MAX - S_MIN) * (vals / pbs_global_max) ** 2.0

df['s_afr'] = pbs_to_size('pbs1_african')
df['s_mel'] = pbs_to_size('pbs3_melanesian')

# ── Genes to always label ──────────────────────────────────────────────────
ALWAYS_LABEL = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R', 'OCA2', 'MLANA',
    'MITF', 'SOX10', 'PAX3', 'TFAP2A', 'KIT', 'KITLG', 'EDNRB',
    'NFKB1', 'AKT1', 'MAPK1',
}

def top_pbs_genes(col, n=8):
    return set(df.nlargest(n, col)['gene'])


# ── Print statistics ───────────────────────────────────────────────────────
print("\n=== Spearman correlations (LOEUF vs. PBS and KEGG) ===")
for pbs_col, label in [('pbs1_african', 'PBS-1 African'), ('pbs3_melanesian', 'PBS-3 Melanesian')]:
    for x_col, x_label in [('LOEUF', 'LOEUF'),
                            ('kegg_pathway_count', 'KEGG count')]:
        sub = df.dropna(subset=[pbs_col, x_col])
        r, p = stats.spearmanr(sub[x_col], sub[pbs_col])
        print(f"  {label} × {x_label}: ρ = {r:.3f}, p = {p:.3e}")

print("\n=== Top 10 genes per PBS scan ===")
for col, label in [('pbs1_african', 'PBS-1 African'), ('pbs3_melanesian', 'PBS-3 Melanesian')]:
    top = df.nlargest(10, col)[['gene', 'functional_category', col,
                                 'LOEUF', 'kegg_pathway_count']]
    print(f"\n{label}:")
    print(top.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════
# Figure
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
fig.subplots_adjust(wspace=0.30, left=0.07, right=0.97, top=0.88, bottom=0.18)

PANEL_SPECS = [
    ('pbs1_african',    's_afr', 'A',
     'PBS-1: African population-specific selection\n(S. Asian outgroup, Melanesian distant)'),
    ('pbs3_melanesian', 's_mel', 'B',
     'PBS-3: Melanesian population-specific selection\n(S. Asian outgroup, African distant)'),
]

for ax, (pbs_col, size_col, letter, subtitle) in zip(axes, PANEL_SPECS):

    label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col, n=8)

    # ── Draw points by category (so legend order is preserved) ────────────
    for cat in CATEGORY_ORDER:
        sub = df[df['functional_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub['x_plot'],
                   sub['y_plot'],
                   s=sub[size_col],
                   c=CATEGORY_COLORS[cat],
                   alpha=0.82,
                   edgecolors='white', linewidths=0.5,
                   zorder=3,
                   label=f'{cat} (n={len(sub)})')

    # ── Gene labels ────────────────────────────────────────────────────────
    texts = []
    for _, row in df.iterrows():
        if row['gene'] not in label_genes:
            continue
        # Bold + larger for top PBS genes; normal for always-label genes
        is_top = row['gene'] in top_pbs_genes(pbs_col, n=8)
        texts.append(
            ax.text(row['x_plot'],
                    row['y_plot'],
                    row['gene'],
                    fontsize=8.5 if is_top else 7.5,
                    fontweight='bold' if is_top else 'normal',
                    style='italic',
                    color='#111111',
                    zorder=5)
        )
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.5),
                force_points=(1.8, 1.8), force_text=(1.8, 1.8),
                expand_points=(2.5, 2.5), expand_text=(2.0, 2.0),
                iter_lim=1000)


    # ── Axis ticks — evenly-spaced in sqrt-plot space for X, regular for Y
    x_ticks_sqrt = np.linspace(0, df['x_plot'].max(), 6)
    x_tick_labels = [f'{v**2:.0f}' for v in x_ticks_sqrt]
    ax.set_xticks(x_ticks_sqrt)
    ax.set_xticklabels(x_tick_labels)

    # Y-axis ticks: evenly spaced LOEUF values
    y_min, y_max = df['y_plot'].min(), df['y_plot'].max()
    y_ticks = np.linspace(y_min, y_max, 6)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{v:.2f}' for v in y_ticks])

    # ── Axis labels and titles ─────────────────────────────────────────────
    ax.set_xlabel('Cross-system connectivity\n(number of KEGG pathways)',
                  fontsize=13)
    ax.set_ylabel('LOEUF (higher = less constrained)',
                  fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)


    # Panel letter
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')

    # ── Size legend (PBS scale) ────────────────────────────────────────────
    pbs_vals_legend = [0.0, 0.2, 0.4, 0.6, 0.8]
    size_legend_handles = []
    for pv in pbs_vals_legend:
        if pv > pbs_global_max + 0.05:
            continue
        s = S_MIN + (S_MAX - S_MIN) * (pv / pbs_global_max) ** 2.0
        h = ax.scatter([], [], s=s, c='#888888', alpha=0.75,
                       edgecolors='white', linewidths=0.5,
                       label=f'PBS = {pv:.1f}')
        size_legend_handles.append(h)

    size_leg = ax.legend(handles=size_legend_handles,
                         title='Node size = PBS',
                         title_fontsize=9, fontsize=8.5,
                         loc='lower right', framealpha=0.92,
                         edgecolor='gray', borderpad=0.8,
                         handletextpad=1.0, labelspacing=0.8)
    ax.add_artist(size_leg)

# ── Shared category legend (bottom, full width) ────────────────────────────
cat_patches = [
    mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85,
                   edgecolor='white', linewidth=0.5)
    for c in CATEGORY_ORDER
]
fig.legend(handles=cat_patches,
           fontsize=9.5, ncol=6,
           loc='lower center', bbox_to_anchor=(0.5, -0.02),
           framealpha=0.92, edgecolor='gray',
           title='Functional category (node color)',
           title_fontsize=9.5)

fig.suptitle(
    'LoF intolerance vs. population-specific selection in the melanogenesis network\n'
    'Node size = PBS value  |  Larger nodes = stronger population-specific selection signal',
    fontsize=12, fontweight='bold', y=0.975)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_network_constraint_loeuf.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {path}")
plt.close(fig)

print("\nDone!")
