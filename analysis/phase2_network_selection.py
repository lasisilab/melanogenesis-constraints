"""
phase2_network_selection.py

Network position vs. population-specific selection — two-panel figure
for the PEQG 2026 poster.

Replicates the Raghunath-style scatter (betweenness centrality × KEGG
pathway count) but encodes PBS value as node size, directly testing:

  "Do population-specific selection signatures concentrate at peripheral
   network positions while core pathway genes remain universally constrained?"

Layout (1 × 2):
  Panel A: Within-pathway centrality (Y) × Cross-system connectivity (X)
           Node size ∝ PBS-1 (African selection, S. Asian outgroup)
  Panel B: Same axes
           Node size ∝ PBS-3 (Melanesian selection, S. Asian outgroup)

Node color = PEQG poster functional category (6 categories).
Labeled genes: key pigment genes + top PBS genes per panel.

Output: output/figure_phase2_network_selection.png / .pdf
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

df = (master[['gene', 'functional_category', 'LOEUF', 'betweenness_centrality']]
      .merge(kegg, on='gene', how='left')
      .merge(pbs,  on='gene', how='left'))

df = df.dropna(subset=['betweenness_centrality', 'kegg_pathway_count'])
print(f"  {len(df)} genes with betweenness + KEGG data")

# ── Square-root transform for plotting (spread skewed distributions) ──────
# Most genes cluster near 0 on both axes; sqrt spreads the mass more evenly.
df['x_plot'] = np.sqrt(df['kegg_pathway_count'].clip(lower=0))
df['y_plot'] = np.sqrt(df['betweenness_centrality'].clip(lower=0))

# ── Node size scaling ──────────────────────────────────────────────────────
# Map PBS → marker area (s parameter in scatter).
# Min size = 60 (PBS = 0, still visible), max size = 900 (top PBS).
S_MIN, S_MAX = 60, 900

def pbs_to_size(pbs_col):
    vals = df[pbs_col].fillna(0)
    pbs_max = vals.max()
    if pbs_max == 0:
        return np.full(len(vals), S_MIN)
    return S_MIN + (S_MAX - S_MIN) * (vals / pbs_max) ** 0.6

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
print("\n=== Spearman correlations (within-network) ===")
for pbs_col, label in [('pbs1_african', 'PBS-1 African'), ('pbs3_melanesian', 'PBS-3 Melanesian')]:
    for x_col, x_label in [('betweenness_centrality', 'betweenness'),
                            ('kegg_pathway_count', 'KEGG count')]:
        sub = df.dropna(subset=[pbs_col, x_col])
        r, p = stats.spearmanr(sub[x_col], sub[pbs_col])
        print(f"  {label} × {x_label}: ρ = {r:.3f}, p = {p:.3e}")

print("\n=== Top 10 genes per PBS scan ===")
for col, label in [('pbs1_african', 'PBS-1 African'), ('pbs3_melanesian', 'PBS-3 Melanesian')]:
    top = df.nlargest(10, col)[['gene', 'functional_category', col,
                                 'betweenness_centrality', 'kegg_pathway_count']]
    print(f"\n{label}:")
    print(top.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════
# Figure
# ══════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(20, 9))
fig.subplots_adjust(wspace=0.30, left=0.07, right=0.97, top=0.88, bottom=0.10)

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

    # ── Quadrant annotations — positions in sqrt-transformed coordinates ───
    xp_max = df['x_plot'].max()   # sqrt(120) ≈ 10.95
    yp_max = df['y_plot'].max()   # sqrt(0.21) ≈ 0.458
    xp_mid = xp_max * 0.45

    quad_style = dict(fontsize=9.5, color='#444444', ha='center',
                      style='italic',
                      bbox=dict(boxstyle='round,pad=0.35', facecolor='#F8F8F8',
                                edgecolor='#CCCCCC', alpha=0.85))

    ax.text(xp_mid * 0.28, yp_max * 0.88,
            'Within-pathway hubs\nPathway-essential; constrained',
            **quad_style)
    ax.text(xp_mid * 2.55, yp_max * 0.88,
            'Pleiotropic hubs\nCross-system essential;\nmost constrained',
            **quad_style)
    ax.text(xp_mid * 0.28, yp_max * 0.10,
            'Peripheral specialists\nPathway-effectors;\nmost variable',
            **quad_style)

    # ── Diagonal arrow (low connectivity → high connectivity) ──────────────
    ax.annotate('', xy=(xp_max * 0.80, yp_max * 0.78),
                xytext=(xp_max * 0.18, yp_max * 0.18),
                arrowprops=dict(arrowstyle='->', color='#AAAAAA',
                                lw=1.5, linestyle='dashed'))

    # ── Axis ticks — show original (untransformed) values ─────────────────
    x_ticks_orig = [0, 10, 20, 40, 60, 80, 100, 120]
    y_ticks_orig = [0.00, 0.02, 0.05, 0.10, 0.15, 0.20]
    ax.set_xticks([np.sqrt(v) for v in x_ticks_orig])
    ax.set_xticklabels(x_ticks_orig)
    ax.set_yticks([np.sqrt(v) for v in y_ticks_orig])
    ax.set_yticklabels(y_ticks_orig)

    # ── Axis labels and titles ─────────────────────────────────────────────
    ax.set_xlabel('Cross-system connectivity\n(number of KEGG pathways)',
                  fontsize=13)
    ax.set_ylabel('Within-pathway centrality\n(betweenness centrality)',
                  fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)

    # HIGH / LOW axis annotations
    ax.text(0.01, 0.98, 'HIGH', transform=ax.transAxes,
            fontsize=9, color='#888888', va='top', ha='left', fontstyle='italic')
    ax.text(0.01, 0.01, 'LOW', transform=ax.transAxes,
            fontsize=9, color='#888888', va='bottom', ha='left', fontstyle='italic')
    ax.text(0.01, -0.06, 'LOW', transform=ax.transAxes,
            fontsize=9, color='#888888', va='top', ha='left', fontstyle='italic')
    ax.text(0.99, -0.06, 'HIGH', transform=ax.transAxes,
            fontsize=9, color='#888888', va='top', ha='right', fontstyle='italic')

    # Panel letter
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')

    # ── Size legend (PBS scale) ────────────────────────────────────────────
    pbs_vals_legend = [0.0, 0.2, 0.4, 0.6, 0.8]
    pbs_max = df[pbs_col].max()
    size_legend_handles = []
    for pv in pbs_vals_legend:
        if pv > pbs_max + 0.05:
            continue
        s = S_MIN + (S_MAX - S_MIN) * (pv / pbs_max) ** 0.6
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
           loc='lower center', bbox_to_anchor=(0.5, 0.005),
           framealpha=0.92, edgecolor='gray',
           title='Functional category (node color)',
           title_fontsize=9.5)

fig.suptitle(
    'Network position vs. population-specific selection in the melanogenesis network\n'
    'Node size = PBS value  |  Larger nodes = stronger population-specific selection signal',
    fontsize=12, fontweight='bold', y=0.975)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_network_selection.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {path}")
plt.close(fig)

print("\nDone!")
