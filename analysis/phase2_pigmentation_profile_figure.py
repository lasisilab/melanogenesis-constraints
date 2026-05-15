"""
phase2_pigmentation_profile_figure.py

Focused replacement for the overwhelming 127-gene heatmap.

Shows only the 34 canonical 'pigmentation network' genes, sorted by their
data-driven cluster (toolkit → hub), with their oriented feature profile.
Adds a small lower panel listing the 8 non-canonical genes that share the
toolkit profile (candidate skin-context genes missed by curated lists).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

m = pd.read_csv(os.path.join(OUT_DIR, 'gene_profile_matrix.csv'))

PROFILE_COLS = ['centrality_inv', 'kegg_inv', 'tau_pct',
                'skin_top_inv', 'skin_z_pct', 'loeuf_pct']
SHORT_LABELS = {
    'centrality_inv': 'low\ncentrality',
    'kegg_inv':       'few KEGG\npathways',
    'tau_pct':        'high tau\n(tissue-specific)',
    'skin_top_inv':   'skin is\ntop tissue',
    'skin_z_pct':     'high skin\nz-score',
    'loeuf_pct':      'high LOEUF\n(LoF-tolerant)',
}
CLUSTER_LABEL = {
    1: 'toolkit',
    2: 'mid-broad',
    3: 'broad, not-skin',
    4: 'central hub',
}
CLUSTER_COLOR = {1: '#D94040', 2: '#E89060', 3: '#7AB6E8', 4: '#2A4878'}

canon = (m[m['canonical']]
         .sort_values(['cluster', 'composite'], ascending=[True, False])
         .reset_index(drop=True))
nontk = (m[(~m['canonical']) & (m['cluster'] == 1)]
         .sort_values('composite', ascending=False)
         .reset_index(drop=True))

print(f"canonical n={len(canon)}, non-canonical toolkit n={len(nontk)}")

def draw_block(ax, df, title, show_xlabels):
    n = len(df)
    heat = df[PROFILE_COLS].values
    im = ax.imshow(heat, aspect='auto', cmap='magma_r', vmin=0, vmax=1)
    ax.set_yticks(range(n))
    ax.set_yticklabels(df['gene'].values, fontsize=9)
    ax.tick_params(axis='y', length=0, pad=2)
    if show_xlabels:
        ax.set_xticks(range(len(PROFILE_COLS)))
        ax.set_xticklabels([SHORT_LABELS[c] for c in PROFILE_COLS],
                           fontsize=10, rotation=0)
        ax.tick_params(axis='x', length=0, pad=4)
    else:
        ax.set_xticks([])
    # cell value annotation
    for i in range(n):
        for j, c in enumerate(PROFILE_COLS):
            v = heat[i, j]
            ax.text(j, i, f'{v:.2f}', ha='center', va='center',
                    fontsize=7,
                    color='white' if v > 0.55 else '#222222')
    # cluster boundary lines
    for i in range(1, n):
        if df.iloc[i]['cluster'] != df.iloc[i - 1]['cluster']:
            ax.axhline(i - 0.5, color='white', lw=2)
    ax.set_title(title, fontsize=11, fontweight='bold', loc='left', pad=6)
    return im

def cluster_strip(ax, df):
    n = len(df)
    for i, c in enumerate(df['cluster']):
        ax.add_patch(Rectangle((0, i - 0.5), 1, 1,
                               color=CLUSTER_COLOR[c], lw=0))
    # cluster label at the centre of each cluster block
    seen = {}
    for i, c in enumerate(df['cluster']):
        seen.setdefault(c, []).append(i)
    for c, idxs in seen.items():
        mid = (min(idxs) + max(idxs)) / 2
        ax.text(0.5, mid, CLUSTER_LABEL[c], ha='center', va='center',
                fontsize=8, color='white', fontweight='bold', rotation=90)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def composite_bar(ax, df):
    n = len(df)
    bar_colors = [CLUSTER_COLOR[c] for c in df['cluster']]
    ax.barh(range(n), df['composite'].values, color=bar_colors,
            edgecolor='black', lw=0.3, height=0.85)
    ax.axvline(0.5, color='gray', ls='--', lw=0.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xticks([0, 0.5, 1.0])
    ax.tick_params(axis='x', labelsize=8)
    ax.set_xlabel('composite\ntoolkit score', fontsize=9)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

n_canon = len(canon)
n_nontk = len(nontk)

fig = plt.figure(figsize=(12, max(11, 0.32 * (n_canon + n_nontk + 4))))
gs = fig.add_gridspec(2, 4,
                      height_ratios=[n_canon, n_nontk + 1],
                      width_ratios=[0.10, 1.0, 0.32, 0.06],
                      hspace=0.10, wspace=0.05)

# Top: canonical pigmentation genes
ax_c1 = fig.add_subplot(gs[0, 0])
cluster_strip(ax_c1, canon)
ax_h1 = fig.add_subplot(gs[0, 1])
im = draw_block(ax_h1, canon,
                f'Canonical "pigmentation network" genes (n={n_canon}) — '
                f'sorted by data-driven cluster (toolkit → hub)',
                show_xlabels=False)
ax_b1 = fig.add_subplot(gs[0, 2])
composite_bar(ax_b1, canon)

# Bottom: non-canonical toolkit-profile genes
ax_c2 = fig.add_subplot(gs[1, 0])
cluster_strip(ax_c2, nontk)
ax_h2 = fig.add_subplot(gs[1, 1])
draw_block(ax_h2, nontk,
           f'Non-canonical genes sharing the toolkit profile (cluster 1, '
           f'n={n_nontk}) — candidate skin-context genes',
           show_xlabels=True)
ax_b2 = fig.add_subplot(gs[1, 2])
composite_bar(ax_b2, nontk)

# Shared colorbar
cax = fig.add_subplot(gs[:, 3])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('percentile across the 127 network genes\n'
             'dark = matches the column header\n'
             '(e.g. dark in "low centrality" = this gene IS low-centrality)',
             fontsize=8.5)
cb.ax.tick_params(labelsize=8)

# Legend for cluster colours along the bottom
from matplotlib.patches import Patch
handles = [Patch(color=CLUSTER_COLOR[c],
                 label=f'cluster {c}: {CLUSTER_LABEL[c]}')
           for c in [1, 2, 3, 4]]
fig.legend(handles=handles, loc='lower center', ncol=4,
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.005))

fig.suptitle(
    'Where canonical "pigmentation" genes actually land in the network profile\n'
    'a row of all-dark cells = a pure "toolkit" gene  (peripheral, narrow, '
    'tissue-specific, LoF-tolerant)',
    fontsize=12, fontweight='bold', y=0.995)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'figure_pigmentation_profile.{ext}'),
                dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("saved figure_pigmentation_profile.png/pdf")
