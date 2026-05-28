"""
phase2_pbs_within_network_smallfigs.py

Compact poster versions of the two within-network forest figures: all four
branch panels (AFR/MEL/EUR/SAS) in a single horizontal row, sized 6" × 3".

Reads the already-computed coefficient tables (no permutations re-run):
  output/phase2_pbs_within_network.csv           → SGDP genome-wide percentile
  output/phase2_pbs_within_network_combined.csv  → combined-dataset within-network

Outputs:
  figure_phase2_pbs_within_network_row.{png,pdf}
  figure_phase2_pbs_within_network_combined_row.{png,pdf}
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

PRED_ORDER = ['tau', 'betweenness_centrality', 'kegg_log1p', 'LOEUF', 'tissue_breadth']
PRED_SHORT = {            # short labels for the tiny left axis
    'tau': 'τ', 'betweenness_centrality': 'Betw.', 'kegg_log1p': 'KEGG',
    'LOEUF': 'LOEUF', 'tissue_breadth': 'Breadth',
}
PRED_COLORS = {
    'tau': '#d62728', 'betweenness_centrality': '#9467bd',
    'kegg_log1p': '#2ca02c', 'LOEUF': '#7f7f7f', 'tissue_breadth': '#1f77b4',
}
BRANCH_ORDER = ['AFR', 'MEL', 'EUR', 'SAS']


def make_row_figure(csv_path, fname, title):
    df = pd.read_csv(csv_path)
    fig, axes = plt.subplots(1, 4, figsize=(11, 2.1), sharey=True)

    for ax, code in zip(axes, BRANCH_ORDER):
        bdf = df[df['branch'] == code].set_index('predictor')
        yvals = list(range(len(PRED_ORDER)))[::-1]
        ax.axvline(0, color='#444', lw=0.6, ls='--', zorder=0)
        for yp, pred in zip(yvals, PRED_ORDER):
            if pred not in bdf.index:
                continue
            r = bdf.loc[pred]
            c = PRED_COLORS[pred]
            ax.errorbar(r['adj_beta'], yp,
                        xerr=[[r['adj_beta'] - r['adj_ci_lo']],
                              [r['adj_ci_hi'] - r['adj_beta']]],
                        fmt='o', color=c, ecolor=c, elinewidth=1.0, capsize=2,
                        markersize=4, markeredgecolor='black', markeredgewidth=0.3)
            if r['perm_p'] < 0.05:
                ax.text(r['adj_ci_hi'], yp, '*', fontsize=7, va='center',
                        ha='left', color='black')
        ax.set_title(code, fontsize=8, fontweight='bold', pad=3)
        ax.tick_params(axis='x', labelsize=5.5, pad=1)
        ax.tick_params(axis='y', length=0)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

    axes[0].set_yticks(list(range(len(PRED_ORDER)))[::-1])
    axes[0].set_yticklabels([PRED_SHORT[p] for p in PRED_ORDER], fontsize=6.5)

    fig.supxlabel('Adjusted β  (Δ percentile per SD, 95% CI;  * perm p<0.05)',
                  fontsize=6.5, y=0.02)
    fig.suptitle(title, fontsize=8, fontweight='bold', y=1.0)
    fig.tight_layout(rect=[0, 0.04, 1, 0.97])

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'{fname}.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f'  Saved → {path}')
    plt.close(fig)


make_row_figure(
    os.path.join(OUT_DIR, 'phase2_pbs_within_network.csv'),
    'figure_phase2_pbs_within_network_row',
    'Within-network predictors of PBS — SGDP genome-wide percentile')

make_row_figure(
    os.path.join(OUT_DIR, 'phase2_pbs_within_network_combined.csv'),
    'figure_phase2_pbs_within_network_combined_row',
    'Within-network predictors of PBS — combined dataset (HGDP+1KGP+SGDP)')
