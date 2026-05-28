"""
phase2_pbs_scatter.py

Per-population PBS scatter plots against a single predictor, so all four
populations can be compared on a consistent x-axis. Produces three figures
(one per predictor), each a 2×2 grid of populations:

  figure_phase2_pbs_scatter_tau.{png,pdf}         x = tissue specificity (τ)
  figure_phase2_pbs_scatter_breadth.{png,pdf}     x = tissue breadth (# tissues)
  figure_phase2_pbs_scatter_betweenness.{png,pdf} x = betweenness centrality

Populations (PBS outcomes):
  African (PBS-1, S.Asian outgroup) · Melanesian (PBS-3, S.Asian outgroup)
  European · South Asian  (both vs the Afr–Mel reference pair; no Eur–SAS FST)

Each panel: points colored by functional category, OLS trend line, Spearman ρ,
and gene-name labels on the outliers (top genes by PBS + extreme-x genes).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy import stats
from adjustText import adjust_text

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data',   'kegg_pathway_counts.csv')
BREADTH_CSV = os.path.join(PROJECT_DIR, 'data',   'gtex_tissue_breadth.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

CATEGORY_ORDER = [
    'Pigment-specific', 'Developmental/NC', 'Generic signaling',
    'Cytokines/growth factors', 'Apoptosis/cell death', 'Other',
]
CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',
    'Developmental/NC':          '#E8907E',
    'Generic signaling':         '#F5C242',
    'Cytokines/growth factors':  '#4878CF',
    'Apoptosis/cell death':      '#6BAD6B',
    'Other':                     '#B0B0B0',
}

# ── Load + merge (mirrors phase2_pbs_regression.py) ────────────────────────
master  = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
kegg    = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()
breadth = pd.read_csv(BREADTH_CSV)[['gene', 'tissue_breadth']]
breadth['gene'] = breadth['gene'].str.upper()
pbs_full = pd.read_csv(PBS_CSV)
pbs_full['gene'] = pbs_full['gene'].str.upper()


def _t(fst):
    fst = np.clip(pd.to_numeric(fst, errors='coerce').fillna(0).clip(lower=0), 0, 0.9999)
    return -np.log(1.0 - fst)


t_afr_eur = _t(pbs_full['fst_african_european'])
t_afr_mel = _t(pbs_full['fst_african_melanesian'])
t_afr_sas = _t(pbs_full['fst_african_southasian'])
t_eur_mel = _t(pbs_full['fst_european_melanesian'])
t_mel_sas = _t(pbs_full['fst_melanesian_southasian'])
pbs_full['pbs_european']   = ((t_afr_eur + t_eur_mel - t_afr_mel) / 2).clip(lower=0)
pbs_full['pbs_southasian'] = ((t_afr_sas + t_mel_sas - t_afr_mel) / 2).clip(lower=0)

pbs = pbs_full[['gene', 'pbs1_african', 'pbs3_melanesian',
                'pbs_european', 'pbs_southasian']].copy()
for c in ['pbs1_african', 'pbs3_melanesian', 'pbs_european', 'pbs_southasian']:
    pbs[c] = pbs[c].clip(lower=0)

df = (master[['gene', 'functional_category', 'betweenness_centrality', 'tau',
              'tissue_breadth']]
      .rename(columns={'tissue_breadth': 'tissue_breadth_master'})
      .merge(kegg,    on='gene', how='left')
      .merge(breadth, on='gene', how='left')   # tissue_breadth from gtex file
      .merge(pbs,     on='gene', how='inner'))

POPULATIONS = [
    ('pbs1_african',    'African  (PBS-1, S.Asian outgroup)'),
    ('pbs3_melanesian', 'Melanesian  (PBS-3, S.Asian outgroup)'),
    ('pbs_european',    'European  (Afr–Mel reference)'),
    ('pbs_southasian',  'South Asian  (Afr–Mel reference)'),
]

PREDICTORS = [
    ('tau',                    'Tissue specificity (τ)',       'tau'),
    ('tissue_breadth',         'Tissue breadth (# tissues TPM>1)', 'breadth'),
    ('betweenness_centrality', 'Betweenness centrality',       'betweenness'),
]

N_TOP_PBS = 8   # label this many top-PBS genes per panel
N_TOP_X   = 3   # plus this many extreme-x genes


def make_scatter_figure(xcol, xlabel, tag):
    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    axes = axes.ravel()

    for ax, (ycol, ytitle) in zip(axes, POPULATIONS):
        sub = df.dropna(subset=[xcol, ycol]).copy()
        x = sub[xcol].values
        y = sub[ycol].values

        # points by functional category
        for cat in CATEGORY_ORDER:
            m = sub['functional_category'] == cat
            if m.sum() == 0:
                continue
            ax.scatter(sub.loc[m, xcol], sub.loc[m, ycol],
                       c=CATEGORY_COLORS[cat], s=46, alpha=0.85,
                       edgecolors='white', linewidths=0.4, zorder=3)

        # OLS trend line + Spearman ρ
        if len(sub) >= 3:
            slope, intercept, *_ = stats.linregress(x, y)
            xs = np.linspace(x.min(), x.max(), 100)
            ax.plot(xs, slope * xs + intercept, '--', color='#555',
                    lw=1.3, alpha=0.7, zorder=2)
            rho, p = stats.spearmanr(x, y)
            stars = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
            ax.text(0.97, 0.96, f'Spearman ρ = {rho:+.3f}\np = {p:.2e}  {stars}\nn = {len(sub)}',
                    transform=ax.transAxes, ha='right', va='top', fontsize=9.5,
                    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#999', lw=0.6))

        # outliers to label: top-PBS ∪ extreme-x
        top_pbs = set(sub.nlargest(N_TOP_PBS, ycol)['gene'])
        top_x   = set(sub.nlargest(N_TOP_X,  xcol)['gene'])
        labelset = top_pbs | top_x
        texts = []
        for _, r in sub.iterrows():
            if r['gene'] in labelset:
                texts.append(ax.text(r[xcol], r[ycol], r['gene'],
                                     fontsize=8, fontweight='bold',
                                     color='#222', style='italic'))
        if texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle='-', color='#aaa', lw=0.5),
                        expand_points=(1.6, 1.8), force_text=(0.4, 0.6))

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel(f'PBS — {ytitle.split("  ")[0]}', fontsize=11)
        ax.set_title(ytitle, fontsize=11.5, fontweight='bold', loc='left', pad=6)
        ax.tick_params(labelsize=9)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)

    # shared category legend
    handles = [mlines.Line2D([], [], color=CATEGORY_COLORS[c], marker='o',
                             linestyle='', markersize=8, markeredgecolor='white',
                             label=c)
               for c in CATEGORY_ORDER if (df['functional_category'] == c).any()]
    fig.legend(handles=handles, loc='lower center', ncol=6, fontsize=9.5,
               frameon=True, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(f'Population-specific PBS vs. {xlabel} — all four populations',
                 fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.99])

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'figure_phase2_pbs_scatter_{tag}.{ext}')
        fig.savefig(path, dpi=220, bbox_inches='tight')
        print(f'  Saved → {path}')
    plt.close(fig)


for xcol, xlabel, tag in PREDICTORS:
    make_scatter_figure(xcol, xlabel, tag)
