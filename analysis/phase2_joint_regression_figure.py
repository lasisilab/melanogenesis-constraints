"""
phase2_joint_regression_figure.py

Poster figure for the joint regression result.

Three panels in one wide figure:
  1. Added-variable plot: residual(LOEUF | KEGG) vs residual(breadth | KEGG)
  2. Added-variable plot: residual(LOEUF | breadth) vs residual(KEGG | breadth)
  3. Forest plot of standardized β ± 95% CI for breadth, KEGG, betweenness
     across three models (breadth+KEGG, τ+KEGG, all three).

Outputs:
  output/figure_phase2_joint_regression.png
  output/figure_phase2_joint_regression.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data',   'kegg_pathway_counts.csv')
BREADTH_CSV = os.path.join(PROJECT_DIR, 'data',   'gtex_tissue_breadth.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',
    'Developmental/NC':          '#E8907E',
    'Generic signaling':         '#F5C242',
    'Cytokines/growth factors':  '#4878CF',
    'Apoptosis/cell death':      '#6BAD6B',
    'Other':                     '#B0B0B0',
}

# ── Load + merge ──────────────────────────────────────────────────────────
master  = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
kegg    = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()
breadth = pd.read_csv(BREADTH_CSV)[['gene', 'tissue_breadth']]
breadth['gene'] = breadth['gene'].str.upper()

df = (master[['gene', 'functional_category', 'LOEUF',
              'betweenness_centrality', 'tau']]
      .merge(kegg,    on='gene', how='left')
      .merge(breadth, on='gene', how='left'))

df['kegg_log1p'] = np.log1p(df['kegg_pathway_count'])
df['betw_sqrt']  = np.sqrt(df['betweenness_centrality'].clip(lower=0))

def _z(s):
    s = pd.to_numeric(s, errors='coerce')
    return (s - s.mean()) / s.std(ddof=0)

# ── Fit the three models we will plot ─────────────────────────────────────
def fit_zscored(df_in, predictors):
    sub = df_in.dropna(subset=['LOEUF'] + list(predictors)).copy()
    X = pd.DataFrame({p: _z(sub[p]) for p in predictors}, index=sub.index)
    Xc = sm.add_constant(X.values)
    res = sm.OLS(sub['LOEUF'].values, Xc).fit()
    return res, sub, X

MODELS = [
    ('Breadth + KEGG',  ['tissue_breadth', 'kegg_log1p']),
    ('Breadth + betw.', ['tissue_breadth', 'betw_sqrt']),
    ('τ + KEGG',        ['tau',            'kegg_log1p']),
    ('τ + betw.',       ['tau',            'betw_sqrt']),
]

fits = []
for label, preds in MODELS:
    res, sub, X = fit_zscored(df, preds)
    fits.append((label, preds, res, sub, X))

# Use the primary model (breadth + KEGG) for the added-variable scatters
primary_label, primary_preds, primary_res, primary_sub, primary_X = fits[0]

# ── Added-variable plot helper ────────────────────────────────────────────
def added_variable(df_sub, target, focal, other):
    """Return residuals for the AV plot: y_resid vs focal_resid, both after
    regressing on `other`. Slope of y_resid on focal_resid equals the focal
    coefficient in the multiple regression."""
    Xo = sm.add_constant(df_sub[other].values)
    y_resid     = df_sub[target].values - sm.OLS(df_sub[target].values, Xo).fit().fittedvalues
    focal_resid = df_sub[focal].values  - sm.OLS(df_sub[focal].values,  Xo).fit().fittedvalues
    # OLS through residuals (intercept ≈ 0 by construction)
    slope, intercept, r, p, se = stats.linregress(focal_resid, y_resid)
    return focal_resid, y_resid, slope, intercept, se, p, r

# Use ORIGINAL-SCALE predictors here (not z-scored) so axes have real units.
# AV plots for breadth + KEGG model
sub_bk = df.dropna(subset=['LOEUF', 'tissue_breadth', 'kegg_log1p']).copy()
fr_b,  yr_b,  slope_b,  int_b,  se_b,  p_b,  r_b  = added_variable(
    sub_bk, target='LOEUF', focal='tissue_breadth', other=['kegg_log1p'])
fr_k,  yr_k,  slope_k,  int_k,  se_k,  p_k,  r_k  = added_variable(
    sub_bk, target='LOEUF', focal='kegg_log1p',     other=['tissue_breadth'])

# AV plots for τ + KEGG model
sub_tk = df.dropna(subset=['LOEUF', 'tau', 'kegg_log1p']).copy()
fr_t,  yr_t,  slope_t,  int_t,  se_t,  p_t,  r_t  = added_variable(
    sub_tk, target='LOEUF', focal='tau',         other=['kegg_log1p'])
fr_k2, yr_k2, slope_k2, int_k2, se_k2, p_k2, r_k2 = added_variable(
    sub_tk, target='LOEUF', focal='kegg_log1p',  other=['tau'])

POINT_COLOR = '#E89A3C'   # uniform orange — category coloring removed
colors_bk = np.array([POINT_COLOR] * len(sub_bk))
colors_tk = np.array([POINT_COLOR] * len(sub_tk))

# ── Figure (wide & short poster size, 5 panels) ────────────────────────────
fig = plt.figure(figsize=(15.5, 2.2))
gs  = fig.add_gridspec(1, 5, width_ratios=[1.0, 1.0, 1.0, 1.0, 1.2], wspace=0.52)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])
ax5 = fig.add_subplot(gs[0, 4])

def _av_panel(ax, fx, fy, colors, slope, intercept, se, p, n, focal_label, other_label):
    ax.axhline(0, color='#999', lw=0.6, zorder=0)
    ax.axvline(0, color='#999', lw=0.6, zorder=0)
    ax.scatter(fx, fy, c=colors, s=14, edgecolor='black', linewidth=0.3,
               alpha=0.85, zorder=2)
    xs = np.linspace(fx.min(), fx.max(), 100)
    ys = intercept + slope * xs
    ax.plot(xs, ys, color='black', lw=1.3, zorder=3)
    # 95% CI band on the line
    t_crit = stats.t.ppf(0.975, df=n - 2)
    band = t_crit * se * np.sqrt(1 / n + (xs - fx.mean()) ** 2 / np.sum((fx - fx.mean()) ** 2))
    ax.fill_between(xs, ys - band, ys + band, color='black', alpha=0.12, zorder=1)
    ax.set_xlabel(f'{focal_label} resid.\n(remove {other_label})', fontsize=8)
    ax.set_ylabel(f'LOEUF resid.\n(remove {other_label})', fontsize=8)
    stars = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
    ax.text(0.04, 0.96,
            f'slope={slope:.3f}\np={p:.1e} {stars}\nn={n}',
            transform=ax.transAxes, fontsize=7, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#666', lw=0.5))
    ax.tick_params(labelsize=7, pad=1)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)

_av_panel(ax1, fr_b, yr_b, colors_bk, slope_b, int_b, se_b, p_b, len(sub_bk),
          focal_label='Tissue breadth', other_label='KEGG')
_av_panel(ax2, fr_k, yr_k, colors_bk, slope_k, int_k, se_k, p_k, len(sub_bk),
          focal_label='KEGG log1p', other_label='breadth')
_av_panel(ax3, fr_t, yr_t, colors_tk, slope_t, int_t, se_t, p_t, len(sub_tk),
          focal_label='τ (tissue spec.)', other_label='KEGG')
_av_panel(ax4, fr_k2, yr_k2, colors_tk, slope_k2, int_k2, se_k2, p_k2, len(sub_tk),
          focal_label='KEGG log1p', other_label='τ')

ax1.set_title('A.  Breadth | KEGG', fontsize=9.5, loc='left', pad=4, weight='bold')
ax2.set_title('B.  KEGG | breadth', fontsize=9.5, loc='left', pad=4, weight='bold')
ax3.set_title('C.  τ | KEGG',       fontsize=9.5, loc='left', pad=4, weight='bold')
ax4.set_title('D.  KEGG | τ',       fontsize=9.5, loc='left', pad=4, weight='bold')

# ── Forest plot of standardized β ± 95% CI ───────────────────────────────
PRED_DISPLAY = {
    'tissue_breadth': 'Tissue breadth',
    'tau':            'τ (tissue specificity)',
    'kegg_log1p':     'KEGG  log1p(pathway count)',
    'betw_sqrt':      'Betweenness  (sqrt)',
}
PRED_COLORS = {
    'tissue_breadth': '#1f77b4',
    'tau':            '#d62728',
    'kegg_log1p':     '#2ca02c',
    'betw_sqrt':      '#9467bd',
}

# Collect rows: one per (model, predictor)
rows = []
for label, preds, res, sub_m, X_m in fits:
    names = ['const'] + list(preds)
    ci = res.conf_int(0.05)
    for i, nm in enumerate(names):
        if nm == 'const':
            continue
        rows.append({
            'model':    label,
            'pred':     nm,
            'beta':     res.params[i],
            'lo':       ci[i, 0],
            'hi':       ci[i, 1],
            'p':        res.pvalues[i],
        })

forest = pd.DataFrame(rows)

# Layout: y position grouped by model, sub-positions per predictor
model_labels = [m[0] for m in MODELS]
y_positions = {}
y = 0
yticks, yticklabels = [], []
for label, preds in MODELS:
    n_preds = len(preds)
    # center the group around an integer model index
    base = y
    for i, nm in enumerate(preds):
        y_positions[(label, nm)] = base + i * 0.7
    yticks.append(base + (n_preds - 1) * 0.35)
    yticklabels.append(label)
    y = base + n_preds * 0.7 + 0.9   # gap between models

ax5.axvline(0, color='#444', lw=0.7, ls='--', zorder=0)
for _, r in forest.iterrows():
    yp = y_positions[(r['model'], r['pred'])]
    c  = PRED_COLORS[r['pred']]
    ax5.errorbar(r['beta'], yp,
                 xerr=[[r['beta'] - r['lo']], [r['hi'] - r['beta']]],
                 fmt='o', color=c, ecolor=c, elinewidth=1.2, capsize=3,
                 markersize=6, markeredgecolor='black', markeredgewidth=0.4)
    stars = '***' if r['p'] < 1e-3 else ('**' if r['p'] < 1e-2 else ('*' if r['p'] < 0.05 else 'ns'))
    ax5.text(r['hi'] + 0.01, yp, stars, fontsize=8, va='center')

ax5.set_yticks(yticks)
ax5.set_yticklabels(yticklabels, fontsize=8)
ax5.invert_yaxis()
ax5.set_xlabel('Standardized β (95% CI)', fontsize=8)
ax5.set_title('E.  Joint-model β', fontsize=9.5, loc='left', pad=4, weight='bold')
ax5.tick_params(labelsize=7, pad=1)
for spine in ('top', 'right'):
    ax5.spines[spine].set_visible(False)

fig.suptitle('Breadth & network connectivity independently predict LOEUF',
             fontsize=10.5, weight='bold', y=1.03)
fig.tight_layout(rect=[0, 0, 1, 0.95])

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_joint_regression.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f'  Saved → {path}')
