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
    ('Breadth + KEGG',         ['tissue_breadth', 'kegg_log1p']),
    ('Breadth + betw.',        ['tissue_breadth', 'betw_sqrt']),
    ('τ + KEGG',               ['tau',            'kegg_log1p']),
    ('τ + betw.',              ['tau',            'betw_sqrt']),
    ('Breadth + KEGG + betw.', ['tissue_breadth', 'kegg_log1p', 'betw_sqrt']),
    ('τ + KEGG + betw.',       ['tau',            'kegg_log1p', 'betw_sqrt']),
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
sub = df.dropna(subset=['LOEUF', 'tissue_breadth', 'kegg_log1p']).copy()

fr_b, yr_b, slope_b, int_b, se_b, p_b, r_b = added_variable(
    sub, target='LOEUF', focal='tissue_breadth', other=['kegg_log1p'])
fr_k, yr_k, slope_k, int_k, se_k, p_k, r_k = added_variable(
    sub, target='LOEUF', focal='kegg_log1p',   other=['tissue_breadth'])

POINT_COLOR = '#E89A3C'   # uniform orange — category coloring removed
colors = np.array([POINT_COLOR] * len(sub))

# ── Figure ────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(17, 7.2))
gs  = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.15], wspace=0.62)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

def _av_panel(ax, fx, fy, slope, intercept, se, p, n, focal_label, other_label):
    ax.axhline(0, color='#999', lw=0.6, zorder=0)
    ax.axvline(0, color='#999', lw=0.6, zorder=0)
    ax.scatter(fx, fy, c=colors, s=42, edgecolor='black', linewidth=0.4,
               alpha=0.85, zorder=2)
    xs = np.linspace(fx.min(), fx.max(), 100)
    ys = intercept + slope * xs
    ax.plot(xs, ys, color='black', lw=2.0, zorder=3)
    # 95% CI band on the line
    t_crit = stats.t.ppf(0.975, df=n - 2)
    band = t_crit * se * np.sqrt(1 / n + (xs - fx.mean()) ** 2 / np.sum((fx - fx.mean()) ** 2))
    ax.fill_between(xs, ys - band, ys + band, color='black', alpha=0.12, zorder=1)
    ax.set_xlabel(f'{focal_label}  residuals\n(after removing {other_label})', fontsize=11)
    ax.set_ylabel(f'LOEUF  residuals\n(after removing {other_label})', fontsize=11)
    stars = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
    ax.text(0.04, 0.96,
            f'slope = {slope:.3f}\np = {p:.2e}  {stars}\nn = {n}',
            transform=ax.transAxes, fontsize=10.5, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#666', lw=0.6))
    ax.tick_params(labelsize=9.5)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)

_av_panel(ax1, fr_b, yr_b, slope_b, int_b, se_b, p_b, len(sub),
          focal_label='Tissue breadth', other_label='KEGG')
_av_panel(ax2, fr_k, yr_k, slope_k, int_k, se_k, p_k, len(sub),
          focal_label='KEGG  log1p(pathway count)', other_label='breadth')

ax1.set_title('A.  Breadth | KEGG controlled',
              fontsize=11.5, loc='left', pad=8, weight='bold')
ax2.set_title('B.  KEGG | breadth controlled',
              fontsize=11.5, loc='left', pad=8, weight='bold')

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

ax3.axvline(0, color='#444', lw=0.8, ls='--', zorder=0)
for _, r in forest.iterrows():
    yp = y_positions[(r['model'], r['pred'])]
    c  = PRED_COLORS[r['pred']]
    ax3.errorbar(r['beta'], yp,
                 xerr=[[r['beta'] - r['lo']], [r['hi'] - r['beta']]],
                 fmt='o', color=c, ecolor=c, elinewidth=1.6, capsize=4,
                 markersize=7, markeredgecolor='black', markeredgewidth=0.5)
    stars = '***' if r['p'] < 1e-3 else ('**' if r['p'] < 1e-2 else ('*' if r['p'] < 0.05 else 'ns'))
    ax3.text(r['hi'] + 0.015, yp, stars, fontsize=10.5, va='center')

# Predictor legend (build manually so it includes all unique predictors used)
import matplotlib.lines as mlines
handles = []
seen = set()
for _, preds in MODELS:
    for nm in preds:
        if nm in seen:
            continue
        seen.add(nm)
        handles.append(mlines.Line2D([], [], color=PRED_COLORS[nm], marker='o',
                                     linestyle='', markersize=7,
                                     markeredgecolor='black', markeredgewidth=0.5,
                                     label=PRED_DISPLAY[nm]))
ax3.legend(handles=handles, fontsize=9, loc='upper left',
           bbox_to_anchor=(1.02, 1.0), frameon=True, borderaxespad=0)

ax3.set_yticks(yticks)
ax3.set_yticklabels(yticklabels, fontsize=10)
ax3.invert_yaxis()
ax3.set_xlabel('Standardized β (95% CI)', fontsize=11)
ax3.set_title('C.  Joint-model β across specifications',
              fontsize=11.5, loc='left', pad=8, weight='bold')
ax3.tick_params(labelsize=9.5)
for spine in ('top', 'right'):
    ax3.spines[spine].set_visible(False)

fig.suptitle('Tissue expression breadth and network connectivity independently predict LOEUF',
             fontsize=13.5, weight='bold', y=1.02)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_joint_regression.{ext}')
    fig.savefig(path, dpi=220, bbox_inches='tight')
    print(f'  Saved → {path}')
