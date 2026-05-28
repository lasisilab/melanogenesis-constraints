"""
phase2_pbs_regression.py

Parallel OLS regressions predicting population-specific PBS from network
centrality and tissue-expression metrics — the selection-side mirror of the
LOEUF joint regression (phase2_joint_regression_figure.py).

Two outcomes, fit separately:
  pbs1_african    — PBS-1, African target, S. Asian outgroup
  pbs3_melanesian — PBS-3, Melanesian target, S. Asian outgroup

Four bivariate specifications (same as the LOEUF regression), predictors z-scored:
  1. tissue_breadth + log1p(KEGG)
  2. tissue_breadth + sqrt(betweenness)
  3. tau            + log1p(KEGG)
  4. tau            + sqrt(betweenness)

Plus an optional 5th "orthogonal-axes" specification per outcome:
  5. tau + log1p(KEGG) + LOEUF
     If tau and KEGG still predict PBS after controlling for LOEUF, selection
     operates independently of constraint.

Right-skew of PBS is checked; if substantial, a log1p(PBS) sensitivity pass is
run and directional changes are noted. Untransformed PBS is the primary model.

Outputs:
  output/phase2_pbs_regression.txt        — full coefficient table
  output/figure_phase2_pbs_regression.png — forest plots (both outcomes) +
  output/figure_phase2_pbs_regression.pdf   added-variable plot per outcome
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data',   'kegg_pathway_counts.csv')
BREADTH_CSV = os.path.join(PROJECT_DIR, 'data',   'gtex_tissue_breadth.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

# ── Load + merge (mirrors phase2_joint_regression_figure.py) ───────────────
master  = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
kegg    = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()
breadth = pd.read_csv(BREADTH_CSV)[['gene', 'tissue_breadth']]
breadth['gene'] = breadth['gene'].str.upper()
pbs_full = pd.read_csv(PBS_CSV)
pbs_full['gene'] = pbs_full['gene'].str.upper()

# European and South Asian PBS are not pre-computed in pbs_per_gene.csv (only
# the African- and Melanesian-target scans are). There is no European–South
# Asian FST, so both must be built from the only complete trios available,
# using African + Melanesian as the reference pair:
#   PBS_european   = (T(afr,eur) + T(eur,mel) − T(afr,mel)) / 2
#   PBS_southasian = (T(afr,sas) + T(sas,mel) − T(afr,mel)) / 2
# (different reference structure than PBS-1/3, which used a S. Asian outgroup,
#  so interpret as European-/South-Asian-specific drift vs the Afr–Mel split.)
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

df = (master[['gene', 'functional_category', 'LOEUF',
              'betweenness_centrality', 'tau']]
      .merge(kegg,    on='gene', how='left')
      .merge(breadth, on='gene', how='left')
      .merge(pbs,     on='gene', how='inner'))

df['kegg_log1p'] = np.log1p(df['kegg_pathway_count'])
df['betw_sqrt']  = np.sqrt(df['betweenness_centrality'].clip(lower=0))


def _z(s):
    s = pd.to_numeric(s, errors='coerce')
    return (s - s.mean()) / s.std(ddof=0)


OUTCOMES = [
    ('pbs1_african',    'PBS-1  African (S. Asian outgroup)'),
    ('pbs3_melanesian', 'PBS-3  Melanesian (S. Asian outgroup)'),
    ('pbs_european',    'PBS  European (Afr–Mel reference)'),
    ('pbs_southasian',  'PBS  South Asian (Afr–Mel reference)'),
]

# Which outcomes go in which figure (each figure is 2 columns wide)
FIGURE_GROUPS = [
    ('figure_phase2_pbs_regression',
     'Network centrality and tissue specificity as predictors of '
     'population-specific PBS',
     ['pbs1_african', 'pbs3_melanesian']),
    ('figure_phase2_pbs_regression_eur_sas',
     'Network centrality and tissue specificity vs. European / South Asian PBS',
     ['pbs_european', 'pbs_southasian']),
]

MODELS = [
    ('Breadth + KEGG',         ['tissue_breadth', 'kegg_log1p']),
    ('Breadth + betw.',        ['tissue_breadth', 'betw_sqrt']),
    ('τ + KEGG',               ['tau',            'kegg_log1p']),
    ('τ + betw.',              ['tau',            'betw_sqrt']),
    ('τ + KEGG + LOEUF',       ['tau',            'kegg_log1p', 'LOEUF']),
]

PRED_DISPLAY = {
    'tissue_breadth': 'Tissue breadth',
    'tau':            'τ (tissue specificity)',
    'kegg_log1p':     'KEGG  log1p(pathway count)',
    'betw_sqrt':      'Betweenness  (sqrt)',
    'LOEUF':          'LOEUF  (constraint covariate)',
}
PRED_COLORS = {
    'tissue_breadth': '#1f77b4',   # blue
    'tau':            '#d62728',   # red
    'kegg_log1p':     '#2ca02c',   # green
    'betw_sqrt':      '#9467bd',   # purple
    'LOEUF':          '#7f7f7f',   # gray (covariate)
}


def fit_zscored(df_in, outcome, predictors, transform=None):
    """Fit OLS of (optionally transformed) outcome on z-scored predictors."""
    sub = df_in.dropna(subset=[outcome] + list(predictors)).copy()
    y = sub[outcome].values.astype(float)
    if transform == 'log1p':
        y = np.log1p(y)
    X = pd.DataFrame({p: _z(sub[p]) for p in predictors}, index=sub.index)
    Xc = sm.add_constant(X.values)
    res = sm.OLS(y, Xc).fit()
    return res, sub, X


def vifs(X):
    """VIF per predictor (with intercept column added)."""
    Xv = np.column_stack([np.ones(len(X)), X.values])
    return {col: variance_inflation_factor(Xv, i + 1)
            for i, col in enumerate(X.columns)}


def added_variable(df_sub, target, focal, others):
    """Residual-vs-residual AV plot data (target & focal each residualized on
    `others`). Slope == focal's multiple-regression coefficient."""
    Xo = sm.add_constant(df_sub[others].values)
    y_resid = df_sub[target].values - sm.OLS(df_sub[target].values, Xo).fit().fittedvalues
    f_resid = df_sub[focal].values  - sm.OLS(df_sub[focal].values,  Xo).fit().fittedvalues
    slope, intercept, r, p, se = stats.linregress(f_resid, y_resid)
    return f_resid, y_resid, slope, intercept, se, p


# ── Fit everything; collect coefficient rows + text table ───────────────────
report = ["phase2_pbs_regression.py — PBS predicted by network & tissue metrics",
          "=" * 74,
          "Predictors z-scored; OLS; β are standardized. VIF<5 = acceptable "
          "collinearity.",
          ""]

fits_by_outcome = {}   # outcome -> list of (label, preds, res, sub, X)
skew_by_outcome = {}

for outcome, out_label in OUTCOMES:
    vals = df[outcome].dropna().values
    sk = stats.skew(vals)
    skew_by_outcome[outcome] = sk
    report.append(f"\n{'#'*74}\n# OUTCOME: {out_label}\n{'#'*74}")
    report.append(f"  n (non-null) = {len(vals)};  skewness = {sk:+.2f} "
                  f"({'right-skewed — log1p sensitivity run below' if sk > 1 else 'approx symmetric'})")

    fits = []
    for label, preds in MODELS:
        res, sub, X = fit_zscored(df, outcome, preds)
        fits.append((label, preds, res, sub, X))
        names = ['const'] + list(preds)
        ci = res.conf_int(0.05)
        vif = vifs(X)
        report.append(f"\n  ── {label}  (n={len(sub)}) ──")
        report.append(f"     R²={res.rsquared:.3f}  adjR²={res.rsquared_adj:.3f}  "
                      f"F p={res.f_pvalue:.3e}")
        report.append(f"     {'predictor':<28}{'β':>8}{'SE':>8}{'t':>8}"
                      f"{'p':>11}{'95% CI':>20}{'VIF':>7}")
        for i, nm in enumerate(names):
            if nm == 'const':
                continue
            report.append(
                f"     {PRED_DISPLAY.get(nm, nm):<28}"
                f"{res.params[i]:>8.3f}{res.bse[i]:>8.3f}{res.tvalues[i]:>8.2f}"
                f"{res.pvalues[i]:>11.2e}"
                f"   [{ci[i,0]:+.3f}, {ci[i,1]:+.3f}]"
                f"{vif[nm]:>7.2f}")
    fits_by_outcome[outcome] = fits

    # log1p sensitivity for the 4 bivariate specs if right-skewed
    if sk > 1:
        report.append(f"\n  ── log1p(PBS) sensitivity (directional check) ──")
        report.append(f"     {'spec':<20}{'predictor':<26}"
                      f"{'β(raw)':>9}{'p(raw)':>11}{'β(log)':>9}{'p(log)':>11}  flag")
        for label, preds in MODELS[:4]:
            res_raw, _, _ = fit_zscored(df, outcome, preds)
            res_log, _, _ = fit_zscored(df, outcome, preds, transform='log1p')
            for i, nm in enumerate(preds):
                br, pr = res_raw.params[i+1], res_raw.pvalues[i+1]
                bl, pl = res_log.params[i+1], res_log.pvalues[i+1]
                sign_flip = (np.sign(br) != np.sign(bl))
                sig_flip  = ((pr < 0.05) != (pl < 0.05))
                flag = ('SIGN-FLIP' if sign_flip else
                        ('sig-change' if sig_flip else 'consistent'))
                report.append(f"     {label:<20}{PRED_DISPLAY.get(nm, nm):<26}"
                              f"{br:>9.3f}{pr:>11.2e}{bl:>9.3f}{pl:>11.2e}  {flag}")

with open(os.path.join(OUT_DIR, 'phase2_pbs_regression.txt'), 'w') as f:
    f.write("\n".join(report) + "\n")
print("  Wrote output/phase2_pbs_regression.txt")


# ── Pick the most-significant bivariate spec per outcome (for AV plot) ──────
def best_spec(fits):
    """Bivariate spec (of the first 4) minimizing the larger of its two
    predictor p-values; returns (label, preds, res, sub, X)."""
    best, best_maxp = None, np.inf
    for f in fits:
        label, preds, res, sub, X = f
        if len(preds) != 2:
            continue
        maxp = max(res.pvalues[1], res.pvalues[2])
        if maxp < best_maxp:
            best_maxp, best = maxp, f
    return best


# ── Figure builder: 2×2 — top row AV plots, bottom row forest plots ────────
OUT_LABELS = dict(OUTCOMES)
PANEL_LETTERS = ['A', 'B', 'C', 'D']


def make_figure(fname, suptitle, outcomes):
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.25], hspace=0.32, wspace=0.32)

    for col, outcome in enumerate(outcomes):
        out_label = OUT_LABELS[outcome]
        fits = fits_by_outcome[outcome]

        # ---- Top row: added-variable plot for the best spec's stronger predictor
        ax_av = fig.add_subplot(gs[0, col])
        label, preds, res, sub, X = best_spec(fits)
        p_by_pred = {preds[i]: res.pvalues[i+1] for i in range(len(preds))}
        focal = min(p_by_pred, key=p_by_pred.get)
        others = [p for p in preds if p != focal]
        fr, yr, slope, intc, se, p = added_variable(sub, outcome, focal, others)

        ax_av.axhline(0, color='#999', lw=0.6, zorder=0)
        ax_av.axvline(0, color='#999', lw=0.6, zorder=0)
        ax_av.scatter(fr, yr, c='#E89A3C', s=40, edgecolor='black', linewidth=0.4,
                      alpha=0.85, zorder=2)
        xs = np.linspace(fr.min(), fr.max(), 100)
        ys = intc + slope * xs
        ax_av.plot(xs, ys, color='black', lw=2.0, zorder=3)
        n = len(sub)
        t_crit = stats.t.ppf(0.975, df=n - 2)
        band = t_crit * se * np.sqrt(1/n + (xs - fr.mean())**2 / np.sum((fr - fr.mean())**2))
        ax_av.fill_between(xs, ys - band, ys + band, color='black', alpha=0.12, zorder=1)
        stars = '***' if p < 1e-3 else ('**' if p < 1e-2 else ('*' if p < 0.05 else 'ns'))
        ax_av.text(0.04, 0.96,
                   f'{PRED_DISPLAY[focal]}\nspec: {label}\nslope={slope:.3f}  p={p:.2e} {stars}\nn={n}',
                   transform=ax_av.transAxes, fontsize=10, va='top', ha='left',
                   bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#666', lw=0.6))
        ax_av.set_xlabel(f'{PRED_DISPLAY[focal]} residuals\n(after removing {", ".join(PRED_DISPLAY[o] for o in others)})',
                         fontsize=10)
        ax_av.set_ylabel(f'{out_label} residuals', fontsize=10)
        ax_av.set_title(f'{PANEL_LETTERS[col]}.  {out_label} — added-variable',
                        fontsize=11.5, loc='left', pad=8, weight='bold')
        ax_av.tick_params(labelsize=9)
        for sp in ('top', 'right'):
            ax_av.spines[sp].set_visible(False)

        # ---- Bottom row: forest plot of standardized β ± 95% CI
        ax_f = fig.add_subplot(gs[1, col])
        rows = []
        for label_m, preds_m, res_m, sub_m, X_m in fits:
            ci = res_m.conf_int(0.05)
            for i, nm in enumerate(preds_m):
                rows.append({'model': label_m, 'pred': nm,
                             'beta': res_m.params[i+1],
                             'lo': ci[i+1, 0], 'hi': ci[i+1, 1],
                             'p': res_m.pvalues[i+1]})
        forest = pd.DataFrame(rows)

        y_positions = {}
        y = 0
        yticks, yticklabels = [], []
        for label_m, preds_m in MODELS:
            base = y
            for i, nm in enumerate(preds_m):
                y_positions[(label_m, nm)] = base + i * 0.6
            yticks.append(base + (len(preds_m) - 1) * 0.3)
            yticklabels.append(label_m)
            y = base + len(preds_m) * 0.6 + 0.8

        ax_f.axvline(0, color='#444', lw=0.8, ls='--', zorder=0)
        for _, r in forest.iterrows():
            yp = y_positions[(r['model'], r['pred'])]
            c = PRED_COLORS[r['pred']]
            ax_f.errorbar(r['beta'], yp,
                          xerr=[[r['beta'] - r['lo']], [r['hi'] - r['beta']]],
                          fmt='o', color=c, ecolor=c, elinewidth=1.6, capsize=4,
                          markersize=7, markeredgecolor='black', markeredgewidth=0.5)
            stars = '***' if r['p'] < 1e-3 else ('**' if r['p'] < 1e-2 else ('*' if r['p'] < 0.05 else 'ns'))
            ax_f.text(r['hi'] + 0.01, yp, stars, fontsize=9.5, va='center')

        ax_f.set_yticks(yticks)
        ax_f.set_yticklabels(yticklabels, fontsize=9.5)
        ax_f.invert_yaxis()
        ax_f.set_xlabel('Standardized β (95% CI)', fontsize=10.5)
        ax_f.set_title(f'{PANEL_LETTERS[col+2]}.  {out_label} — coefficients',
                       fontsize=11.5, loc='left', pad=8, weight='bold')
        ax_f.tick_params(labelsize=9)
        for sp in ('top', 'right'):
            ax_f.spines[sp].set_visible(False)

    # Shared predictor legend (bottom)
    seen, handles = set(), []
    for _, preds in MODELS:
        for nm in preds:
            if nm in seen:
                continue
            seen.add(nm)
            handles.append(mlines.Line2D([], [], color=PRED_COLORS[nm], marker='o',
                                         linestyle='', markersize=7,
                                         markeredgecolor='black', markeredgewidth=0.5,
                                         label=PRED_DISPLAY[nm]))
    fig.legend(handles=handles, fontsize=9.5, loc='lower center', ncol=5,
               frameon=True, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(suptitle, fontsize=14, weight='bold', y=0.995)

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'{fname}.{ext}')
        fig.savefig(path, dpi=220, bbox_inches='tight')
        print(f'  Saved → {path}')
    plt.close(fig)


for fname, suptitle, outcomes in FIGURE_GROUPS:
    make_figure(fname, suptitle, outcomes)
