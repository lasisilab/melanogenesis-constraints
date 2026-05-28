"""
phase2_pbs_within_network_combined.py

Within-network PBS association on the COMBINED-dataset network PBS
(pbs_per_gene.csv: gnomAD HGDP+1KGP for African/S.Asian/European, HGDP+SGDP for
Melanesian). This is the better-powered dataset than the SGDP-only genome-wide
scan, but it is network-scale only — there is no genome-wide background — so the
outcome is the WITHIN-NETWORK PBS percentile (rank among the network genes,
0–100), not a genome-wide percentile.

Mirrors phase2_pbs_within_network.py (same predictors, covariates, tests, and
2×2 forest layout) so the two are directly comparable:
  - that script: SGDP-only genome-wide percentile  → null everywhere
  - this script: combined-dataset within-network percentile → better powered

NOTE on interpretation: Spearman ρ here is identical to ρ on raw combined PBS
(rank-invariant). Using within-network percentile only sets a common 0–100 scale
across branches; it does NOT make the combined-dataset signal more robust than
the SGDP null — it reports the dataset where the signal exists. Caption both.

Outputs (output/):
  phase2_pbs_within_network_combined.{txt,csv}
  figure_phase2_pbs_within_network_combined.{png,pdf}
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
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_counts.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')

N_PERM = 10000
RNG = np.random.default_rng(0)
PAD = 20000   # ±10 kb flanking added by the gene-window BED


# ── Load combined-dataset network PBS ──────────────────────────────────────
pbs = pd.read_csv(PBS_CSV)
pbs['gene'] = pbs['gene'].str.upper()


def _t(fst):
    fst = np.clip(pd.to_numeric(fst, errors='coerce').fillna(0).clip(lower=0), 0, 0.9999)
    return -np.log(1.0 - fst)


t_afr_eur = _t(pbs['fst_african_european'])
t_afr_mel = _t(pbs['fst_african_melanesian'])
t_afr_sas = _t(pbs['fst_african_southasian'])
t_eur_mel = _t(pbs['fst_european_melanesian'])
t_mel_sas = _t(pbs['fst_melanesian_southasian'])

pbs['pbs_afr'] = pbs['pbs1_african'].clip(lower=0)
pbs['pbs_mel'] = pbs['pbs3_melanesian'].clip(lower=0)
pbs['pbs_eur'] = ((t_afr_eur + t_eur_mel - t_afr_mel) / 2).clip(lower=0)
pbs['pbs_sas'] = ((t_afr_sas + t_mel_sas - t_afr_mel) / 2).clip(lower=0)

BRANCHES = [
    ('pbs_afr', 'AFR', 'African (PBS-1, S.Asian outgroup)'),
    ('pbs_mel', 'MEL', 'Melanesian (PBS-3, S.Asian outgroup)'),
    ('pbs_eur', 'EUR', 'European (Afr–Mel reference)'),
    ('pbs_sas', 'SAS', 'South Asian (Afr–Mel reference)'),
]


# ── Merge network features ─────────────────────────────────────────────────
master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
kegg = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()
feat = master[['gene', 'functional_category', 'tau', 'betweenness_centrality',
               'LOEUF', 'tissue_breadth']].merge(kegg, on='gene', how='left')

net = pbs.merge(feat, on='gene', how='inner')
print(f"{len(net)} network genes with combined-dataset PBS")


# ── Within-network percentile (rank among the network genes) per branch ────
for col, code, _ in BRANCHES:
    pct = np.full(len(net), np.nan)
    valid = net[col].notna().values & (net['n_sites_shared'] > 0).values
    vals = net.loc[valid, col].values
    ranks = stats.rankdata(vals, method='average')
    pct[valid] = 100.0 * (ranks - 1) / (len(vals) - 1)
    net[f'{code}_pctile'] = pct


# ── Predictors / covariates / transforms ───────────────────────────────────
net['kegg_log1p']  = np.log1p(net['kegg_pathway_count'])
net['nsnps_log1p'] = np.log1p(net['n_sites_shared'])
net['body_length'] = (net['end'] - net['start'] - PAD).clip(lower=1)
net['glen_log10']  = np.log10(net['body_length'])

PREDICTORS = [
    ('tau',                    'τ (tissue specificity)'),
    ('betweenness_centrality', 'Betweenness centrality'),
    ('kegg_log1p',             'KEGG  log1p(pathway count)'),
    ('LOEUF',                  'LOEUF (constraint)'),
    ('tissue_breadth',         'Tissue breadth'),
]
COVARS = ['nsnps_log1p', 'glen_log10']
PRED_COLORS = {
    'tau': '#d62728', 'betweenness_centrality': '#9467bd',
    'kegg_log1p': '#2ca02c', 'LOEUF': '#7f7f7f', 'tissue_breadth': '#1f77b4',
}
PRED_DISPLAY = dict(PREDICTORS)


def _z(s):
    s = pd.to_numeric(s, errors='coerce')
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)


def _dir(beta, pred):
    if pred == 'LOEUF':
        return 'less-constrained genes higher' if beta > 0 else 'more-constrained genes higher'
    if pred == 'tau':
        return 'tissue-specific genes higher' if beta > 0 else 'broadly-expressed genes higher'
    if pred == 'betweenness_centrality':
        return 'central genes higher' if beta > 0 else 'peripheral genes higher'
    if pred == 'kegg_log1p':
        return 'more-pathway genes higher' if beta > 0 else 'fewer-pathway genes higher'
    if pred == 'tissue_breadth':
        return 'broadly-expressed genes higher' if beta > 0 else 'tissue-restricted genes higher'
    return ''


# ── Tests per branch × predictor ───────────────────────────────────────────
rows = []
detail = ["Within-network PBS association — COMBINED dataset (HGDP+1KGP+SGDP)",
          "=" * 78,
          "Outcome: WITHIN-network PBS percentile (rank among network genes, 0–100).",
          "Source: pbs_per_gene.csv (network-scale; no genome-wide background).",
          "Adjusted models: pctile ~ z(predictor) + z(log1p n_sites) + z(log10 body_len).",
          f"Permutation: {N_PERM:,} shuffles of the percentile among network genes.",
          "NOTE: ρ identical to ρ on raw combined PBS (rank-invariant). This reports",
          "the better-powered dataset; it does NOT resolve the SGDP-vs-combined",
          "robustness gap — the SGDP genome-wide version of this analysis is null.",
          ""]

for col, code, blabel in BRANCHES:
    ycol = f'{code}_pctile'
    detail.append(f"\n{'#'*78}\n# BRANCH: {code} — {blabel}\n{'#'*78}")
    sub_all = net.dropna(subset=[ycol]).copy()
    detail.append(f"  network genes with {code} percentile: {len(sub_all)}")

    for pred, pdisp in PREDICTORS:
        sub = sub_all.dropna(subset=[pred] + COVARS).copy()
        n = len(sub)
        if n < 10:
            detail.append(f"  {pdisp}: n={n} too few, skipped")
            continue
        y = sub[ycol].values
        rho, p_sp = stats.spearmanr(sub[pred].values, y)

        Xcols = [pred] + COVARS
        X = pd.DataFrame({c: _z(sub[c]) for c in Xcols}, index=sub.index)
        Xc = sm.add_constant(X.values)
        res = sm.OLS(y, Xc).fit()
        beta, se, p_b = res.params[1], res.bse[1], res.pvalues[1]
        ci = res.conf_int(0.05)[1]
        Xv = np.column_stack([np.ones(len(X)), X.values])
        vif = variance_inflation_factor(Xv, 1)

        null = np.empty(N_PERM)
        for i in range(N_PERM):
            null[i] = sm.OLS(RNG.permutation(y), Xc).fit().params[1]
        p_perm = (np.sum(np.abs(null) >= abs(beta)) + 1) / (N_PERM + 1)

        rows.append({'branch': code, 'predictor': pred, 'n': n,
                     'spearman_rho': rho, 'spearman_p': p_sp,
                     'adj_beta': beta, 'adj_se': se, 'adj_ci_lo': ci[0],
                     'adj_ci_hi': ci[1], 'adj_p': p_b, 'perm_p': p_perm,
                     'vif': vif, 'direction': _dir(beta, pred)})
        stars = ('***' if p_perm < 1e-3 else '**' if p_perm < 1e-2
                 else '*' if p_perm < 0.05 else 'ns')
        detail.append(
            f"  {pdisp:<28} ρ={rho:+.3f} (p={p_sp:.2e}) | "
            f"adjβ={beta:+.2f} [{ci[0]:+.2f},{ci[1]:+.2f}] p={p_b:.2e} | "
            f"perm p={p_perm:.4f} {stars} | VIF={vif:.2f} | {_dir(beta, pred)}")

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, 'phase2_pbs_within_network_combined.csv'), index=False)

detail.append(f"\n\n{'='*78}\nSUMMARY TABLE  (β on 0–100 within-network percentile scale)\n{'='*78}")
detail.append(f"{'branch':<7}{'predictor':<26}{'rho':>7}{'adjβ':>8}{'adjp':>11}{'permp':>9}  direction")
for _, r in summary.iterrows():
    detail.append(f"{r['branch']:<7}{PRED_DISPLAY[r['predictor']][:25]:<26}"
                  f"{r['spearman_rho']:>+7.3f}{r['adj_beta']:>+8.2f}"
                  f"{r['adj_p']:>11.2e}{r['perm_p']:>9.4f}  {r['direction']}")

detail += [
    "", "-" * 78, "CAVEATS:",
    "  • Combined dataset: gnomAD HGDP+1KGP (AFR/SAS/EUR) + HGDP+SGDP (MEL, n=47).",
    "  • Network-scale only → outcome is WITHIN-network percentile, not genome-wide.",
    "  • The SGDP-only genome-wide version of this analysis is NULL; this signal",
    "    does not replicate at SGDP sample sizes. Report as dataset-dependent.",
    "  • EUR/SAS use Afr–Mel reference (no eur–sas FST). EUR = depigmentation contrast.",
    "  • PBS is a frequency-differentiation statistic, NOT proof of selection.",
]
with open(os.path.join(OUT_DIR, 'phase2_pbs_within_network_combined.txt'), 'w') as f:
    f.write("\n".join(detail) + "\n")
print("Wrote phase2_pbs_within_network_combined.{txt,csv}")


# ── Figure: 2×2 forest, one panel per branch ───────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
axes = axes.ravel()
pred_order = [p for p, _ in PREDICTORS]

for ax, (col, code, blabel) in zip(axes, BRANCHES):
    bdf = summary[summary['branch'] == code].set_index('predictor')
    yvals = list(range(len(pred_order)))[::-1]
    ax.axvline(0, color='#444', lw=0.8, ls='--', zorder=0)
    for yp, pred in zip(yvals, pred_order):
        if pred not in bdf.index:
            continue
        r = bdf.loc[pred]
        c = PRED_COLORS[pred]
        ax.errorbar(r['adj_beta'], yp,
                    xerr=[[r['adj_beta'] - r['adj_ci_lo']],
                          [r['adj_ci_hi'] - r['adj_beta']]],
                    fmt='o', color=c, ecolor=c, elinewidth=1.6, capsize=4,
                    markersize=8, markeredgecolor='black', markeredgewidth=0.5)
        stars = ('***' if r['perm_p'] < 1e-3 else '**' if r['perm_p'] < 1e-2
                 else '*' if r['perm_p'] < 0.05 else 'ns')
        ax.text(r['adj_ci_hi'], yp, f"  {stars}", fontsize=10, va='center')
    ax.set_yticks(yvals)
    ax.set_yticklabels([PRED_DISPLAY[p] for p in pred_order], fontsize=10)
    ax.set_xlabel('Adjusted β  (Δ within-network percentile per SD; 95% CI)', fontsize=10.5)
    ax.set_title(f'{code} — {blabel}', fontsize=11.5, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=9)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

handles = [mlines.Line2D([], [], color=PRED_COLORS[p], marker='o', linestyle='',
                         markersize=8, markeredgecolor='black', markeredgewidth=0.5,
                         label=PRED_DISPLAY[p]) for p in pred_order]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=9.5,
           frameon=True, bbox_to_anchor=(0.5, -0.01))
fig.suptitle('Within-network predictors of branch-specific PBS — combined dataset '
             '(HGDP+1KGP+SGDP)\n'
             'outcome = within-network percentile; * = permutation p < 0.05  '
             '(SGDP-only genome-wide version is null — see caveats)',
             fontsize=13, fontweight='bold', y=0.995)
fig.tight_layout(rect=[0, 0.03, 1, 0.99])

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_within_network_combined.{ext}')
    fig.savefig(path, dpi=220, bbox_inches='tight')
    print(f"  Saved → {path}")
plt.close(fig)
