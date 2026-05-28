"""
phase2_pbs_within_network.py

Within-network PBS association analysis (per the PEQG steering note).

Among the 129 melanogenesis-network genes, tests whether gene-level features
predict branch-specific, genome-wide-percentile PBS for each focal branch
(AFR, MEL, EUR, SAS). Genome-wide data is used only to normalize raw PBS into
branch-specific percentiles — this is a WITHIN-network association analysis,
not a network-vs-genome test.

Pipeline (see analysis/PBS_WITHIN_NETWORK_PLAN.md):
  1. Concatenate genome-wide SGDP PBS (pbs_sgdp_genomewide_chr{1..22}.csv)
  2. Gene-level PBS per branch:
       AFR = pbs1_african, MEL = pbs3_melanesian (precomputed)
       EUR = (T(afr,eur)+T(eur,mel)-T(afr,mel))/2
       SAS = (T(afr,sas)+T(sas,mel)-T(afr,mel))/2   (Afr-Mel reference)
  3. Branch-specific genome-wide percentile (rank vs all genes, n_snps>0)
  4. Subset to 129 network genes
  5. Predictors: tau, betweenness, log1p(KEGG), LOEUF, tissue_breadth (z-scored)
     Covariates: log1p(n_snps), log10(body_length)  (z-scored)
  6. Per branch x predictor: Spearman rho; adjusted single-predictor OLS
     (beta, 95% CI, p, VIF); within-network permutation (10k) -> empirical p
  7. Outputs: phase2_pbs_within_network.{txt,csv} + figure (2x2 forest)

Outputs (output/):
  phase2_pbs_within_network.txt        — full table + per-model detail
  phase2_pbs_within_network.csv        — machine-readable summary
  figure_phase2_pbs_within_network.png — 2x2 forest (one panel per branch)
  figure_phase2_pbs_within_network.pdf
"""

import os
import glob
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
GW_GLOB     = os.path.join(OUT_DIR, 'pbs_sgdp_genomewide_chr*.csv')

N_PERM = 10000
RNG = np.random.default_rng(0)
PAD = 20000   # ±10 kb flanking added by the gene-window BED


# ── Step 1: assemble genome-wide scan ──────────────────────────────────────
gw_files = sorted(glob.glob(GW_GLOB))
if not gw_files:
    raise FileNotFoundError(f"No genome-wide CSVs matching {GW_GLOB}")
gw = pd.concat([pd.read_csv(f) for f in gw_files], ignore_index=True)
gw['gene_name'] = gw['gene_name'].str.upper()
print(f"Step 1: {len(gw):,} genome-wide gene rows from {len(gw_files)} files")


# ── Step 2: gene-level PBS for 4 branches ──────────────────────────────────
def _t(fst):
    fst = np.clip(pd.to_numeric(fst, errors='coerce').fillna(0).clip(lower=0), 0, 0.9999)
    return -np.log(1.0 - fst)

t_afr_eur = _t(gw['fst_afr_eur'])
t_afr_mel = _t(gw['fst_afr_mel'])
t_afr_sas = _t(gw['fst_afr_sas'])
t_eur_mel = _t(gw['fst_eur_mel'])
t_mel_sas = _t(gw['fst_mel_sas'])

gw['pbs_afr'] = gw['pbs1_african'].clip(lower=0)
gw['pbs_mel'] = gw['pbs3_melanesian'].clip(lower=0)
gw['pbs_eur'] = ((t_afr_eur + t_eur_mel - t_afr_mel) / 2).clip(lower=0)
gw['pbs_sas'] = ((t_afr_sas + t_mel_sas - t_afr_mel) / 2).clip(lower=0)

BRANCHES = [
    ('pbs_afr', 'AFR', 'African (PBS-1, S.Asian outgroup)'),
    ('pbs_mel', 'MEL', 'Melanesian (PBS-3, S.Asian outgroup)'),
    ('pbs_eur', 'EUR', 'European (Afr–Mel reference)'),
    ('pbs_sas', 'SAS', 'South Asian (Afr–Mel reference)'),
]
print("Step 2: gene-level PBS computed for AFR, MEL, EUR, SAS")


# ── Step 3: branch-specific genome-wide percentiles ────────────────────────
for col, code, _ in BRANCHES:
    pct = np.full(len(gw), np.nan)
    valid = gw[col].notna().values & (gw['n_snps'] > 0).values
    vals = gw.loc[valid, col].values
    # percentile rank 0–100 among valid genome-wide genes
    ranks = stats.rankdata(vals, method='average')
    pct[valid] = 100.0 * (ranks - 1) / (len(vals) - 1)
    gw[f'{code}_pctile'] = pct
print("Step 3: branch-specific genome-wide percentiles assigned")


# ── Step 4: subset to network genes + merge features ───────────────────────
master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
kegg = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()
feat = master[['gene', 'functional_category', 'tau', 'betweenness_centrality',
               'LOEUF', 'tissue_breadth']].merge(kegg, on='gene', how='left')

net = gw.merge(feat, left_on='gene_name', right_on='gene', how='inner')
print(f"Step 4: {len(net)} network genes matched in the genome-wide scan")


# ── Step 5: predictors, covariates, transforms ─────────────────────────────
net['kegg_log1p']    = np.log1p(net['kegg_pathway_count'])
net['nsnps_log1p']   = np.log1p(net['n_snps'])
net['body_length']   = (net['end'] - net['start'] - PAD).clip(lower=1)
net['glen_log10']    = np.log10(net['body_length'])

PREDICTORS = [
    ('tau',                    'τ (tissue specificity)'),
    ('betweenness_centrality', 'Betweenness centrality'),
    ('kegg_log1p',             'KEGG  log1p(pathway count)'),
    ('LOEUF',                  'LOEUF (constraint)'),
    ('tissue_breadth',         'Tissue breadth'),
]
COVARS = ['nsnps_log1p', 'glen_log10']

PRED_COLORS = {
    'tau':                    '#d62728',   # red
    'betweenness_centrality': '#9467bd',   # purple
    'kegg_log1p':             '#2ca02c',   # green
    'LOEUF':                  '#7f7f7f',   # gray
    'tissue_breadth':         '#1f77b4',   # blue
}
PRED_DISPLAY = dict(PREDICTORS)


def _z(s):
    s = pd.to_numeric(s, errors='coerce')
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)


def _dir(beta, pred):
    """Human-readable direction for the table."""
    if pred == 'LOEUF':
        return ('less-constrained genes higher' if beta > 0
                else 'more-constrained genes higher')
    if pred == 'tau':
        return ('tissue-specific genes higher' if beta > 0
                else 'broadly-expressed genes higher')
    if pred == 'betweenness_centrality':
        return ('central genes higher' if beta > 0 else 'peripheral genes higher')
    if pred == 'kegg_log1p':
        return ('more-pathway genes higher' if beta > 0 else 'fewer-pathway genes higher')
    if pred == 'tissue_breadth':
        return ('broadly-expressed genes higher' if beta > 0
                else 'tissue-restricted genes higher')
    return ''


# ── Step 6: tests per branch × predictor ───────────────────────────────────
rows = []
detail = ["Within-network PBS association — branch-specific genome-wide percentile",
          "=" * 78,
          "Outcome: branch-specific genome-wide PBS percentile (0–100), among the "
          "network genes.",
          "Adjusted models: pctile ~ z(predictor) + z(log1p n_snps) + z(log10 body_len).",
          f"Permutation: {N_PERM:,} shuffles of the percentile among network genes.",
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

        # (a) Spearman (unadjusted)
        rho, p_sp = stats.spearmanr(sub[pred].values, y)

        # (b) adjusted single-predictor OLS
        Xcols = [pred] + COVARS
        X = pd.DataFrame({c: _z(sub[c]) for c in Xcols}, index=sub.index)
        Xc = sm.add_constant(X.values)
        res = sm.OLS(y, Xc).fit()
        beta = res.params[1]
        se   = res.bse[1]
        p_b  = res.pvalues[1]
        ci   = res.conf_int(0.05)[1]
        # VIF on predictor+covariates
        Xv = np.column_stack([np.ones(len(X)), X.values])
        vif = variance_inflation_factor(Xv, 1)

        # (c) within-network permutation of the adjusted beta
        Xc_perm = sm.add_constant(X.values)
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            yp = RNG.permutation(y)
            null[i] = sm.OLS(yp, Xc_perm).fit().params[1]
        p_perm = (np.sum(np.abs(null) >= abs(beta)) + 1) / (N_PERM + 1)

        rows.append({
            'branch': code, 'predictor': pred, 'n': n,
            'spearman_rho': rho, 'spearman_p': p_sp,
            'adj_beta': beta, 'adj_se': se, 'adj_ci_lo': ci[0], 'adj_ci_hi': ci[1],
            'adj_p': p_b, 'perm_p': p_perm, 'vif': vif,
            'direction': _dir(beta, pred),
        })
        stars = ('***' if p_perm < 1e-3 else '**' if p_perm < 1e-2
                 else '*' if p_perm < 0.05 else 'ns')
        detail.append(
            f"  {pdisp:<28} ρ={rho:+.3f} (p={p_sp:.2e}) | "
            f"adjβ={beta:+.2f} [{ci[0]:+.2f},{ci[1]:+.2f}] p={p_b:.2e} | "
            f"perm p={p_perm:.4f} {stars} | VIF={vif:.2f} | {_dir(beta, pred)}")

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT_DIR, 'phase2_pbs_within_network.csv'), index=False)

# Suggested table (compact) appended to the text report
detail.append(f"\n\n{'='*78}\nSUMMARY TABLE  (β on 0–100 percentile scale)\n{'='*78}")
detail.append(f"{'branch':<7}{'predictor':<26}{'rho':>7}{'adjβ':>8}{'adjp':>11}{'permp':>9}  direction")
for _, r in summary.iterrows():
    detail.append(f"{r['branch']:<7}{PRED_DISPLAY[r['predictor']][:25]:<26}"
                  f"{r['spearman_rho']:>+7.3f}{r['adj_beta']:>+8.2f}"
                  f"{r['adj_p']:>11.2e}{r['perm_p']:>9.4f}  {r['direction']}")

detail += [
    "", "-" * 78,
    "CAVEATS:",
    "  • EUR/SAS PBS use the Afr–Mel reference pair (no eur–sas FST); valid",
    "    lineage-specific drift but not perfectly parallel to AFR/MEL scans.",
    "  • Melanesian PBS rests on 17 SGDP samples — underpowered, drift/admixture.",
    "  • PBS is a frequency-differentiation statistic, NOT proof of selection.",
    "  • EUR is a contrast / positive control (high EUR PBS at pigment genes =",
    "    European depigmentation, not convergent darkening).",
    "  • Covariates omitted (not on disk): recombination rate, GC, mappability,",
    "    background-selection score.",
]

with open(os.path.join(OUT_DIR, 'phase2_pbs_within_network.txt'), 'w') as f:
    f.write("\n".join(detail) + "\n")
print("Step 7: wrote phase2_pbs_within_network.{txt,csv}")


# ── Step 7: figure — 2×2 forest, one panel per branch ──────────────────────
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
    ax.set_xlabel('Adjusted β  (Δ percentile per SD; 95% CI)', fontsize=10.5)
    ax.set_title(f'{code} — {blabel}', fontsize=11.5, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=9)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

handles = [mlines.Line2D([], [], color=PRED_COLORS[p], marker='o', linestyle='',
                         markersize=8, markeredgecolor='black', markeredgewidth=0.5,
                         label=PRED_DISPLAY[p]) for p in pred_order]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=9.5,
           frameon=True, bbox_to_anchor=(0.5, -0.01))
fig.suptitle('Within-network predictors of branch-specific genome-wide-percentile PBS\n'
             '(* = within-network permutation p < 0.05)',
             fontsize=14, fontweight='bold', y=0.995)
fig.tight_layout(rect=[0, 0.03, 1, 0.99])

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_within_network.{ext}')
    fig.savefig(path, dpi=220, bbox_inches='tight')
    print(f"  Saved → {path}")
plt.close(fig)
