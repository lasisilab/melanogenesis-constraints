"""
phase2_peripheral_pigmentation.py

Hypothesis-driven integration of the network-level features for the
melanogenesis network:

    Are peripheral, narrow-function genes
        (low betweenness centrality + few KEGG pathways)
    also tissue-specific to skin
        (high tau + skin-enriched GTEx expression)
    and do those genes carry distinct constraint / selection signals
        (LOEUF, PBS in MEL/AFR, nucleotide diversity pi)?

Outputs (in output/):
  - peripheral_pigmentation_pergene.csv     unified per-gene matrix + cohort
  - peripheral_pigmentation_cohort_stats.csv  Mann-Whitney peripheral vs core
  - peripheral_pigmentation_corr.csv        Spearman matrix + BH FDR
  - figure_peripheral_pigmentation.{png,pdf}

Cohorts are terciles of (betweenness_centrality, kegg_pathway_count):
  peripheral = bottom tercile of BOTH
  core       = top tercile of BOTH
  mixed      = remainder
"""

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import false_discovery_control

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(PROJECT_DIR, 'data')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

SKIN_TISSUES = [
    'Skin - Sun Exposed (Lower leg)',
    'Skin - Not Sun Exposed (Suprapubic)',
]
FIBROBLAST_TISSUE = 'Cells - Cultured fibroblasts'

# ── Load per-gene tables ───────────────────────────────────────────────────
cat  = pd.read_csv(os.path.join(DATA_DIR, 'network_constraint_categorized.csv'))
kegg = pd.read_csv(os.path.join(DATA_DIR, 'kegg_pathway_counts.csv'))
pbs  = pd.read_csv(os.path.join(DATA_DIR, 'pbs_per_gene.csv'))
pi   = pd.read_csv(os.path.join(DATA_DIR, 'pi_per_gene.csv'))

for d in (cat, kegg, pbs, pi):
    d['gene'] = d['gene'].str.upper()

# ── Build skin-specificity from raw GTEx ───────────────────────────────────
print("Parsing GTEx GCT for per-tissue TPM...")
with gzip.open(os.path.join(DATA_DIR, 'GTEx_v8_gene_median_tpm.gct.gz'), 'rt') as f:
    f.readline(); f.readline()
    gtex = pd.read_csv(f, sep='\t')

TISSUES = [c for c in gtex.columns if c not in ('Name', 'Description')]
gtex['gene'] = gtex['Description'].str.upper()
tpm = gtex.groupby('gene')[TISSUES].max().reset_index()

# Subset to network genes
tpm = tpm[tpm['gene'].isin(cat['gene'])].reset_index(drop=True)
log = np.log2(tpm[TISSUES].values + 1)

skin_idx     = [TISSUES.index(t) for t in SKIN_TISSUES]
fibro_idx    = TISSUES.index(FIBROBLAST_TISSUE)
nonskin_idx  = [i for i in range(len(TISSUES))
                if i not in skin_idx and i != fibro_idx]

skin_max         = log[:, skin_idx].max(axis=1)
nonskin_mean     = log[:, nonskin_idx].mean(axis=1)
nonskin_sd       = log[:, nonskin_idx].std(axis=1)
nonskin_sd_safe  = np.where(nonskin_sd == 0, np.nan, nonskin_sd)
skin_z           = (skin_max - nonskin_mean) / nonskin_sd_safe

# Rank of best skin tissue among all 54 tissues (1 = highest expression)
ranks = (-log).argsort(axis=1).argsort(axis=1) + 1   # rank of each col per row
skin_rank_min    = ranks[:, skin_idx].min(axis=1)

# Is the top tissue a skin tissue?
top_tissue_idx   = log.argmax(axis=1)
is_skin_top      = np.isin(top_tissue_idx, skin_idx)

skin_df = pd.DataFrame({
    'gene':          tpm['gene'],
    'skin_log':      skin_max,
    'nonskin_mean':  nonskin_mean,
    'skin_z':        skin_z,
    'skin_rank':     skin_rank_min,
    'is_skin_top':   is_skin_top,
})

# ── Merge everything ───────────────────────────────────────────────────────
df = (cat[['gene', 'functional_category', 'LOEUF',
           'betweenness_centrality', 'tau', 'tissue_breadth', 'max_tpm']]
      .merge(kegg[['gene', 'kegg_pathway_count']], on='gene', how='left')
      .merge(skin_df, on='gene', how='left')
      .merge(pbs[['gene', 'pbs1_african', 'pbs3_melanesian']],
             on='gene', how='left')
      .merge(pi[['gene', 'pi_african', 'pi_melanesian']],
             on='gene', how='left'))

df['kegg_pathway_count'] = df['kegg_pathway_count'].fillna(0)

print(f"  merged matrix: {len(df)} genes × {df.shape[1]} features")

# ── Define cohorts (terciles of centrality & KEGG count) ───────────────────
c_low, c_high = df['betweenness_centrality'].quantile([1/3, 2/3])
k_low, k_high = df['kegg_pathway_count'].quantile([1/3, 2/3])

def assign(row):
    c, k = row['betweenness_centrality'], row['kegg_pathway_count']
    if c <= c_low and k <= k_low:
        return 'peripheral'
    if c >= c_high and k >= k_high:
        return 'core'
    return 'mixed'

df['cohort'] = df.apply(assign, axis=1)
n_per  = (df['cohort'] == 'peripheral').sum()
n_core = (df['cohort'] == 'core').sum()
n_mix  = (df['cohort'] == 'mixed').sum()
print(f"  cohorts → peripheral={n_per}, core={n_core}, mixed={n_mix}")

# ── Cohort comparisons (peripheral vs core) ────────────────────────────────
FEATURES = [
    ('tau',                    'tau (tissue specificity)'),
    ('skin_z',                 'skin z-score'),
    ('skin_rank',              'skin rank (1=highest)'),
    ('LOEUF',                  'LOEUF'),
    ('pbs1_african',           'PBS African'),
    ('pbs3_melanesian',        'PBS Melanesian'),
    ('pi_african',             'pi African'),
    ('pi_melanesian',          'pi Melanesian'),
]

stat_rows = []
per_vals  = df[df['cohort'] == 'peripheral']
core_vals = df[df['cohort'] == 'core']
for col, label in FEATURES:
    a = per_vals[col].dropna().values
    b = core_vals[col].dropna().values
    if len(a) < 3 or len(b) < 3:
        u, p, r = np.nan, np.nan, np.nan
    else:
        u, p = stats.mannwhitneyu(a, b, alternative='two-sided')
        # rank-biserial effect size
        r = 1 - 2 * u / (len(a) * len(b))
    stat_rows.append({
        'feature':           col,
        'label':             label,
        'n_peripheral':      len(a),
        'n_core':            len(b),
        'median_peripheral': np.median(a) if len(a) else np.nan,
        'median_core':       np.median(b) if len(b) else np.nan,
        'mw_U':              u,
        'mw_p':              p,
        'rank_biserial_r':   r,
    })
stat_df = pd.DataFrame(stat_rows)
stat_df['mw_p_bh'] = false_discovery_control(stat_df['mw_p'].fillna(1))
stat_df.to_csv(os.path.join(OUT_DIR,
               'peripheral_pigmentation_cohort_stats.csv'), index=False)

print("\nPeripheral vs Core (Mann-Whitney, BH-corrected):")
for _, r in stat_df.iterrows():
    sig = '***' if r['mw_p_bh'] < 0.001 else '**' if r['mw_p_bh'] < 0.01 \
          else '*' if r['mw_p_bh'] < 0.05 else ''
    print(f"  {r['label']:32s} peri med={r['median_peripheral']:.3g}  "
          f"core med={r['median_core']:.3g}  p={r['mw_p']:.2e}  "
          f"q={r['mw_p_bh']:.2e}  r={r['rank_biserial_r']:+.2f} {sig}")

# ── Spearman correlation matrix + BH FDR ───────────────────────────────────
CORR_FEATURES = ['betweenness_centrality', 'kegg_pathway_count', 'LOEUF',
                 'tau', 'skin_z', 'skin_rank',
                 'pbs1_african', 'pbs3_melanesian',
                 'pi_african', 'pi_melanesian']
sub = df[CORR_FEATURES]
rho  = pd.DataFrame(np.nan, index=CORR_FEATURES, columns=CORR_FEATURES)
pval = pd.DataFrame(np.nan, index=CORR_FEATURES, columns=CORR_FEATURES)
for i, a in enumerate(CORR_FEATURES):
    for j, b in enumerate(CORR_FEATURES):
        if i == j:
            rho.loc[a, b] = 1.0
            pval.loc[a, b] = 0.0
            continue
        ok = sub[[a, b]].dropna()
        if len(ok) < 5:
            continue
        r_, p_ = stats.spearmanr(ok.iloc[:, 0], ok.iloc[:, 1])
        rho.loc[a, b]  = float(r_)
        pval.loc[a, b] = float(p_)

# BH across off-diagonal upper-triangle p-values
flat_ps, idx = [], []
for i in range(len(CORR_FEATURES)):
    for j in range(i + 1, len(CORR_FEATURES)):
        a, b = CORR_FEATURES[i], CORR_FEATURES[j]
        if not np.isnan(pval.loc[a, b]):
            flat_ps.append(pval.loc[a, b]); idx.append((a, b))
qs = false_discovery_control(np.asarray(flat_ps))
qmat = pd.DataFrame(np.nan, index=CORR_FEATURES, columns=CORR_FEATURES)
for (a, b), q in zip(idx, qs):
    qmat.loc[a, b] = q
    qmat.loc[b, a] = q

corr_out = rho.copy()
corr_out.to_csv(os.path.join(OUT_DIR, 'peripheral_pigmentation_corr.csv'))
qmat.to_csv(os.path.join(OUT_DIR, 'peripheral_pigmentation_corr_qvals.csv'))

# ── Partial Spearman: feature vs skin_z, controlling LOEUF ─────────────────
def partial_spearman(x, y, z):
    ok = pd.concat([x, y, z], axis=1).dropna()
    if len(ok) < 5:
        return np.nan, np.nan
    rxy, _ = stats.spearmanr(ok.iloc[:, 0], ok.iloc[:, 1])
    rxz, _ = stats.spearmanr(ok.iloc[:, 0], ok.iloc[:, 2])
    ryz, _ = stats.spearmanr(ok.iloc[:, 1], ok.iloc[:, 2])
    denom = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    if denom == 0:
        return np.nan, np.nan
    pr = (rxy - rxz * ryz) / denom
    # approximate t-test
    n = len(ok)
    tstat = pr * np.sqrt((n - 3) / (1 - pr ** 2)) if abs(pr) < 1 else np.inf
    p_ = 2 * (1 - stats.t.cdf(abs(tstat), df=n - 3))
    return pr, p_

print("\nPartial Spearman (controlling LOEUF):")
for f in ['betweenness_centrality', 'kegg_pathway_count', 'tau']:
    pr, p_ = partial_spearman(df[f], df['skin_z'], df['LOEUF'])
    print(f"  {f:28s} vs skin_z | LOEUF :  rho={pr:+.3f}  p={p_:.2e}")

# ── Save merged per-gene table ─────────────────────────────────────────────
df.to_csv(os.path.join(OUT_DIR, 'peripheral_pigmentation_pergene.csv'),
          index=False)

# ── Figure: 3-panel integrative view ───────────────────────────────────────
print("\nDrawing figure...")
COHORT_COLOR = {'peripheral': '#D94040',
                'core':       '#4878CF',
                'mixed':      '#B0B0B0'}

fig = plt.figure(figsize=(16, 11))
gs  = fig.add_gridspec(2, 3, width_ratios=[1.2, 1.2, 1.0],
                       height_ratios=[1.1, 1.0],
                       hspace=0.40, wspace=0.32)

# ---- A : centrality vs KEGG, color = skin_z, peripheral box shaded --------
axA = fig.add_subplot(gs[0, :2])
sc = axA.scatter(df['betweenness_centrality'], df['kegg_pathway_count'],
                 c=df['skin_z'], cmap='RdBu_r',
                 vmin=-2, vmax=2, s=60, edgecolor='k', linewidth=0.4,
                 alpha=0.9, zorder=3)
axA.axvline(c_low,  color='gray', ls='--', lw=0.6)
axA.axhline(k_low,  color='gray', ls='--', lw=0.6)
axA.axvspan(df['betweenness_centrality'].min() - 1e-4, c_low,
            ymin=0, ymax=k_low / axA.get_ylim()[1] if axA.get_ylim()[1] else 1,
            color='#D94040', alpha=0.05, zorder=1)

LABEL_GENES = {'TYR', 'TYRP1', 'DCT', 'OCA2', 'MC1R', 'SLC24A5', 'SLC45A2',
               'MFSD12', 'PMEL', 'MITF', 'SOX10', 'PAX3', 'KIT', 'KITLG',
               'EDNRB', 'TP53', 'AKT1', 'MAPK1'}
from adjustText import adjust_text
texts = []
for _, r in df.iterrows():
    if r['gene'] in LABEL_GENES:
        texts.append(axA.text(r['betweenness_centrality'],
                              r['kegg_pathway_count'], r['gene'],
                              fontsize=9, fontweight='bold', style='italic'))
adjust_text(texts, ax=axA,
            arrowprops=dict(arrowstyle='-', color='#888888', lw=0.6))

axA.set_xlabel('Betweenness centrality')
axA.set_ylabel('KEGG pathway count')
axA.set_title('A. Network position vs. skin-specificity\n'
              '(point colour = skin z-score; lower-left tercile = peripheral)',
              fontsize=12, fontweight='bold', loc='left')
cb = fig.colorbar(sc, ax=axA, fraction=0.04, pad=0.02)
cb.set_label('skin z-score', fontsize=10)

# ---- B : Spearman heatmap ------------------------------------------------
axB = fig.add_subplot(gs[0, 2])
mat = rho.values.astype(float)
mask = np.zeros_like(mat, dtype=bool)
mask[np.triu_indices_from(mask, k=1)] = True
heat = np.where(mask, np.nan, mat)
im = axB.imshow(heat, cmap='RdBu_r', vmin=-1, vmax=1)
axB.set_xticks(range(len(CORR_FEATURES)))
axB.set_yticks(range(len(CORR_FEATURES)))
short = {'betweenness_centrality': 'centrality',
         'kegg_pathway_count':     'KEGG count',
         'LOEUF':                  'LOEUF',
         'tau':                    'tau',
         'skin_z':                 'skin z',
         'skin_rank':              'skin rank',
         'pbs1_african':           'PBS AFR',
         'pbs3_melanesian':        'PBS MEL',
         'pi_african':             'pi AFR',
         'pi_melanesian':          'pi MEL'}
labels = [short[f] for f in CORR_FEATURES]
axB.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
axB.set_yticklabels(labels, fontsize=9)
for i in range(len(CORR_FEATURES)):
    for j in range(len(CORR_FEATURES)):
        if mask[i, j] or np.isnan(mat[i, j]):
            continue
        q = qmat.iloc[i, j]
        star = '***' if q < 0.001 else '**' if q < 0.01 \
               else '*' if q < 0.05 else ''
        axB.text(j, i, f'{mat[i,j]:+.2f}\n{star}',
                 ha='center', va='center', fontsize=7,
                 color='white' if abs(mat[i, j]) > 0.5 else 'black')
fig.colorbar(im, ax=axB, fraction=0.045, pad=0.02).set_label(
    'Spearman ρ', fontsize=10)
axB.set_title('B. Feature correlations (BH-FDR)',
              fontsize=12, fontweight='bold', loc='left')

# ---- C : peripheral vs core strip+box for key features --------------------
axC = fig.add_subplot(gs[1, :])
PANEL_FEATURES = ['tau', 'skin_z', 'LOEUF',
                  'pbs1_african', 'pbs3_melanesian',
                  'pi_african', 'pi_melanesian']
labels_panel = ['tau', 'skin z', 'LOEUF',
                'PBS AFR', 'PBS MEL', 'pi AFR', 'pi MEL']
positions = np.arange(len(PANEL_FEATURES))
width = 0.32

# Z-normalise each feature for visual comparability
norm = df[PANEL_FEATURES].copy()
for c in PANEL_FEATURES:
    v = norm[c]
    norm[c] = (v - v.mean()) / v.std()

for k, (col, lab) in enumerate(zip(PANEL_FEATURES, labels_panel)):
    for offset, coh in zip([-width / 2, width / 2], ['peripheral', 'core']):
        vals = norm.loc[df['cohort'] == coh, col].dropna().values
        x = np.full(len(vals), positions[k] + offset) \
            + np.random.uniform(-0.06, 0.06, len(vals))
        axC.scatter(x, vals, s=28, alpha=0.6,
                    color=COHORT_COLOR[coh],
                    edgecolor='k', linewidth=0.3, zorder=2)
        # median line
        if len(vals):
            axC.plot([positions[k] + offset - width / 3,
                      positions[k] + offset + width / 3],
                     [np.median(vals)] * 2,
                     color='black', lw=1.6, zorder=3)
    # significance marker from stat_df
    q = stat_df.loc[stat_df['feature'] == col, 'mw_p_bh'].values[0]
    star = '***' if q < 0.001 else '**' if q < 0.01 \
           else '*' if q < 0.05 else 'ns'
    axC.text(positions[k], axC.get_ylim()[1] * 0.92 if False else 3.2,
             star, ha='center', fontsize=11, fontweight='bold')

axC.axhline(0, color='gray', lw=0.5)
axC.set_xticks(positions)
axC.set_xticklabels(labels_panel, fontsize=10)
axC.set_ylabel('z-scored value')
axC.set_title('C. Peripheral (red) vs. core (blue) cohorts — BH-FDR stars',
              fontsize=12, fontweight='bold', loc='left')

# Cohort legend
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker='o', color='w',
                  markerfacecolor=COHORT_COLOR[c], markeredgecolor='k',
                  markersize=8, label=f'{c} (n={(df["cohort"]==c).sum()})')
           for c in ['peripheral', 'core', 'mixed']]
axC.legend(handles=handles, loc='lower right', fontsize=9, frameon=True)

fig.suptitle(
    'Peripheral, narrow-function genes vs. skin-specific expression and selection',
    fontsize=14, fontweight='bold', y=0.995)

for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR,
                f'figure_peripheral_pigmentation.{ext}'),
                dpi=180, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  saved figure_peripheral_pigmentation.png/pdf")
print("\nDone.")
