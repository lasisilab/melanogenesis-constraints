"""
phase1_gtex_analysis.py

Phase 1.2: GTEx tissue expression breadth for the Raghunath melanogenesis network.

Downloads GTEx v8 gene-level median TPM (if not already present), computes
per-gene tissue breadth (number of tissues with median TPM > 1), and produces:

  1. data/gtex_tissue_breadth.csv       — per-gene tissue breadth scores
  2. data/network_constraint_gtex.csv   — merged dataset (LOEUF + tissue breadth)
  3. output/figure_phase1_gtex.png/pdf  — two-panel figure:
       Panel A: Tissue breadth vs. LOEUF scatter
       Panel B: Tissue breadth by functional category (boxplot)

Hypothesis 2 test (regression): LOEUF ~ tissue_breadth + functional_category

GTEx v8 source:
  https://gtexportal.org/home/datasets
  File: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz (~80 MB)
"""

import os
import sys
import urllib.request
import gzip
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf
import scikit_posthocs as sp
from adjustText import adjust_text

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOEUF_FILE   = os.path.join(PROJECT_DIR, 'LOEUF_by_functional_category.xlsx')
GTEX_GCT     = os.path.join(PROJECT_DIR, 'data', 'GTEx_v8_gene_median_tpm.gct.gz')
BREADTH_CSV  = os.path.join(PROJECT_DIR, 'data', 'gtex_tissue_breadth.csv')
OUT_CSV      = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
OUT_DIR      = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
TPM_THRESHOLD = 1.0   # tissue breadth: TPM > this value

# ── Step 1: Download GTEx data ─────────────────────────────────────────────
if not os.path.exists(GTEX_GCT):
    print(f"Downloading GTEx v8 median TPM (~80 MB)...")
    print(f"  URL: {GTEX_URL}")
    try:
        urllib.request.urlretrieve(GTEX_URL, GTEX_GCT)
        print(f"  Saved → {GTEX_GCT}")
    except Exception as e:
        print(f"\nERROR: Could not download GTEx file: {e}")
        print(
            "\nPlease download the file manually from the GTEx portal:\n"
            "  https://gtexportal.org/home/datasets\n"
            "  File: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz\n"
            f"  Save to: {GTEX_GCT}\n"
            "Then re-run this script."
        )
        sys.exit(1)
else:
    print(f"GTEx file already present: {GTEX_GCT}")

# ── Step 2: Parse GCT file and compute tissue breadth ─────────────────────
print("\nParsing GTEx GCT file...")
with gzip.open(GTEX_GCT, 'rt') as f:
    f.readline()          # #1.2 header
    dims_line = f.readline().strip().split('\t')
    n_genes, n_samples = int(dims_line[0]), int(dims_line[1])
    print(f"  {n_genes} genes × {n_samples} tissues")
    gtex_df = pd.read_csv(f, sep='\t')

# GCT columns: Name (Ensembl ID), Description (gene symbol), tissue1, tissue2, ...
tissue_cols = [c for c in gtex_df.columns if c not in ('Name', 'Description')]
print(f"  {len(tissue_cols)} tissues")

# Gene symbol = Description column; Ensembl ID stripped of version
gtex_df['gene'] = gtex_df['Description'].str.upper()
gtex_df['ensembl_id'] = gtex_df['Name'].str.split('.').str[0]

# Compute tissue breadth
tpm_matrix = gtex_df[tissue_cols].values.astype(float)
gtex_df['tissue_breadth']      = (tpm_matrix > TPM_THRESHOLD).sum(axis=1)
gtex_df['max_tpm']             = tpm_matrix.max(axis=1)
gtex_df['median_tpm_all']      = np.median(tpm_matrix, axis=1)

breadth_df = gtex_df[['gene', 'ensembl_id', 'tissue_breadth',
                        'max_tpm', 'median_tpm_all']].copy()

# Keep one row per gene symbol (take max tissue_breadth if duplicates)
breadth_df = (breadth_df.sort_values('tissue_breadth', ascending=False)
              .drop_duplicates(subset='gene', keep='first')
              .reset_index(drop=True))

breadth_df.to_csv(BREADTH_CSV, index=False)
print(f"Tissue breadth saved → {BREADTH_CSV}")
print(f"  Median breadth (all {len(breadth_df)} genes): "
      f"{breadth_df['tissue_breadth'].median():.0f} tissues")

# ── Step 3: Merge with LOEUF network data ─────────────────────────────────
print("\nMerging with LOEUF network data...")
loeuf = pd.read_excel(LOEUF_FILE, sheet_name='All Genes by Category')
loeuf.columns = ['gene', 'functional_category', 'disease_class',
                  'LOEUF', 'pLI', 'betweenness_centrality']
loeuf['gene'] = loeuf['gene'].str.upper()

df = loeuf.merge(breadth_df[['gene', 'tissue_breadth', 'max_tpm', 'median_tpm_all']],
                 on='gene', how='left')
n_missing = df['tissue_breadth'].isna().sum()
print(f"  {len(df)} network genes, {n_missing} without GTEx data")
if n_missing > 0:
    print("  Missing:", df.loc[df['tissue_breadth'].isna(), 'gene'].tolist())

df.to_csv(OUT_CSV, index=False)
print(f"Merged data saved → {OUT_CSV}")

# ── Step 4: Statistics ─────────────────────────────────────────────────────
df_clean = df.dropna(subset=['tissue_breadth', 'LOEUF'])
rho, pval = stats.spearmanr(df_clean['tissue_breadth'], df_clean['LOEUF'])
print(f"\nSpearman ρ(tissue_breadth, LOEUF) = {rho:.3f}, p = {pval:.3e}")

# Kruskal-Wallis across categories
CATEGORY_ORDER = [
    'Pigment-specific',
    'Developmental/NC',
    'Generic signaling',
    'Cytokines/growth factors',
    'Apoptosis/cell death',
    'Other',
]
groups = [df_clean.loc[df_clean['functional_category'] == c, 'tissue_breadth'].values
          for c in CATEGORY_ORDER]
kw_stat, kw_p = stats.kruskal(*[g for g in groups if len(g) > 0])
print(f"Kruskal-Wallis tissue_breadth across categories: H = {kw_stat:.2f}, p = {kw_p:.3e}")

dunn = sp.posthoc_dunn(
    df_clean, val_col='tissue_breadth', group_col='functional_category',
    p_adjust='bonferroni')
print("\nDunn's posthoc (Bonferroni):")
cats_present = [c for c in CATEGORY_ORDER if c in df_clean['functional_category'].values]
for i in range(len(cats_present)):
    for j in range(i + 1, len(cats_present)):
        p = dunn.loc[cats_present[i], cats_present[j]]
        if p < 0.05:
            print(f"  {cats_present[i]} vs {cats_present[j]}: p = {p:.4f}")

# Hypothesis 2: Regression LOEUF ~ tissue_breadth + functional_category
print("\nHypothesis 2 — Multiple regression: LOEUF ~ tissue_breadth + functional_category")
reg_df = df_clean.copy()
reg_df['functional_category_r'] = pd.Categorical(
    reg_df['functional_category'], categories=CATEGORY_ORDER)
model = smf.ols(
    'LOEUF ~ tissue_breadth + C(functional_category, Treatment("Pigment-specific"))',
    data=reg_df).fit()
print(model.summary().tables[1])

# ── Step 5: Figure ─────────────────────────────────────────────────────────
CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',
    'Developmental/NC':          '#E8907E',
    'Generic signaling':         '#F5C242',
    'Cytokines/growth factors':  '#4878CF',
    'Apoptosis/cell death':      '#6BAD6B',
    'Other':                     '#B0B0B0',
}
CATEGORY_SHORT = {
    'Pigment-specific':          'Pigment-\nspecific',
    'Developmental/NC':          'Developmental\n/NC',
    'Generic signaling':         'Generic\nsignaling',
    'Cytokines/growth factors':  'Cytokines/\nGF',
    'Apoptosis/cell death':      'Apoptosis/\ncell death',
    'Other':                     'Other',
}
LABEL_GENES = {'TYR', 'TYRP1', 'DCT', 'OCA2', 'MC1R', 'SOX10', 'MITF',
               'PAX3', 'TFAP2A', 'AKT1', 'TP53', 'MAPK1', 'NFKB1', 'STAT3'}

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.subplots_adjust(wspace=0.32, left=0.07, right=0.98, top=0.88, bottom=0.15)

# ── Panel A: Tissue breadth vs. LOEUF scatter ─────────────────────────────
ax = axes[0]
for cat in CATEGORY_ORDER:
    sub = df_clean[df_clean['functional_category'] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(sub['tissue_breadth'], sub['LOEUF'],
               c=CATEGORY_COLORS[cat], s=65, alpha=0.85,
               edgecolors='white', linewidths=0.4,
               label=f'{cat} (n={len(sub)})', zorder=3)

slope, intercept, *_ = stats.linregress(df_clean['tissue_breadth'], df_clean['LOEUF'])
x_line = np.linspace(df_clean['tissue_breadth'].min(),
                     df_clean['tissue_breadth'].max(), 100)
ax.plot(x_line, slope * x_line + intercept, '--',
        color='#555555', alpha=0.6, lw=1, zorder=1)

texts = []
for _, row in df_clean.iterrows():
    if row['gene'] in LABEL_GENES:
        texts.append(ax.text(row['tissue_breadth'], row['LOEUF'],
                             row['gene'], fontsize=13, fontweight='bold',
                             color='#333333', style='italic'))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7),
            force_points=(2.5, 2.5), force_text=(2.5, 2.5),
            expand_points=(3.0, 3.0), expand_text=(2.0, 2.0),
            iter_lim=1000)

ax.text(0.97, 0.97, f'Spearman ρ = {rho:.3f}\np = {pval:.2e}',
        transform=ax.transAxes, fontsize=14, ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.9))

ax.set_xlabel(f'Tissue breadth (# tissues with TPM > {TPM_THRESHOLD:.0f})', fontsize=16)
ax.set_ylabel('LOEUF (lower = more constrained)', fontsize=16)
ax.set_title('Tissue expression breadth vs.\nLoF intolerance', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=12, loc='upper left', framealpha=0.9,
          edgecolor='gray', handletextpad=0.4, borderpad=0.5, markerscale=1.3)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

# ── Panel B: Tissue breadth boxplot by category ───────────────────────────
ax = axes[1]
bp_data, colors, ns = [], [], []
for cat in CATEGORY_ORDER:
    vals = df_clean.loc[df_clean['functional_category'] == cat, 'tissue_breadth'].values
    bp_data.append(vals)
    colors.append(CATEGORY_COLORS[cat])
    ns.append(len(vals))

positions = list(range(len(CATEGORY_ORDER)))
bp = ax.boxplot(bp_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=False, zorder=2,
                medianprops=dict(color='black', linewidth=1.5),
                whiskerprops=dict(color='gray'), capprops=dict(color='gray'))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.85)
    patch.set_edgecolor('gray')

rng = np.random.default_rng(42)
for i, (vals, cat) in enumerate(zip(bp_data, CATEGORY_ORDER)):
    if len(vals) == 0:
        continue
    jitter = rng.uniform(-0.15, 0.15, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals,
               c=CATEGORY_COLORS[cat], s=25, alpha=0.8,
               edgecolors='white', linewidths=0.3, zorder=3)

for i, n in enumerate(ns):
    ax.text(i, -2, f'n={n}', ha='center', fontsize=11, color='#555555',
            fontstyle='italic', zorder=4)


def draw_bracket(ax, x1, x2, y, p_text, lw=1.0, fontsize=14):
    h = 1.0
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h], color='black', lw=lw, zorder=5)
    ax.text((x1 + x2) / 2, y + 0.3, p_text,
            ha='center', va='bottom', fontsize=fontsize, zorder=5)


sig_pairs = []
for i in range(len(CATEGORY_ORDER)):
    for j in range(i + 1, len(CATEGORY_ORDER)):
        c1, c2 = CATEGORY_ORDER[i], CATEGORY_ORDER[j]
        if c1 in dunn.index and c2 in dunn.columns:
            p = dunn.loc[c1, c2]
            if p < 0.05:
                sig_pairs.append((i, j, p))
sig_pairs.sort(key=lambda x: x[1] - x[0])

y_vals_all = [v for vals in bp_data for v in vals]
y_max = max(y_vals_all) if y_vals_all else 60
bracket_y_start = y_max + 3
bracket_spacing = 6
for k, (i, j, p) in enumerate(sig_pairs):
    y_bracket = bracket_y_start + k * bracket_spacing
    p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
    draw_bracket(ax, i, j, y_bracket, p_text)

top_y = bracket_y_start + max(len(sig_pairs), 1) * bracket_spacing + 3
ax.text(len(CATEGORY_ORDER) / 2 - 0.5, top_y,
        f'Kruskal-Wallis p = {kw_p:.2e}',
        fontsize=13, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.9))

ax.set_xticks(positions)
ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                   fontsize=13, rotation=30, ha='right')
ax.set_ylabel(f'Tissue breadth (# tissues with TPM > {TPM_THRESHOLD:.0f})', fontsize=16)
ax.set_title('Expression breadth by\nfunctional category', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

for ext in ('png', 'pdf'):
    out_path = os.path.join(OUT_DIR, f'figure_phase1_gtex.{ext}')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {out_path}")
plt.close(fig)

# ── Summary table ──────────────────────────────────────────────────────────
print("\nTissue breadth summary by functional category:")
summary = (df_clean.groupby('functional_category')['tissue_breadth']
           .agg(['count', 'median', 'mean'])
           .round(1)
           .rename(columns={'count': 'N', 'median': 'Median breadth',
                            'mean': 'Mean breadth'}))
summary = summary.loc[[c for c in CATEGORY_ORDER if c in summary.index]]
print(summary.to_string())

print("\nDone!")
