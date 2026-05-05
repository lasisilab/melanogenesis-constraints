"""
phase1_gtex_extras.py

Phase 1.2 supplementary: three additional GTEx tissue-pattern visualizations
of the Raghunath melanogenesis network:

  1. Tissue specificity (Tau) vs. LOEUF scatter        → figure_phase1_gtex_tau.png
  2. Per-tissue effect plot (ΔLOEUF expr vs not-expr)  → figure_phase1_gtex_per_tissue.png
  3. Clustered heatmap of genes × tissues with LOEUF   → figure_phase1_gtex_heatmap.png

Tau is the standard tissue-specificity index (Yanai et al. 2005):
  τ = Σ (1 - x_i / max(x))  /  (n - 1),   where x_i = log2(TPM+1)
  τ = 1 → expressed in only one tissue; τ = 0 → uniform expression.
"""

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GTEX_GCT    = os.path.join(PROJECT_DIR, 'data', 'GTEx_v8_gene_median_tpm.gct.gz')
NETWORK_CSV = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

TPM_THRESHOLD = 1.0

CATEGORY_COLORS = {
    'Pigment-specific':         '#D94040',
    'Developmental/NC':         '#E8907E',
    'Generic signaling':        '#F5C242',
    'Cytokines/growth factors': '#4878CF',
    'Apoptosis/cell death':     '#6BAD6B',
    'Other':                    '#B0B0B0',
}

# ── Load GTEx ──────────────────────────────────────────────────────────────
print("Parsing GTEx GCT...")
with gzip.open(GTEX_GCT, 'rt') as f:
    f.readline(); f.readline()
    gtex_df = pd.read_csv(f, sep='\t')

TISSUES = [c for c in gtex_df.columns if c not in ('Name', 'Description')]
gtex_df['gene'] = gtex_df['Description'].str.upper()

# Collapse duplicate gene symbols (max TPM across rows for each tissue)
tpm = gtex_df.groupby('gene')[TISSUES].max().reset_index()

# ── Compute Tau on log2(TPM+1) ─────────────────────────────────────────────
log_tpm = np.log2(tpm[TISSUES].values + 1)
row_max = log_tpm.max(axis=1, keepdims=True)
# Avoid divide-by-zero for genes with all-zero expression
safe_max = np.where(row_max == 0, 1, row_max)
tau = np.sum(1 - log_tpm / safe_max, axis=1) / (len(TISSUES) - 1)
tau = np.where(row_max.flatten() == 0, np.nan, tau)
tpm['tau'] = tau

# ── Merge with network LOEUF ───────────────────────────────────────────────
network = pd.read_csv(NETWORK_CSV)
network['gene'] = network['gene'].str.upper()
df = network.merge(tpm[['gene', 'tau'] + TISSUES], on='gene', how='inner')
df = df.dropna(subset=['LOEUF', 'tau']).reset_index(drop=True)
print(f"  {len(df)} network genes with LOEUF + GTEx")

# ===========================================================================
# Figure 1: Tau (tissue specificity) vs. LOEUF scatter
# ===========================================================================
print("\nFigure 1: Tau vs. LOEUF...")
rho_t, p_t = stats.spearmanr(df['tau'], df['LOEUF'])
print(f"  Spearman ρ(Tau, LOEUF) = {rho_t:.3f}, p = {p_t:.3e}")

LABEL_GENES = {'TYR', 'TYRP1', 'DCT', 'OCA2', 'MC1R', 'SOX10', 'MITF',
               'PAX3', 'TFAP2A', 'AKT1', 'TP53', 'MAPK1', 'NFKB1', 'STAT3',
               'KIT', 'EDNRB', 'KITLG'}

fig, ax = plt.subplots(figsize=(11, 8))
for cat, color in CATEGORY_COLORS.items():
    sub = df[df['functional_category'] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(sub['tau'], sub['LOEUF'], c=color, s=70, alpha=0.85,
               edgecolors='white', linewidths=0.5,
               label=f'{cat} (n={len(sub)})', zorder=3)

slope, intercept, *_ = stats.linregress(df['tau'], df['LOEUF'])
xs = np.linspace(df['tau'].min(), df['tau'].max(), 100)
ax.plot(xs, slope * xs + intercept, '--', color='#555555', alpha=0.6, lw=1, zorder=1)

# Label key genes
from adjustText import adjust_text
texts = []
for _, r in df.iterrows():
    if r['gene'] in LABEL_GENES:
        texts.append(ax.text(r['tau'], r['LOEUF'], r['gene'],
                             fontsize=12, fontweight='bold',
                             color='#333333', style='italic'))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7))

ax.text(0.02, 0.98, f'Spearman ρ = {rho_t:.3f},  p = {p_t:.2e}\n'
                    f'n = {len(df)} genes',
        transform=ax.transAxes, fontsize=12, va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                  edgecolor='gray', alpha=0.9))

ax.set_xlabel('Tissue specificity (τ)\n← uniform · · · · · tissue-specific →', fontsize=14)
ax.set_ylabel('LOEUF (higher = less constrained)', fontsize=14)
ax.set_title('Tissue specificity (Tau) vs. LoF intolerance',
             fontsize=15, fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=10, loc='lower right', framealpha=0.9, edgecolor='gray',
          ncol=2)
ax.tick_params(labelsize=12)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'figure_phase1_gtex_tau.{ext}'),
                dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  Saved figure_phase1_gtex_tau.png/pdf")

# ===========================================================================
# Figure 2: Per-tissue effect — ΔLOEUF (median expressed - not expressed)
# ===========================================================================
print("\nFigure 2: Per-tissue effect on LOEUF...")
records = []
for t in TISSUES:
    expr_mask = df[t] > TPM_THRESHOLD
    expr_loeuf = df.loc[expr_mask, 'LOEUF']
    nonexpr_loeuf = df.loc[~expr_mask, 'LOEUF']
    if len(expr_loeuf) < 3 or len(nonexpr_loeuf) < 3:
        u, p = np.nan, np.nan
    else:
        u, p = stats.mannwhitneyu(expr_loeuf, nonexpr_loeuf, alternative='two-sided')
    records.append({
        'tissue': t,
        'n_expressed': int(expr_mask.sum()),
        'n_not_expressed': int((~expr_mask).sum()),
        'median_loeuf_expressed': expr_loeuf.median() if len(expr_loeuf) else np.nan,
        'median_loeuf_not_expressed': nonexpr_loeuf.median() if len(nonexpr_loeuf) else np.nan,
        'delta_loeuf': (expr_loeuf.median() - nonexpr_loeuf.median()
                        if len(expr_loeuf) and len(nonexpr_loeuf) else np.nan),
        'mw_p': p,
    })
per_tissue = pd.DataFrame(records).sort_values('delta_loeuf')
per_tissue.to_csv(os.path.join(OUT_DIR, 'table_phase1_gtex_per_tissue.csv'),
                  index=False)

fig, ax = plt.subplots(figsize=(11, 14))
y_pos = np.arange(len(per_tissue))
colors_bar = ['#4878CF' if d < 0 else '#D94040' for d in per_tissue['delta_loeuf']]
bars = ax.barh(y_pos, per_tissue['delta_loeuf'], color=colors_bar, alpha=0.85,
               edgecolor='gray', linewidth=0.4)

# Significance markers
for i, (_, row) in enumerate(per_tissue.iterrows()):
    if pd.notna(row['mw_p']) and row['mw_p'] < 0.05:
        marker = '***' if row['mw_p'] < 0.001 else '**' if row['mw_p'] < 0.01 else '*'
        x_off = 0.02 if row['delta_loeuf'] >= 0 else -0.02
        ha = 'left' if row['delta_loeuf'] >= 0 else 'right'
        ax.text(row['delta_loeuf'] + x_off, i, marker,
                ha=ha, va='center', fontsize=11, fontweight='bold')

# Add n=expressed under each bar (right side, italic)
for i, (_, row) in enumerate(per_tissue.iterrows()):
    ax.text(ax.get_xlim()[1] if False else 0, i - 0.35,
            f'n={row["n_expressed"]}', fontsize=7, color='#888888',
            ha='center', va='top', fontstyle='italic')

ax.axvline(0, color='black', lw=0.6)
ax.set_yticks(y_pos)
ax.set_yticklabels(per_tissue['tissue'], fontsize=9)
ax.set_xlabel('ΔLOEUF (median expressed − median not expressed)', fontsize=13)
ax.set_title('Per-tissue effect on LOEUF\n← expressed → more constrained · · · · · expressed → less constrained →',
             fontsize=13, fontweight='bold', loc='left', pad=10)
ax.text(0.99, 0.01, 'Bonferroni-uncorrected MW: * p<0.05  ** p<0.01  *** p<0.001',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
        color='#666666', style='italic')
ax.tick_params(axis='x', labelsize=11)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(OUT_DIR, f'figure_phase1_gtex_per_tissue.{ext}'),
                dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  Saved figure_phase1_gtex_per_tissue.png/pdf")

# ===========================================================================
# Figure 3: Clustered heatmap (genes × tissues) with LOEUF side annotation
# ===========================================================================
print("\nFigure 3: Clustered heatmap...")
expr_df = df.set_index('gene')[TISSUES]

# Sort genes by LOEUF (low = constrained at top)
gene_order_idx = df['LOEUF'].sort_values().index
gene_order    = df.loc[gene_order_idx, 'gene'].values
loeuf_order   = df.loc[gene_order_idx, 'LOEUF'].values

loeuf_min, loeuf_max = loeuf_order.min(), loeuf_order.max()
scale_vals = [v for v in [0.0, 1.0, 1.5, 2.0] if loeuf_min <= v <= loeuf_max]
scale_pos  = [(v - loeuf_min) / (loeuf_max - loeuf_min) * (len(gene_order) - 1)
              for v in scale_vals]

EXCLUDE_TISSUES = {
    'Cells - EBV-transformed lymphocytes',
    'Brain - Caudate (basal ganglia)',
    'Brain - Nucleus accumbens (basal ganglia)',
    'Brain - Putamen (basal ganglia)',
    'Brain - Anterior cingulate cortex (BA24)',
    'Minor Salivary Gland',
    'Small Intestine - Terminal Ileum',
    'Brain - Cerebellar Hemisphere',
    'Esophagus - Muscularis',
    'Esophagus - Mucosa',
    'Esophagus - Gastroesophageal Junction',
}


def draw_heatmap(tissue_subset, fname, n_tissues_label):
    sub_expr = expr_df[tissue_subset].loc[gene_order]
    log_data = np.log2(sub_expr.values + 1)

    # Hierarchical clustering of tissues (columns)
    t_link  = linkage(pdist(log_data.T, metric='correlation'), method='average')
    t_order = leaves_list(t_link)
    heatmap_log    = log_data[:, t_order]
    tissue_labels  = [tissue_subset[i] for i in t_order]

    n = len(gene_order)
    fig = plt.figure(figsize=(22, 22))
    left         = 0.06
    bottom       = 0.13
    height       = 0.77
    names_width  = 0.08
    strip_width  = 0.013
    scale_gap    = 0.045   # space for LOEUF scale tick labels
    heatmap_width = 0.67
    gap          = 0.005

    ax_names = fig.add_axes([left, bottom, names_width, height])
    ax_loeuf = fig.add_axes([left + names_width + gap, bottom, strip_width, height])
    ax_h     = fig.add_axes([left + names_width + gap + strip_width + scale_gap,
                             bottom, heatmap_width, height])

    # Gene names
    ax_names.set_yticks(np.arange(n))
    ax_names.set_yticklabels(gene_order, fontsize=7.5)
    ax_names.set_ylim(-0.5, n - 0.5)
    ax_names.invert_yaxis()
    ax_names.set_xticks([])
    ax_names.tick_params(axis='y', length=0, pad=3)
    for sp in ax_names.spines.values():
        sp.set_visible(False)

    # LOEUF bar
    ax_loeuf.imshow(loeuf_order.reshape(-1, 1), aspect='auto', cmap='viridis_r')
    ax_loeuf.set_xticks([])
    ax_loeuf.set_yticks([])
    ax_loeuf.set_frame_on(False)

    # LOEUF scale on right side (within scale_gap space)
    ax_loeuf_r = ax_loeuf.twinx()
    ax_loeuf_r.set_ylim(ax_loeuf.get_ylim())
    ax_loeuf_r.set_yticks(scale_pos)
    ax_loeuf_r.set_yticklabels([f'{v:.1f}' for v in scale_vals], fontsize=9)
    ax_loeuf_r.tick_params(axis='y', length=4, pad=3)
    for sp in ax_loeuf_r.spines.values():
        sp.set_visible(False)

    # Heatmap
    im = ax_h.imshow(heatmap_log, aspect='auto', cmap='magma',
                     vmin=0, vmax=np.percentile(heatmap_log, 99))
    ax_h.set_xticks(range(len(tissue_labels)))
    ax_h.set_xticklabels(tissue_labels, rotation=90, fontsize=11)
    ax_h.set_yticks([])
    ax_h.set_frame_on(False)

    # Expression colorbar
    cbar_x = left + names_width + gap + strip_width + scale_gap + heatmap_width + gap
    cbar_ax = fig.add_axes([cbar_x, bottom, strip_width, height])
    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.set_label('log2(TPM + 1)', fontsize=11, rotation=270, labelpad=20)

    fig.suptitle(
        f'GTEx expression heatmap — {n} network genes × {n_tissues_label}\n'
        f'(genes sorted by LOEUF; tissues hierarchically clustered)',
        fontsize=13, fontweight='bold', y=0.995)

    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT_DIR, f'{fname}.{ext}'),
                    dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved {fname}.png/pdf")


# Full heatmap (all 54 tissues)
draw_heatmap(TISSUES, 'figure_phase1_gtex_heatmap',
             f'{len(TISSUES)} tissues')

# Limited heatmap (curated tissue subset)
limited_tissues = [t for t in TISSUES if t not in EXCLUDE_TISSUES]
draw_heatmap(limited_tissues, 'figure_phase1_gtex_heatmap_limited',
             f'{len(limited_tissues)} tissues (curated)')

print("\nDone!")
