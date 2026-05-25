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

# Master gene list — actual row order is computed per-heatmap via hierarchical
# clustering of the normalized expression profiles (so genes with similar
# tissue patterns cluster into modules).
gene_order = df['gene'].values

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

# Manual tissue groupings for the limited (poster) heatmap.
# Order within each group is preserved; groups are drawn left-to-right.
TISSUE_GROUPS = [
    ('Skin', [
        'Skin - Sun Exposed (Lower leg)',
        'Skin - Not Sun Exposed (Suprapubic)',
    ]),
    ('Neural / Pigment', [
        'Brain - Substantia nigra',
        'Cells - Cultured fibroblasts',
        'Brain - Cerebellum',
        'Brain - Cortex',
        'Brain - Frontal Cortex (BA9)',
        'Brain - Hippocampus',
        'Brain - Hypothalamus',
        'Brain - Amygdala',
        'Brain - Spinal cord (cervical c-1)',
        'Pituitary',
        'Nerve - Tibial',
    ]),
    ('Cardiovascular', [
        'Heart - Atrial Appendage',
        'Heart - Left Ventricle',
        'Artery - Aorta',
        'Artery - Coronary',
        'Artery - Tibial',
    ]),
    ('Endocrine / Metabolic', [
        'Thyroid',
        'Adrenal Gland',
        'Pancreas',
        'Liver',
        'Adipose - Subcutaneous',
        'Adipose - Visceral (Omentum)',
    ]),
    ('Respiratory / Excretory', [
        'Lung',
        'Kidney - Cortex',
        'Kidney - Medulla',
        'Bladder',
    ]),
    ('Female Reproductive', [
        'Ovary',
        'Uterus',
        'Vagina',
        'Cervix - Ectocervix',
        'Cervix - Endocervix',
        'Fallopian Tube',
        'Breast - Mammary Tissue',
    ]),
    ('Male Reproductive', [
        'Testis',
        'Prostate',
    ]),
    ('Digestive', [
        'Stomach',
        'Colon - Sigmoid',
        'Colon - Transverse',
    ]),
    ('Immune', [
        'Whole Blood',
        'Spleen',
    ]),
    ('Muscular', [
        'Muscle - Skeletal',
    ]),
]

BOLD_TISSUES = {
    'Skin - Sun Exposed (Lower leg)',
    'Skin - Not Sun Exposed (Suprapubic)',
}
BOX_TISSUES = {
    'Skin - Sun Exposed (Lower leg)',
    'Skin - Not Sun Exposed (Suprapubic)',
}


def draw_heatmap(tissue_subset, fname, n_tissues_label, groups=None):
    sub_expr = expr_df[tissue_subset].loc[gene_order]
    log_data = np.log2(sub_expr.values + 1)
    # Row-normalize: each cell = x_i / max_i(x). Range 0–1. This is the
    # "tau-style" expression view — tau = mean(1 - normalized) across tissues.
    row_max  = log_data.max(axis=1, keepdims=True)
    safe_max = np.where(row_max == 0, 1.0, row_max)
    norm_data = log_data / safe_max
    norm_data = np.where(row_max == 0, np.nan, norm_data)

    # Hierarchical clustering of rows (genes) on the normalized profiles, so
    # genes with similar tissue-expression patterns cluster into modules.
    g_link  = linkage(pdist(np.nan_to_num(norm_data), metric='correlation'),
                      method='average')
    g_order = leaves_list(g_link)
    norm_data    = norm_data[g_order]
    gene_labels  = gene_order[g_order]
    gene_sort_method = 'rows clustered by expression profile'

    if groups is None:
        # Hierarchical clustering of tissues (columns) — cluster on the
        # normalized data so column ordering reflects the same view shown.
        t_link  = linkage(pdist(np.nan_to_num(norm_data).T, metric='correlation'),
                          method='average')
        t_order = leaves_list(t_link)
        heatmap_data  = norm_data[:, t_order]
        tissue_labels = [tissue_subset[i] for i in t_order]
        group_boundaries = []   # (start_idx, end_idx, label)
        cluster_method   = 'hierarchically clustered'
    else:
        # Manual ordering by groups
        manual_order = []
        group_boundaries = []
        idx_lookup = {t: i for i, t in enumerate(tissue_subset)}
        cursor = 0
        for label, tissues in groups:
            present = [t for t in tissues if t in idx_lookup]
            if not present:
                continue
            start = cursor
            for t in present:
                manual_order.append(idx_lookup[t])
            cursor += len(present)
            group_boundaries.append((start, cursor, label))
        heatmap_data  = norm_data[:, manual_order]
        tissue_labels = [tissue_subset[i] for i in manual_order]
        cluster_method = 'grouped by system'

    n = len(gene_labels)
    fig_width = max(14, round(22 * len(tissue_labels) / 54))
    fig = plt.figure(figsize=(fig_width, 22))
    left            = 0.06
    bottom          = 0.13
    # leave a little more room at the top when category headers are drawn
    # (the limited heatmap staggers headers across two rows)
    height          = 0.66 if groups is not None else 0.77
    strip_width     = 0.013
    names_width     = 0.045   # gene names sit at the left edge now
    gap_names_heat  = 0.003
    heatmap_width   = 0.62
    gap_heat_cbar   = 0.008

    # Order from LEFT → RIGHT: gene names | heatmap | cbar
    ax_names = fig.add_axes([left, bottom, names_width, height])
    ax_h     = fig.add_axes([left + names_width + gap_names_heat,
                             bottom, heatmap_width, height])

    # Gene names — right-aligned to the heatmap edge
    ax_names.set_xlim(0, 1)
    ax_names.set_ylim(-0.5, n - 0.5)
    ax_names.invert_yaxis()
    ax_names.set_xticks([])
    ax_names.set_yticks([])
    for sp in ax_names.spines.values():
        sp.set_visible(False)
    for i, gene in enumerate(gene_labels):
        ax_names.text(1.0, i, gene, fontsize=7.5, ha='right', va='center')

    # Heatmap (row-normalized expression, 0–1)
    im = ax_h.imshow(heatmap_data, aspect='auto', cmap='magma',
                     vmin=0, vmax=1)
    ax_h.set_xticks(range(len(tissue_labels)))
    ax_h.set_xticklabels(tissue_labels, rotation=90, fontsize=11)
    ax_h.set_yticks([])
    ax_h.set_frame_on(False)

    # Bold tissue labels for skin + pigment-relevant tissues
    for ticklabel, tname in zip(ax_h.get_xticklabels(), tissue_labels):
        if tname in BOLD_TISSUES:
            ticklabel.set_fontweight('bold')

    # Group decorations: vertical separators, headers above heatmap, box around skin
    if group_boundaries:
        import matplotlib.patches as mpatches
        for start, end, label in group_boundaries:
            # Vertical separator at right edge of each group (except last)
            if end < len(tissue_labels):
                ax_h.axvline(end - 0.5, color='black', lw=0.8, alpha=0.55, zorder=4)
            # Category header — rotated 90° so all labels share the same
            # baseline regardless of group width.
            center = (start + end - 1) / 2
            ax_h.text(center, 1.14, label, ha='center', va='bottom',
                      rotation=90, rotation_mode='anchor',
                      fontsize=14, fontweight='bold', color='#222',
                      transform=ax_h.get_xaxis_transform())
            # Box around the Skin group (red, bold)
            if label == 'Skin':
                rect = mpatches.Rectangle(
                    (start - 0.5, -0.5), end - start, n,
                    linewidth=2.0, edgecolor='#c0282c', facecolor='none',
                    zorder=5, clip_on=False)
                ax_h.add_patch(rect)

    # Expression colorbar (row-normalized expression → τ component)
    cbar_x = (left + names_width + gap_names_heat + heatmap_width + gap_heat_cbar)
    cbar_ax = fig.add_axes([cbar_x, bottom, strip_width, height])
    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.set_label('Relative expression  (x / max per gene)\n'
                 'τ  =  mean(1 − this)  across tissues',
                 fontsize=12, fontweight='bold', rotation=270, labelpad=32)

    fig.suptitle(
        f'GTEx expression heatmap — {n} network genes × {n_tissues_label}\n'
        f'({gene_sort_method}; cells = τ-style row-normalized expression; '
        f'tissues {cluster_method})',
        fontsize=13, fontweight='bold', y=0.995)

    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT_DIR, f'{fname}.{ext}'),
                    dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved {fname}.png/pdf")


# Full heatmap (all 54 tissues)
draw_heatmap(TISSUES, 'figure_phase1_gtex_heatmap',
             f'{len(TISSUES)} tissues')

# Limited heatmap (curated tissue subset, grouped by system)
limited_tissues = [t for t in TISSUES if t not in EXCLUDE_TISSUES]
draw_heatmap(limited_tissues, 'figure_phase1_gtex_heatmap_limited',
             f'{len(limited_tissues)} tissues (curated)',
             groups=TISSUE_GROUPS)

# ===========================================================================
# Figure 4: Skin-specificity sorted pair (all 54 tissues, two orderings)
# ===========================================================================
# Two side-by-side heatmaps showing the same data but ordered by two different
# per-gene skin-specificity scores:
#   LEFT  — Option 2: mean(norm in skin) − mean(norm in non-skin)
#   RIGHT — Welch's t-statistic on log2(TPM+1), skin vs non-skin
# Same tissue grouping/labels/separators as the limited heatmap, but uses all
# 54 GTEx tissues. Keeps the curated limited heatmap (above) unchanged.

SKIN_SET = {
    'Skin - Sun Exposed (Lower leg)',
    'Skin - Not Sun Exposed (Suprapubic)',
}

# Full tissue groupings — covers ALL 54 GTEx v8 tissues, including the ones
# excluded from the curated limited heatmap.
TISSUE_GROUPS_FULL = [
    ('Skin', [
        'Skin - Sun Exposed (Lower leg)',
        'Skin - Not Sun Exposed (Suprapubic)',
    ]),
    ('Neural / Pigment', [
        'Brain - Substantia nigra',
        'Cells - Cultured fibroblasts',
        'Brain - Cerebellum',
        'Brain - Cerebellar Hemisphere',
        'Brain - Cortex',
        'Brain - Frontal Cortex (BA9)',
        'Brain - Anterior cingulate cortex (BA24)',
        'Brain - Hippocampus',
        'Brain - Amygdala',
        'Brain - Hypothalamus',
        'Brain - Caudate (basal ganglia)',
        'Brain - Nucleus accumbens (basal ganglia)',
        'Brain - Putamen (basal ganglia)',
        'Brain - Spinal cord (cervical c-1)',
        'Pituitary',
        'Nerve - Tibial',
    ]),
    ('Cardiovascular', [
        'Heart - Atrial Appendage',
        'Heart - Left Ventricle',
        'Artery - Aorta',
        'Artery - Coronary',
        'Artery - Tibial',
    ]),
    ('Endocrine / Metabolic', [
        'Thyroid',
        'Adrenal Gland',
        'Pancreas',
        'Liver',
        'Adipose - Subcutaneous',
        'Adipose - Visceral (Omentum)',
    ]),
    ('Respiratory / Excretory', [
        'Lung',
        'Kidney - Cortex',
        'Kidney - Medulla',
        'Bladder',
    ]),
    ('Female Reproductive', [
        'Ovary',
        'Uterus',
        'Vagina',
        'Cervix - Ectocervix',
        'Cervix - Endocervix',
        'Fallopian Tube',
        'Breast - Mammary Tissue',
    ]),
    ('Male Reproductive', [
        'Testis',
        'Prostate',
    ]),
    ('Digestive', [
        'Stomach',
        'Esophagus - Mucosa',
        'Esophagus - Gastroesophageal Junction',
        'Esophagus - Muscularis',
        'Small Intestine - Terminal Ileum',
        'Colon - Sigmoid',
        'Colon - Transverse',
        'Minor Salivary Gland',
    ]),
    ('Immune', [
        'Whole Blood',
        'Spleen',
        'Cells - EBV-transformed lymphocytes',
    ]),
    ('Muscular', [
        'Muscle - Skeletal',
    ]),
]


def make_skin_sorted_panels_figure(fname='figure_phase1_gtex_heatmap_skin_sorted',
                                   tau_threshold=0.7):
    """Three side-by-side heatmaps over all 54 tissues, all rows sharing the
    same τ-style normalized cells but ordered three different ways:
        LEFT   — Welch's t-statistic (skin vs non-skin), descending
        MIDDLE — τ-split: τ ≥ tau_threshold on top (sorted by Welch's t),
                 τ < tau_threshold below (also sorted by Welch's t); a
                 horizontal divider marks the split
        RIGHT  — Composite score: Welch's t × τ, descending
    """
    print(f"\nFigure 4: Skin-specificity sorted panels ({fname})...")
    tissue_subset = TISSUES
    groups = TISSUE_GROUPS_FULL

    # Per-gene τ aligned to gene_order
    tau_map = dict(zip(df['gene'].values, df['tau'].values))
    tau_values = np.array([tau_map[g] for g in gene_order])

    sub_expr = expr_df[tissue_subset].loc[gene_order]
    log_data = np.log2(sub_expr.values + 1)
    row_max  = log_data.max(axis=1, keepdims=True)
    safe_max = np.where(row_max == 0, 1.0, row_max)
    norm_data = log_data / safe_max
    norm_data = np.where(row_max == 0, np.nan, norm_data)

    skin_mask = np.array([t in SKIN_SET for t in tissue_subset])

    # Welch's t-statistic on log2(TPM+1), skin vs non-skin (per gene)
    with np.errstate(invalid='ignore', divide='ignore'):
        t_stat, _ = stats.ttest_ind(log_data[:, skin_mask],
                                    log_data[:, ~skin_mask],
                                    axis=1, equal_var=False)
    t_stat = np.where(np.isnan(t_stat), -np.inf, t_stat)

    # Ordering 1 — Welch's t descending
    welch_order = np.argsort(-t_stat)

    # Ordering 2 — τ split: high-τ block on top (sorted by Welch's t),
    #             low-τ block below (also sorted by Welch's t).
    high_tau_idx = np.where(tau_values >= tau_threshold)[0]
    low_tau_idx  = np.where(tau_values <  tau_threshold)[0]
    high_sorted  = high_tau_idx[np.argsort(-t_stat[high_tau_idx])]
    low_sorted   = low_tau_idx[np.argsort(-t_stat[low_tau_idx])]
    split_order   = np.concatenate([high_sorted, low_sorted])
    split_boundary = len(high_sorted)   # first row index of low-τ block

    # Ordering 3 — Composite = Welch's t × τ
    finite_t   = np.where(np.isfinite(t_stat), t_stat, np.nan)
    composite  = finite_t * tau_values
    composite_order = np.argsort(-np.where(np.isnan(composite), -np.inf, composite))

    # Column ordering by tissue group
    manual_order = []
    group_boundaries = []
    idx_lookup = {t: i for i, t in enumerate(tissue_subset)}
    cursor = 0
    for label, tissues in groups:
        present = [t for t in tissues if t in idx_lookup]
        if not present:
            continue
        start = cursor
        for t in present:
            manual_order.append(idx_lookup[t])
        cursor += len(present)
        group_boundaries.append((start, cursor, label))

    norm_grouped  = norm_data[:, manual_order]
    tissue_labels = [tissue_subset[i] for i in manual_order]
    n = len(gene_order)

    finite_t_vals = t_stat[np.isfinite(t_stat)]
    finite_comp   = composite[np.isfinite(composite)]

    panels = [
        {
            'order': welch_order,
            'title': "Welch's t-statistic\n(skin vs non-skin)",
            'subtitle': f"range: {finite_t_vals.min():+.1f} → {finite_t_vals.max():+.1f}",
            'divider': None,
        },
        {
            'order': split_order,
            'title': f"τ-split  (τ ≥ {tau_threshold} on top,  τ < {tau_threshold} below)\n"
                     "rows within each block sorted by Welch's t",
            'subtitle': f"{len(high_sorted)} high-τ  ·  {len(low_sorted)} low-τ",
            'divider': split_boundary,
        },
        {
            'order': composite_order,
            'title': "Composite:  Welch's t × τ\n(combines skin-enrichment and overall tissue-specificity)",
            'subtitle': f"range: {finite_comp.min():+.1f} → {finite_comp.max():+.1f}",
            'divider': None,
        },
    ]

    # ── Figure layout (figure-fraction coords) ─────────────────────────────
    fig = plt.figure(figsize=(36, 22))
    bottom        = 0.10
    height        = 0.62
    left_margin   = 0.022
    names_w       = 0.026
    gap_nh        = 0.002
    heatmap_w     = 0.270
    panel_gap     = 0.022
    gap_cbar      = 0.006
    cbar_w        = 0.008

    import matplotlib.patches as mpatches

    def _draw_panel(ax_h, ax_names, data, genes, divider, panel_idx):
        im = ax_h.imshow(data, aspect='auto', cmap='magma', vmin=0, vmax=1)
        ax_h.set_xticks(range(len(tissue_labels)))
        ax_h.set_xticklabels(tissue_labels, rotation=90, fontsize=8)
        ax_h.set_yticks([])
        ax_h.set_frame_on(False)
        for ticklabel, tname in zip(ax_h.get_xticklabels(), tissue_labels):
            if tname in BOLD_TISSUES:
                ticklabel.set_fontweight('bold')
        for start, end, label in group_boundaries:
            if end < len(tissue_labels):
                ax_h.axvline(end - 0.5, color='black', lw=0.8, alpha=0.55,
                             zorder=4)
            center = (start + end - 1) / 2
            ax_h.text(center, 1.14, label, ha='center', va='bottom',
                      rotation=90, rotation_mode='anchor',
                      fontsize=12, fontweight='bold', color='#222',
                      transform=ax_h.get_xaxis_transform())
            if label == 'Skin':
                rect = mpatches.Rectangle(
                    (start - 0.5, -0.5), end - start, n,
                    linewidth=2.0, edgecolor='#c0282c', facecolor='none',
                    zorder=5, clip_on=False)
                ax_h.add_patch(rect)
        # Horizontal divider for τ-split panel
        if divider is not None:
            ax_h.axhline(divider - 0.5, color='#1f6feb', lw=2.0, zorder=6)
            ax_h.text(-0.5, divider - 0.5, ' τ-split',
                      ha='left', va='center', fontsize=10, fontweight='bold',
                      color='#1f6feb', zorder=6,
                      transform=ax_h.get_yaxis_transform())
        # Gene names (on the left of each panel)
        ax_names.set_xlim(0, 1)
        ax_names.set_ylim(-0.5, n - 0.5)
        ax_names.invert_yaxis()
        ax_names.set_xticks([])
        ax_names.set_yticks([])
        for sp in ax_names.spines.values():
            sp.set_visible(False)
        for i, gene in enumerate(genes):
            ax_names.text(1.0, i, gene, fontsize=7, ha='right', va='center')
        return im

    im = None
    panel_centers = []
    cursor_x = left_margin
    for idx, panel in enumerate(panels):
        ax_names = fig.add_axes([cursor_x, bottom, names_w, height])
        ax_h     = fig.add_axes([cursor_x + names_w + gap_nh, bottom,
                                 heatmap_w, height])
        data  = norm_grouped[panel['order']]
        genes = gene_order[panel['order']]
        im_now = _draw_panel(ax_h, ax_names, data, genes,
                             panel['divider'], idx)
        if im is None:
            im = im_now
        panel_centers.append(cursor_x + names_w + gap_nh + heatmap_w / 2)
        cursor_x += names_w + gap_nh + heatmap_w + panel_gap

    cbar_x  = cursor_x - panel_gap + gap_cbar
    cbar_ax = fig.add_axes([cbar_x, bottom, cbar_w, height])

    # Panel titles ABOVE the rotated group labels
    panel_title_y    = 0.945
    panel_subtitle_y = 0.915
    for cx, panel in zip(panel_centers, panels):
        fig.text(cx, panel_title_y, panel['title'],
                 ha='center', va='bottom', fontsize=14, fontweight='bold',
                 linespacing=1.15)
        fig.text(cx, panel_subtitle_y, panel['subtitle'],
                 ha='center', va='top', fontsize=11, color='#555')

    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.set_label('Relative expression  (x / max per gene)\n'
                 'τ  =  mean(1 − this)  across tissues',
                 fontsize=12, fontweight='bold', rotation=270, labelpad=32)

    fig.suptitle(
        f'Skin-specificity orderings — {n} network genes × '
        f'{len(tissue_labels)} GTEx tissues\n'
        f'(cells = τ-style row-normalized expression; tissues grouped by system)',
        fontsize=16, fontweight='bold', y=0.985)

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'{fname}.{ext}')
        fig.savefig(path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved {fname}.png/pdf")


make_skin_sorted_panels_figure()


def make_skin_split_only_figure(fname='figure_phase1_gtex_heatmap_skin_split',
                                tau_threshold=0.7):
    """Standalone single-panel version of the τ-split heatmap: high-τ block
    (tissue-specific) on top sorted by Welch's t, low-τ block (broadly
    expressed) below also sorted by Welch's t, with a horizontal divider
    marking the τ cutoff."""
    print(f"\nFigure 5: τ-split heatmap only ({fname})...")
    tissue_subset = TISSUES
    groups = TISSUE_GROUPS_FULL

    tau_map = dict(zip(df['gene'].values, df['tau'].values))
    tau_values = np.array([tau_map[g] for g in gene_order])

    sub_expr = expr_df[tissue_subset].loc[gene_order]
    log_data = np.log2(sub_expr.values + 1)
    row_max  = log_data.max(axis=1, keepdims=True)
    safe_max = np.where(row_max == 0, 1.0, row_max)
    norm_data = log_data / safe_max
    norm_data = np.where(row_max == 0, np.nan, norm_data)

    skin_mask = np.array([t in SKIN_SET for t in tissue_subset])
    with np.errstate(invalid='ignore', divide='ignore'):
        t_stat, _ = stats.ttest_ind(log_data[:, skin_mask],
                                    log_data[:, ~skin_mask],
                                    axis=1, equal_var=False)
    t_stat = np.where(np.isnan(t_stat), -np.inf, t_stat)

    high_tau_idx = np.where(tau_values >= tau_threshold)[0]
    low_tau_idx  = np.where(tau_values <  tau_threshold)[0]
    high_sorted  = high_tau_idx[np.argsort(-t_stat[high_tau_idx])]
    low_sorted   = low_tau_idx[np.argsort(-t_stat[low_tau_idx])]
    split_order  = np.concatenate([high_sorted, low_sorted])
    split_boundary = len(high_sorted)

    # Column ordering
    manual_order = []
    group_boundaries = []
    idx_lookup = {t: i for i, t in enumerate(tissue_subset)}
    cursor = 0
    for label, tissues in groups:
        present = [t for t in tissues if t in idx_lookup]
        if not present:
            continue
        start = cursor
        for t in present:
            manual_order.append(idx_lookup[t])
        cursor += len(present)
        group_boundaries.append((start, cursor, label))

    norm_grouped  = norm_data[:, manual_order]
    tissue_labels = [tissue_subset[i] for i in manual_order]
    n = len(gene_order)

    data  = norm_grouped[split_order]
    genes = gene_order[split_order]

    # Layout: gene names | heatmap | cbar
    fig = plt.figure(figsize=(16, 22))
    bottom        = 0.10
    height        = 0.62
    left          = 0.06
    names_w       = 0.045
    gap_nh        = 0.003
    heatmap_w     = 0.78
    gap_cbar      = 0.008
    cbar_w        = 0.012

    ax_names = fig.add_axes([left, bottom, names_w, height])
    ax_h     = fig.add_axes([left + names_w + gap_nh, bottom,
                             heatmap_w, height])
    cbar_x   = left + names_w + gap_nh + heatmap_w + gap_cbar
    cbar_ax  = fig.add_axes([cbar_x, bottom, cbar_w, height])

    import matplotlib.patches as mpatches

    im = ax_h.imshow(data, aspect='auto', cmap='magma', vmin=0, vmax=1)
    ax_h.set_xticks([])
    ax_h.set_yticks([])
    ax_h.set_frame_on(False)

    for start, end, label in group_boundaries:
        if end < len(tissue_labels):
            ax_h.axvline(end - 0.5, color='black', lw=0.8, alpha=0.55, zorder=4)
        center = (start + end - 1) / 2
        # Category label below the heatmap (replaces per-tissue names),
        # rotated to read top-to-bottom (anchor at top, text extends down).
        ax_h.text(center, -0.02, label, ha='right', va='center',
                  rotation=90, rotation_mode='anchor',
                  fontsize=14, fontweight='bold', color='#222',
                  transform=ax_h.get_xaxis_transform())
        if label == 'Skin':
            rect = mpatches.Rectangle(
                (start - 0.5, -0.5), end - start, n,
                linewidth=2.0, edgecolor='#c0282c', facecolor='none',
                zorder=5, clip_on=False)
            ax_h.add_patch(rect)

    # τ-split divider line
    ax_h.axhline(split_boundary - 0.5, color='#1f6feb', lw=2.5, zorder=6)

    # Gene names
    ax_names.set_xlim(0, 1)
    ax_names.set_ylim(-0.5, n - 0.5)
    ax_names.invert_yaxis()
    ax_names.set_xticks([])
    ax_names.set_yticks([])
    for sp in ax_names.spines.values():
        sp.set_visible(False)
    for i, gene in enumerate(genes):
        ax_names.text(1.0, i, gene, fontsize=8, ha='right', va='center')

    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.set_label('Relative expression  (x / max per gene)\n'
                 'τ  =  mean(1 − this)  across tissues',
                 fontsize=12, fontweight='bold', rotation=270, labelpad=32)

    fig.suptitle(
        f'τ-split skin-specificity heatmap — {n} network genes × '
        f'{len(tissue_labels)} GTEx tissues\n'
        f'tissue-specific (τ ≥ {tau_threshold}) on top, broadly expressed below; '
        "within each block rows are sorted by Welch's t (skin vs non-skin)",
        fontsize=14, fontweight='bold', y=0.985)

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'{fname}.{ext}')
        fig.savefig(path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved {fname}.png/pdf")


make_skin_split_only_figure()


# ===========================================================================
# Figure 6: τ-split heatmap with selection-quadrant side bar (overlay)
# ===========================================================================
# Overlays the population-specific-selection axis from phase2_pop_difference.py
# directly onto the τ-split heatmap. Each gene row gets a colored cell in a
# slim column between the gene names and the heatmap, with the same purple /
# yellow-green / teal / grey palette as the PBS quadrant scatter. Lets the
# reader see at a glance whether Melanesian-specific selection (yellow-green)
# concentrates in the upper, tissue-restricted block and African-specific
# selection (purple) concentrates in the lower, broadly-expressed block.
COLOR_AFR  = '#6600cc'    # purple — African-specific
COLOR_MEL  = '#ccff00'    # yellow-green — Melanesian-specific
COLOR_BOTH = '#00daa7'    # teal — both
COLOR_NONE = '#cccccc'    # light gray — neither
QUAD_COLORS = {'African': COLOR_AFR, 'Melanesian': COLOR_MEL,
               'Both': COLOR_BOTH, 'Neither': COLOR_NONE,
               'No PBS': '#ffffff'}


def make_skin_split_with_selection_figure(
        fname='figure_phase1_gtex_heatmap_skin_split_selection',
        tau_threshold=0.7):
    """τ-split heatmap with a selection-quadrant side bar to overlay the
    PBS-quadrant analysis on the tissue-specificity heatmap."""
    print(f"\nFigure 6: τ-split heatmap with selection overlay ({fname})...")
    tissue_subset = TISSUES
    groups = TISSUE_GROUPS_FULL

    # Per-gene τ
    tau_map = dict(zip(df['gene'].values, df['tau'].values))
    tau_values = np.array([tau_map[g] for g in gene_order])

    # Load PBS data and assign quadrants (same logic as phase2_pop_difference.py)
    pbs_csv = os.path.join(PROJECT_DIR, 'data', 'pbs_per_gene.csv')
    pbs = pd.read_csv(pbs_csv)[['gene', 'pbs1_african', 'pbs3_melanesian']]
    pbs['gene'] = pbs['gene'].str.upper()
    pbs_for_genes = (pd.DataFrame({'gene': gene_order})
                     .merge(pbs, on='gene', how='left'))
    thr_a = pbs_for_genes['pbs1_african'].quantile(0.75)
    thr_m = pbs_for_genes['pbs3_melanesian'].quantile(0.75)

    def _q(row):
        if pd.isna(row['pbs1_african']) or pd.isna(row['pbs3_melanesian']):
            return 'No PBS'
        hi_a = row['pbs1_african']    >= thr_a
        hi_m = row['pbs3_melanesian'] >= thr_m
        if hi_a and hi_m: return 'Both'
        if hi_a:          return 'African'
        if hi_m:          return 'Melanesian'
        return 'Neither'

    quadrants = pbs_for_genes.apply(_q, axis=1).values

    # ── Stats: contingency between τ-block and selection quadrant ────────
    tau_block = np.where(tau_values >= tau_threshold, 'high-τ', 'low-τ')

    print(f"\n  Contingency (τ-block × selection quadrant), all genes:")
    contingency_full = pd.crosstab(pd.Series(tau_block, name='τ-block'),
                                   pd.Series(quadrants, name='quadrant'))
    print(contingency_full.to_string())

    afr_mel_mask = np.isin(quadrants, ['African', 'Melanesian'])
    table_2x2 = pd.crosstab(pd.Series(tau_block[afr_mel_mask], name='τ-block'),
                            pd.Series(quadrants[afr_mel_mask], name='quadrant'))
    # Force consistent row/column order for the test
    for col in ['African', 'Melanesian']:
        if col not in table_2x2.columns:
            table_2x2[col] = 0
    for row in ['high-τ', 'low-τ']:
        if row not in table_2x2.index:
            table_2x2.loc[row] = 0
    table_2x2 = table_2x2.loc[['high-τ', 'low-τ'], ['African', 'Melanesian']]
    print("\n  Restricted to African-specific and Melanesian-specific only:")
    print(table_2x2.to_string())
    odds, p_fisher = stats.fisher_exact(table_2x2.values)
    print(f"  Fisher's exact:  OR = {odds:.3f},  p = {p_fisher:.3e}")

    # ── Expression data + Welch's t for row ordering (mirrors split-only) ─
    sub_expr = expr_df[tissue_subset].loc[gene_order]
    log_data = np.log2(sub_expr.values + 1)
    row_max  = log_data.max(axis=1, keepdims=True)
    safe_max = np.where(row_max == 0, 1.0, row_max)
    norm_data = log_data / safe_max
    norm_data = np.where(row_max == 0, np.nan, norm_data)

    skin_mask = np.array([t in SKIN_SET for t in tissue_subset])
    with np.errstate(invalid='ignore', divide='ignore'):
        t_stat, _ = stats.ttest_ind(log_data[:, skin_mask],
                                    log_data[:, ~skin_mask],
                                    axis=1, equal_var=False)
    t_stat = np.where(np.isnan(t_stat), -np.inf, t_stat)

    high_tau_idx = np.where(tau_values >= tau_threshold)[0]
    low_tau_idx  = np.where(tau_values <  tau_threshold)[0]
    high_sorted  = high_tau_idx[np.argsort(-t_stat[high_tau_idx])]
    low_sorted   = low_tau_idx[np.argsort(-t_stat[low_tau_idx])]
    split_order  = np.concatenate([high_sorted, low_sorted])
    split_boundary = len(high_sorted)

    manual_order = []
    group_boundaries = []
    idx_lookup = {t: i for i, t in enumerate(tissue_subset)}
    cursor = 0
    for label, tissues in groups:
        present = [t for t in tissues if t in idx_lookup]
        if not present:
            continue
        start = cursor
        for t in present:
            manual_order.append(idx_lookup[t])
        cursor += len(present)
        group_boundaries.append((start, cursor, label))

    norm_grouped  = norm_data[:, manual_order]
    tissue_labels = [tissue_subset[i] for i in manual_order]
    n = len(gene_order)

    data       = norm_grouped[split_order]
    genes      = gene_order[split_order]
    quad_order = quadrants[split_order]

    # ── Layout: gene names | quadrant bar | heatmap | cbar ───────────────
    fig = plt.figure(figsize=(16, 22))
    bottom    = 0.10
    height    = 0.62
    left      = 0.04
    names_w   = 0.045
    gap_nq    = 0.002
    quad_w    = 0.014
    gap_qh    = 0.005
    heatmap_w = 0.74
    gap_cbar  = 0.008
    cbar_w    = 0.012

    ax_names = fig.add_axes([left, bottom, names_w, height])
    ax_quad  = fig.add_axes([left + names_w + gap_nq, bottom,
                             quad_w, height])
    ax_h     = fig.add_axes([left + names_w + gap_nq + quad_w + gap_qh,
                             bottom, heatmap_w, height])
    cbar_x   = (left + names_w + gap_nq + quad_w + gap_qh
                + heatmap_w + gap_cbar)
    cbar_ax  = fig.add_axes([cbar_x, bottom, cbar_w, height])

    import matplotlib.patches as mpatches

    im = ax_h.imshow(data, aspect='auto', cmap='magma', vmin=0, vmax=1)
    ax_h.set_xticks([])
    ax_h.set_yticks([])
    ax_h.set_frame_on(False)
    for start, end, label in group_boundaries:
        if end < len(tissue_labels):
            ax_h.axvline(end - 0.5, color='black', lw=0.8, alpha=0.55, zorder=4)
        center = (start + end - 1) / 2
        ax_h.text(center, -0.02, label, ha='right', va='center',
                  rotation=90, rotation_mode='anchor',
                  fontsize=14, fontweight='bold', color='#222',
                  transform=ax_h.get_xaxis_transform())
        if label == 'Skin':
            rect = mpatches.Rectangle(
                (start - 0.5, -0.5), end - start, n,
                linewidth=2.0, edgecolor='#c0282c', facecolor='none',
                zorder=5, clip_on=False)
            ax_h.add_patch(rect)
    # τ-split divider on the heatmap
    ax_h.axhline(split_boundary - 0.5, color='#1f6feb', lw=2.5, zorder=6)

    # Quadrant side-bar — one colored cell per gene row
    ax_quad.set_xlim(0, 1)
    ax_quad.set_ylim(-0.5, n - 0.5)
    ax_quad.invert_yaxis()
    ax_quad.set_xticks([])
    ax_quad.set_yticks([])
    for sp in ax_quad.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
        sp.set_color('#888')
    for i, q in enumerate(quad_order):
        ax_quad.add_patch(plt.Rectangle((0, i - 0.5), 1, 1,
                                        facecolor=QUAD_COLORS[q],
                                        edgecolor='none'))
    ax_quad.axhline(split_boundary - 0.5, color='#1f6feb', lw=2.0, zorder=6)
    ax_quad.set_title('Sel.', fontsize=10, fontweight='bold', pad=5)

    # Gene names
    ax_names.set_xlim(0, 1)
    ax_names.set_ylim(-0.5, n - 0.5)
    ax_names.invert_yaxis()
    ax_names.set_xticks([])
    ax_names.set_yticks([])
    for sp in ax_names.spines.values():
        sp.set_visible(False)
    for i, gene in enumerate(genes):
        ax_names.text(1.0, i, gene, fontsize=8, ha='right', va='center')

    # Expression colorbar
    cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
    cb.set_label('Relative expression  (x / max per gene)\n'
                 'τ  =  mean(1 − this)  across tissues',
                 fontsize=12, fontweight='bold', rotation=270, labelpad=32)

    # Quadrant legend (top-right of figure)
    legend_handles = [
        mpatches.Patch(facecolor=COLOR_MEL,  edgecolor='black', lw=0.6,
                       label=f'Melanesian-specific  '
                             f'(n={int((quadrants == "Melanesian").sum())})'),
        mpatches.Patch(facecolor=COLOR_AFR,  edgecolor='black', lw=0.6,
                       label=f'African-specific  '
                             f'(n={int((quadrants == "African").sum())})'),
        mpatches.Patch(facecolor=COLOR_BOTH, edgecolor='black', lw=0.6,
                       label=f'Both  (n={int((quadrants == "Both").sum())})'),
        mpatches.Patch(facecolor=COLOR_NONE, edgecolor='black', lw=0.6,
                       label=f'Neither  (n={int((quadrants == "Neither").sum())})'),
        mpatches.Patch(facecolor='#ffffff',  edgecolor='black', lw=0.6,
                       label=f'No PBS data  '
                             f'(n={int((quadrants == "No PBS").sum())})'),
    ]
    fig.legend(handles=legend_handles, loc='upper right',
               bbox_to_anchor=(0.99, 0.96), frameon=True, fontsize=10,
               title='PBS selection quadrant', title_fontsize=11)

    # Fisher result inset (top-left, mirrors the legend on the right)
    fisher_txt = (f"τ-block × selection (Afr vs. Mel only)\n"
                  f"OR = {odds:.2f}   ·   Fisher p = {p_fisher:.2e}\n"
                  f"high-τ:  Afr = {table_2x2.loc['high-τ','African']},  "
                  f"Mel = {table_2x2.loc['high-τ','Melanesian']}\n"
                  f"low-τ:   Afr = {table_2x2.loc['low-τ','African']},  "
                  f"Mel = {table_2x2.loc['low-τ','Melanesian']}")
    fig.text(0.04, 0.96, fisher_txt, ha='left', va='top',
             fontsize=10, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f8f8',
                       edgecolor='#888', lw=0.6))

    fig.suptitle(
        f'τ-split skin-specificity heatmap with selection overlay — '
        f'{n} network genes × {len(tissue_labels)} GTEx tissues\n'
        f'tissue-specific (τ ≥ {tau_threshold}) on top, broadly expressed '
        "below; coloured side-bar = PBS quadrant",
        fontsize=14, fontweight='bold', y=0.985)

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'{fname}.{ext}')
        fig.savefig(path, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved {fname}.png/pdf")


make_skin_split_with_selection_figure()

print("\nDone!")
