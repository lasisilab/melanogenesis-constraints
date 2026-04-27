"""
phase0_network_constraint.py

Phase 0: Network position vs. evolutionary constraint — three-panel figure
matching the layout of the Raghunath original (Tina's main-branch analysis)
but using the PEQG poster's 6 functional categories and color scheme.

Layout (1 × 3):
  Panel A: Betweenness centrality vs. LOEUF scatter
  Panel B: KEGG pathway count vs. LOEUF scatter  [requires data/kegg_pathway_counts.csv]
  Panel C: LOEUF by functional category (boxplot, median = white open circle)

Falls back to 1 × 2 (A + C) if KEGG data has not yet been fetched.

Statistics:
  - Spearman ρ(betweenness_centrality, LOEUF)
  - Spearman ρ(kegg_pathway_count, LOEUF)
  - Kruskal-Wallis + Dunn's posthoc (Bonferroni) on LOEUF by category
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import scikit_posthocs as sp
from adjustText import adjust_text

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_counts.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared constants (consistent with all Phase 1 figures) ────────────────
CATEGORY_ORDER = [
    'Pigment-specific',
    'Developmental/NC',
    'Generic signaling',
    'Cytokines/growth factors',
    'Apoptosis/cell death',
    'Other',
]

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

# Genes to label in scatter panels (matches original figure + additions)
LABEL_GENES = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MLANA', 'MC1R',
    'OCA2', 'MITF', 'SOX10', 'PAX3', 'TFAP2A',
    'KIT', 'EDNRB',
    'NFKB1', 'TP53', 'AKT1', 'MAPK1',
}

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading master dataset...")
df = pd.read_csv(MASTER_CSV)
df['gene'] = df['gene'].str.upper()
print(f"  {len(df)} genes loaded")

df_bc = df.dropna(subset=['betweenness_centrality', 'LOEUF']).copy()
print(f"  {len(df_bc)} genes with betweenness_centrality and LOEUF")

missing_bc = df.loc[df['betweenness_centrality'].isna(), 'gene'].tolist()
if missing_bc:
    print(f"  Missing betweenness_centrality: {missing_bc}")

# ── Load KEGG pathway counts (optional) ───────────────────────────────────
has_kegg = False
if os.path.exists(KEGG_CSV):
    print(f"\nLoading KEGG pathway counts...")
    kegg_df = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
    kegg_df['gene'] = kegg_df['gene'].str.upper()
    df_kegg = df.dropna(subset=['LOEUF']).merge(kegg_df, on='gene', how='inner')
    df_kegg = df_kegg.dropna(subset=['kegg_pathway_count'])
    print(f"  {len(df_kegg)} genes with KEGG pathway counts and LOEUF")
    has_kegg = True
else:
    df_kegg = None
    print(f"\nNote: {KEGG_CSV} not found — Panel B skipped.")
    print("       Run analysis/fetch_kegg_pathways.py to generate it.")

# ── Statistics: Betweenness centrality ────────────────────────────────────
print("\n=== Betweenness Centrality ===")
rho_bc, pval_bc = stats.spearmanr(df_bc['betweenness_centrality'], df_bc['LOEUF'])
print(f"  Spearman ρ = {rho_bc:.3f}, p = {pval_bc:.3e}")

print("\n  Per-category Spearman ρ(betweenness_centrality, LOEUF):")
for cat in CATEGORY_ORDER:
    sub = df_bc[df_bc['functional_category'] == cat]
    if len(sub) < 5:
        continue
    r, p = stats.spearmanr(sub['betweenness_centrality'], sub['LOEUF'])
    print(f"    {cat} (n={len(sub)}): ρ = {r:.3f}, p = {p:.3e}")

# ── Statistics: KEGG pathway count ────────────────────────────────────────
if has_kegg:
    print(f"\n=== KEGG Pathway Count (n={len(df_kegg)}) ===")
    rho_kegg, pval_kegg = stats.spearmanr(
        df_kegg['kegg_pathway_count'], df_kegg['LOEUF'])
    print(f"  Spearman ρ = {rho_kegg:.3f}, p = {pval_kegg:.3e}")

# ── Statistics: LOEUF by functional category ──────────────────────────────
print("\n=== LOEUF by Functional Category ===")
cats_present = [c for c in CATEGORY_ORDER if c in df_bc['functional_category'].values]
groups_loeuf = [df_bc.loc[df_bc['functional_category'] == c, 'LOEUF'].values
                for c in cats_present]
kw_loeuf, kw_p_loeuf = stats.kruskal(*[g for g in groups_loeuf if len(g) > 0])
print(f"  Kruskal-Wallis: H = {kw_loeuf:.2f}, p = {kw_p_loeuf:.3e}")

dunn_loeuf = sp.posthoc_dunn(
    df_bc, val_col='LOEUF', group_col='functional_category',
    p_adjust='bonferroni')
print("\n  Dunn's posthoc (Bonferroni), significant pairs:")
for i in range(len(cats_present)):
    for j in range(i + 1, len(cats_present)):
        p = dunn_loeuf.loc[cats_present[i], cats_present[j]]
        if p < 0.05:
            print(f"    {cats_present[i]} vs {cats_present[j]}: p = {p:.4f}")

print("\nLOEUF summary by functional category:")
summary = (df_bc.groupby('functional_category')['LOEUF']
           .agg(['count', 'median', 'mean'])
           .round(3)
           .rename(columns={'count': 'N', 'median': 'Median', 'mean': 'Mean'}))
summary = summary.loc[[c for c in CATEGORY_ORDER if c in summary.index]]
print(summary.to_string())


# ══════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════

def draw_bracket(ax, x1, x2, y, p_text, y_range, lw=1.0, fontsize=12):
    h = y_range * 0.025
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h],
            color='black', lw=lw, zorder=5)
    ax.text((x1 + x2) / 2, y + h * 0.3, p_text,
            ha='center', va='bottom', fontsize=fontsize, zorder=5)


def scatter_panel(ax, data, x_col, rho, pval, xlabel, title, panel_letter,
                  x_label_note=''):
    """Scatter of x_col vs. LOEUF, colored by functional category."""
    for cat in CATEGORY_ORDER:
        sub = data[data['functional_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub[x_col], sub['LOEUF'],
                   c=CATEGORY_COLORS[cat], s=55, alpha=0.85,
                   edgecolors='white', linewidths=0.4,
                   label=f'{cat} (n={len(sub)})', zorder=3)

    # OLS trend line
    slope, intercept, *_ = stats.linregress(data[x_col], data['LOEUF'])
    x_line = np.linspace(data[x_col].min(), data[x_col].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, '--',
            color='#555555', alpha=0.55, lw=1.2, zorder=1)

    # Gene labels
    texts = []
    for _, row in data.iterrows():
        if row['gene'] in LABEL_GENES:
            texts.append(ax.text(row[x_col], row['LOEUF'],
                                 row['gene'], fontsize=9, fontweight='bold',
                                 color='#222222', style='italic'))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.6),
                force_points=(2.0, 2.0), force_text=(2.0, 2.0),
                expand_points=(2.5, 2.5), expand_text=(1.8, 1.8),
                iter_lim=1000)

    # Spearman annotation box (inside lower-right)
    ax.text(0.97, 0.04,
            f'Spearman \u03c1 = {rho:.3f},  p = {pval:.3f}',
            transform=ax.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#999999', alpha=0.9))

    # Legend (upper right, inside)
    ax.legend(fontsize=8.5, loc='upper right',
              framealpha=0.9, edgecolor='gray',
              handletextpad=0.4, borderpad=0.5, markerscale=1.1)

    x_full = xlabel
    if x_label_note:
        x_full += f'\n{x_label_note}'
    ax.set_xlabel(x_full, fontsize=12)
    ax.set_ylabel('LOEUF (higher = less constrained)', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=11)
    ax.text(-0.10, 1.06, panel_letter, transform=ax.transAxes, fontsize=20,
            fontweight='bold', va='top')


def loeuf_boxplot_panel(ax, data, dunn, kw_p, title, panel_letter):
    """LOEUF distributions by functional category with white-circle medians."""
    bp_data, colors, ns = [], [], []
    for cat in CATEGORY_ORDER:
        vals = data.loc[data['functional_category'] == cat, 'LOEUF'].values
        bp_data.append(vals)
        colors.append(CATEGORY_COLORS[cat])
        ns.append(len(vals))

    positions = list(range(len(CATEGORY_ORDER)))

    # Draw boxplot with black median line
    bp = ax.boxplot(bp_data, positions=positions, widths=0.58,
                    patch_artist=True, showfliers=False,
                    zorder=2,
                    whiskerprops=dict(color='#666666', linewidth=1.2),
                    capprops=dict(color='#666666', linewidth=1.2),
                    boxprops=dict(linewidth=1.0),
                    medianprops=dict(color='black', linewidth=1.5))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.80)
        patch.set_edgecolor('#444444')

    # Individual data points (jittered)
    rng = np.random.default_rng(42)
    for i, (vals, cat) in enumerate(zip(bp_data, CATEGORY_ORDER)):
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.16, 0.16, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c='white', s=22, alpha=0.75,
                   edgecolors=CATEGORY_COLORS[cat], linewidths=0.7, zorder=3)


    # n= labels below boxes
    all_vals = [v for vals in bp_data for v in vals]
    y_min, y_max = min(all_vals), max(all_vals)
    y_range = y_max - y_min
    y_label_pos = y_min - y_range * 0.09
    for i, n in enumerate(ns):
        ax.text(i, y_label_pos, f'n={n}', ha='center', fontsize=9.5,
                color='#555555', fontstyle='italic', zorder=4)

    # Significance brackets (Dunn's posthoc)
    local_cats = [c for c in CATEGORY_ORDER if c in data['functional_category'].values]
    sig_pairs = []
    for i in range(len(local_cats)):
        for j in range(i + 1, len(local_cats)):
            c1, c2 = local_cats[i], local_cats[j]
            if c1 in dunn.index and c2 in dunn.columns:
                p = dunn.loc[c1, c2]
                if p < 0.05:
                    idx1 = CATEGORY_ORDER.index(c1)
                    idx2 = CATEGORY_ORDER.index(c2)
                    sig_pairs.append((idx1, idx2, p))
    sig_pairs.sort(key=lambda x: x[1] - x[0])

    bracket_y_start  = y_max + y_range * 0.08
    bracket_spacing  = y_range * 0.11
    for k, (i, j, p) in enumerate(sig_pairs):
        y_bracket = bracket_y_start + k * bracket_spacing
        p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
        draw_bracket(ax, i, j, y_bracket, p_text, y_range)

    top_y = bracket_y_start + max(len(sig_pairs), 1) * bracket_spacing
    ax.set_ylim(y_label_pos - abs(y_label_pos) * 0.1,
                top_y + bracket_spacing * 0.8)

    # KW annotation (top centre)
    ax.text(0.50, 0.98, f'Kruskal-Wallis p = {kw_p:.2e}',
            transform=ax.transAxes, fontsize=10, ha='center', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#999999', alpha=0.9))

    ax.set_xticks(positions)
    ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                       fontsize=10.5, rotation=30, ha='right')
    ax.set_ylabel('LOEUF (higher = less constrained)', fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=11)
    ax.text(-0.10, 1.06, panel_letter, transform=ax.transAxes, fontsize=20,
            fontweight='bold', va='top')


# ══════════════════════════════════════════════════════════════════════════
# Build figure
# ══════════════════════════════════════════════════════════════════════════
if has_kegg:
    fig, (ax_A, ax_B, ax_C) = plt.subplots(1, 3, figsize=(21, 7.5))
    fig.subplots_adjust(wspace=0.32, left=0.06, right=0.98,
                        top=0.88, bottom=0.25)
else:
    fig, (ax_A, ax_C) = plt.subplots(1, 2, figsize=(14, 7.5))
    fig.subplots_adjust(wspace=0.32, left=0.08, right=0.98,
                        top=0.88, bottom=0.25)

# ── Panel A: Betweenness centrality vs. LOEUF ─────────────────────────────
scatter_panel(
    ax_A, df_bc,
    x_col='betweenness_centrality',
    rho=rho_bc, pval=pval_bc,
    xlabel='Betweenness centrality',
    x_label_note='(within-pathway integration)',
    title='Network position vs.\nevolutionary constraint',
    panel_letter='A',
)

# ── Panel B: KEGG pathway count vs. LOEUF ─────────────────────────────────
if has_kegg:
    scatter_panel(
        ax_B, df_kegg,
        x_col='kegg_pathway_count',
        rho=rho_kegg, pval=pval_kegg,
        xlabel='Number of KEGG pathways',
        x_label_note='(cross-system connectivity)',
        title='Pathway involvement vs.\nevolutionary constraint',
        panel_letter='B',
    )

# ── Panel C: LOEUF by functional category ─────────────────────────────────
loeuf_boxplot_panel(
    ax_C, df_bc,
    dunn=dunn_loeuf, kw_p=kw_p_loeuf,
    title='Constraint by\nfunctional category',
    panel_letter='B' if not has_kegg else 'C',
)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ('png', 'pdf'):
    out_path = os.path.join(OUT_DIR, f'figure_phase0_network_constraint.{ext}')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {out_path}")
plt.close(fig)

print("\nDone!")
