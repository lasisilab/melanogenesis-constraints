"""
Generate preliminary results figure for R35 grant proposal:
Melanogenesis network constraint and evolutionary selection.

Composite Figure: 3-panel layout
  A) Betweenness centrality vs. LOEUF (colored by functional category)
  B) KEGG pathway count vs. LOEUF (colored by functional category)
  C) LOEUF by functional category (boxplot with Dunn's posthoc)

All panels share a unified color scheme based on 5 functional categories
and a common y-axis scale (LOEUF 0–3.0).

Data sources:
  - network_constraint_data.csv: 130 genes from Raghunath et al. (2015)
    melanocyte network (134 genes minus 4 encoding cAMP which lacks
    a gnomAD LOEUF score).
  - LOEUF scores: gnomAD v2.1.1 (Karczewski et al. 2020)
  - Betweenness centrality: computed from Raghunath et al. (2015)
    directed network (429 edges)
  - KEGG pathway counts: queried from KEGG REST API (hsa pathways)

Functional category derivation:
  Categories are based on KEGG melanogenesis pathway (hsa04916)
  membership, KEGG pathway count, and known gene function:

  1. "Pigmentation enzymes" (n=7): Genes whose protein products
     directly participate in melanin biosynthesis or melanosome
     structure/function. All are members of the KEGG melanogenesis
     pathway (hsa04916) and function primarily in pigmentation.
     Genes: TYR, TYRP1, DCT (melanin biosynthesis enzymes),
            OCA2 (melanosome membrane transporter),
            PMEL (melanosome structural protein),
            MLANA (melanosome component/melanocyte marker),
            MC1R (melanocortin-1 receptor, pigmentation signaling)

  2. "Melanocyte regulatory" (n=5): Transcription factors and
     receptors that specify melanocyte fate from neural crest
     progenitors. These have dual roles in melanocyte biology AND
     broader neural crest/developmental processes, explaining their
     stronger evolutionary constraint.
     Genes: MITF (master melanocyte TF; also in KEGG melanogenesis),
            EDNRB (endothelin receptor B; also in KEGG melanogenesis),
            SOX10 (neural crest TF, melanocyte specification),
            PAX3 (neural crest TF, melanocyte development),
            TFAP2A (neural crest TF)
     SOX10, PAX3, and TFAP2A are not in the KEGG melanogenesis pathway
     but are well-established melanocyte specification factors
     (Baxter et al. 2004; Seberg et al. 2017).

  3. "Shared melanogenesis signaling" (n=17): Genes that appear in
     the KEGG melanogenesis pathway (hsa04916) BUT are generic
     signaling molecules present in many other KEGG pathways
     (typically >8 pathways). These are shared signaling
     infrastructure, not melanocyte-specific.
     Examples: MAPK1 (119 pathways), RAF1 (83), MAP2K1 (91),
              CREB1 (43), CTNNB1 (33), etc.

  4. "Pleiotropic signaling" (n=58): Genes NOT in KEGG melanogenesis
     but present in multiple KEGG pathways. These are broadly
     pleiotropic signaling genes included in the Raghunath network
     through their interactions with melanogenesis components.
     Examples: AKT1, NFKB1, PIK3CA, TP53, STAT3, etc.

  5. "Specialized/other" (n=43): Genes NOT in KEGG melanogenesis,
     typically in few KEGG pathways. These have specialized functions
     in apoptosis, sphingolipid metabolism, coagulation, prostaglandin
     signaling, or other processes connected to the melanocyte network.
     Examples: CFLAR, MCL1, SMPD1, PTGER1, etc.

Author: Generated with Claude for Tina Lasisi
Date: February 2026
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from scipy import stats
from adjustText import adjust_text
import scikit_posthocs as sp
import os

# ── Configuration ──
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FILE = os.path.join(PROJECT_DIR, 'output', 'network_constraint_data.csv')
OUTPUT_DIR = os.path.join(PROJECT_DIR, 'output')

# ── Load data ──
df = pd.read_csv(DATA_FILE)

# ══════════════════════════════════════════════════════════════
# FUNCTIONAL CATEGORY ASSIGNMENT
# ══════════════════════════════════════════════════════════════

PIGMENTATION_ENZYMES = {'TYR', 'TYRP1', 'DCT', 'OCA2', 'PMEL', 'MLANA', 'MC1R'}
MELANOCYTE_REGULATORY = {'MITF', 'EDNRB', 'SOX10', 'PAX3', 'TFAP2A'}


def assign_refined_category(row):
    """
    Assign refined functional category based on KEGG pathway membership
    and known gene function.
    """
    gene = row['gene']
    if gene in PIGMENTATION_ENZYMES:
        return 'Pigmentation enzymes'
    elif gene in MELANOCYTE_REGULATORY:
        return 'Melanocyte regulatory'
    elif row['kegg_functional_category'] == 'Shared melanogenesis signaling':
        return 'Shared melanogenesis signaling'
    elif row['kegg_functional_category'] == 'Pleiotropic signaling':
        return 'Pleiotropic signaling'
    else:
        return 'Specialized/other'


df['refined_category'] = df.apply(assign_refined_category, axis=1)

print("Refined functional category counts:")
print(df['refined_category'].value_counts().sort_index())
print()

# ══════════════════════════════════════════════════════════════
# UNIFIED COLOR SCHEME (used across ALL panels)
# ══════════════════════════════════════════════════════════════

CATEGORY_ORDER = [
    'Pigmentation enzymes',
    'Melanocyte regulatory',
    'Shared melanogenesis signaling',
    'Pleiotropic signaling',
    'Specialized/other',
]

CATEGORY_COLORS = {
    'Pigmentation enzymes':            '#D94040',   # Red
    'Melanocyte regulatory':           '#E8907E',   # Salmon
    'Shared melanogenesis signaling':  '#F5C242',   # Gold/amber
    'Pleiotropic signaling':           '#4878CF',   # Blue
    'Specialized/other':               '#B0B0B0',   # Gray
}

CATEGORY_SHORT_LABELS = {
    'Pigmentation enzymes':            'Pigmentation\nenzymes',
    'Melanocyte regulatory':           'Melanocyte\nregulatory',
    'Shared melanogenesis signaling':  'Shared melan.\nsignaling',
    'Pleiotropic signaling':           'Pleiotropic\nsignaling',
    'Specialized/other':               'Specialized/\nother',
}

# Common y-axis limits
YLIM = (-0.1, 3.0)

# ══════════════════════════════════════════════════════════════
# PLOTTING FUNCTIONS
# ══════════════════════════════════════════════════════════════

def plot_scatter(ax, x_col, y_col, plot_df, xlabel, ylabel, title,
                 label_genes, legend_loc='upper right'):
    """
    Scatter plot colored by refined functional category.
    All genes belong to one of 5 categories — no uncategorized dots.
    """
    # Plot each category in order so legend matches
    for cat in CATEGORY_ORDER:
        sub = plot_df[plot_df['refined_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub[x_col], sub[y_col],
                   c=CATEGORY_COLORS[cat], s=65, alpha=0.85,
                   edgecolors='white', linewidths=0.4,
                   label=f'{cat} (n={len(sub)})', zorder=3)

    # Trend line
    valid = plot_df.dropna(subset=[x_col, y_col])
    if x_col == 'n_kegg_pathways':
        valid = valid[valid[x_col] > 0]
    slope, intercept, _, _, _ = stats.linregress(valid[x_col], valid[y_col])
    x_line = np.linspace(valid[x_col].min(), valid[x_col].max(), 100)
    ax.plot(x_line, slope * x_line + intercept, '--',
            color='#555555', alpha=0.6, lw=1, zorder=1)

    # Gene labels (bold italic)
    texts = []
    for _, row in plot_df.iterrows():
        if (row['gene'] in label_genes
                and not pd.isna(row[x_col])
                and not pd.isna(row[y_col])):
            texts.append(ax.text(
                row[x_col], row[y_col], row['gene'],
                fontsize=16, fontweight='bold', color='#333333',
                style='italic'))

    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7),
                force_points=(2.5, 2.5), force_text=(2.5, 2.5),
                expand_points=(3.0, 3.0), expand_text=(2.0, 2.0),
                iter_lim=1000)

    # Spearman correlation — positioned below the legend
    rho, pval = stats.spearmanr(valid[x_col], valid[y_col])
    ax.text(0.97, 0.58, f'Spearman ρ = {rho:.3f}, p = {pval:.3f}',
            transform=ax.transAxes, fontsize=16, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.9))

    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=19, fontweight='bold', loc='left', pad=14)
    ax.legend(fontsize=16, loc=legend_loc, framealpha=0.9,
              edgecolor='gray', handletextpad=0.4, borderpad=0.5,
              markerscale=1.5)
    ax.tick_params(labelsize=18)
    ax.set_ylim(YLIM)


def draw_bracket(ax, x1, x2, y, p_text, lw=1.0, fontsize=18):
    """Draw significance bracket between positions x1 and x2."""
    h = 0.04
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h],
            color='black', lw=lw, zorder=5)
    ax.text((x1 + x2) / 2, y + 0.02, p_text,
            ha='center', va='bottom', fontsize=fontsize, zorder=5)


def plot_boxplot_with_pairwise(ax, plot_df, title):
    """
    Boxplot of LOEUF by refined functional category with
    Dunn's posthoc pairwise significance brackets.
    """
    bp_data, bp_color_list, ns = [], [], []
    for cat in CATEGORY_ORDER:
        vals = plot_df[plot_df['refined_category'] == cat]['LOEUF'].dropna()
        bp_data.append(vals.values)
        bp_color_list.append(CATEGORY_COLORS[cat])
        ns.append(len(vals))

    positions = list(range(len(CATEGORY_ORDER)))
    bp = ax.boxplot(bp_data, positions=positions, widths=0.6,
                    patch_artist=True, showfliers=False, zorder=2,
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(color='gray'),
                    capprops=dict(color='gray'))
    for patch, color in zip(bp['boxes'], bp_color_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)
        patch.set_edgecolor('gray')

    # Jittered individual points
    for i, (cat, vals) in enumerate(zip(CATEGORY_ORDER, bp_data)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c=CATEGORY_COLORS[cat], s=25, alpha=0.8,
                   edgecolors='white', linewidths=0.3, zorder=3)

    for i, n in enumerate(ns):
        ax.text(i, 0.03, f'n={n}', ha='center', fontsize=15, color='#555555',
                fontstyle='italic', zorder=4)

    # Kruskal-Wallis
    valid_data = [d for d in bp_data if len(d) > 0]
    _, kw_p = stats.kruskal(*valid_data)

    # Dunn's posthoc (Bonferroni)
    dunn_df = plot_df[plot_df['refined_category'].isin(CATEGORY_ORDER)][
        ['LOEUF', 'refined_category']].dropna()
    dunn_result = sp.posthoc_dunn(
        dunn_df, val_col='LOEUF', group_col='refined_category',
        p_adjust='bonferroni')

    sig_pairs = []
    for i_idx in range(len(CATEGORY_ORDER)):
        for j_idx in range(i_idx + 1, len(CATEGORY_ORDER)):
            c1, c2 = CATEGORY_ORDER[i_idx], CATEGORY_ORDER[j_idx]
            p = dunn_result.loc[c1, c2]
            if p < 0.05:
                sig_pairs.append((i_idx, j_idx, p))
    sig_pairs.sort(key=lambda x: x[1] - x[0])  # narrower first

    print(f"\n{title}")
    print(f"  Kruskal-Wallis p = {kw_p:.2e}")
    print(f"  Dunn's posthoc (Bonferroni):")
    for i_idx, j_idx, p in sig_pairs:
        print(f"    {CATEGORY_ORDER[i_idx]} vs {CATEGORY_ORDER[j_idx]}: "
              f"p = {p:.4f}")

    # Draw brackets
    y_max = max(max(d) for d in bp_data if len(d) > 0)
    bracket_y_start = y_max + 0.12
    bracket_spacing = 0.18
    for k, (i_idx, j_idx, p) in enumerate(sig_pairs):
        y_bracket = bracket_y_start + k * bracket_spacing
        p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
        draw_bracket(ax, i_idx, j_idx, y_bracket, p_text)

    top_y = bracket_y_start + max(len(sig_pairs), 1) * bracket_spacing + 0.05
    ax.text(len(CATEGORY_ORDER) / 2 - 0.5, top_y,
            f'Kruskal-Wallis p = {kw_p:.2e}',
            fontsize=16, ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.9))

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [CATEGORY_SHORT_LABELS[c] for c in CATEGORY_ORDER],
        fontsize=16, rotation=35, ha='right')
    ax.set_ylabel('LOEUF (higher = less constrained)', fontsize=20)
    ax.set_title(title, fontsize=19, fontweight='bold', loc='left', pad=14)
    ax.tick_params(labelsize=18)
    ax.set_ylim(YLIM)


# ══════════════════════════════════════════════════════════════
# COMPOSITE FIGURE
# ══════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 3, figsize=(26, 9))
fig.subplots_adjust(wspace=0.32, left=0.05, right=0.98, top=0.88, bottom=0.15)

# ── Panel A: Betweenness centrality vs. LOEUF ──
aim1_labels = (PIGMENTATION_ENZYMES | MELANOCYTE_REGULATORY
               | {'NFKB1', 'AKT1', 'KIT'})
plot_scatter(axes[0], 'BetweennessCentrality', 'LOEUF', df,
             'Betweenness centrality\n(within-pathway integration)',
             'LOEUF (higher = less constrained)',
             'Network position vs. evolutionary constraint',
             aim1_labels, legend_loc='upper right')
axes[0].text(-0.08, 1.10, 'A', transform=axes[0].transAxes,
             fontsize=28, fontweight='bold', va='top')

# ── Panel B: KEGG pathway count vs. LOEUF ──
# Exclude genes not in KEGG
df_kegg = df[df['n_kegg_pathways'] > 0].copy()
pigment_in_kegg = PIGMENTATION_ENZYMES - {'OCA2', 'PMEL', 'MLANA'}
regulatory_in_kegg = MELANOCYTE_REGULATORY & set(df_kegg['gene'])
aim2_labels = (pigment_in_kegg | regulatory_in_kegg
               | {'NFKB1', 'AKT1', 'MAPK1'})
plot_scatter(axes[1], 'n_kegg_pathways', 'LOEUF', df_kegg,
             'Number of KEGG pathways\n(cross-system connectivity)',
             'LOEUF (higher = less constrained)',
             'Pathway involvement vs. evolutionary constraint',
             aim2_labels, legend_loc='upper right')
axes[1].text(-0.08, 1.10, 'B', transform=axes[1].transAxes,
             fontsize=28, fontweight='bold', va='top')

# ── Panel C: Boxplot by refined category ──
plot_boxplot_with_pairwise(axes[2], df, 'Constraint by functional category')
axes[2].text(-0.08, 1.10, 'C', transform=axes[2].transAxes,
             fontsize=28, fontweight='bold', va='top')

# Save
fig.savefig(os.path.join(OUTPUT_DIR, 'figure_preliminary_results_composite.png'),
            dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(OUTPUT_DIR, 'figure_preliminary_results_composite.pdf'),
            dpi=300, bbox_inches='tight', facecolor='white')
print("\nComposite figure saved!")
plt.close(fig)
