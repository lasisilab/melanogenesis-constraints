"""
phase2_pi_pbs_analysis.py

Phase 2: Nucleotide diversity (π) and PBS visualisation for the Raghunath
melanogenesis network (129 genes).

Reads:
  output/pi_per_gene.csv   — per-gene π for 5 populations
  output/pbs_per_gene.csv  — per-gene PBS for 4 approved scans
  data/network_constraint_gtex.csv — functional category + LOEUF

Produces:
  output/figure_phase2_pi.png/.pdf     — 3-panel π figure
      Panel A: π by functional category, all populations (grouped boxplot)
      Panel B: π_african vs π_melanesian scatter
      Panel C: MC1R π across populations (bar chart, highlighted)

  output/figure_phase2_pbs.png/.pdf    — 3-panel PBS figure
      Panel A: PBS-1 (African, SAS outgroup) by functional category
      Panel B: PBS-3 (Melanesian, SAS outgroup) by functional category
      Panel C: Top-20 genes by PBS-1 and PBS-3 (lollipop)

Statistics printed:
  - Kruskal-Wallis + Dunn's posthoc on π_african by category
  - Kruskal-Wallis + Dunn's posthoc on π_melanesian by category
  - Spearman ρ(π_african, LOEUF) and ρ(π_melanesian, LOEUF)
  - MC1R π and PBS values vs. dataset median
  - Top genes by each PBS scan
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import scikit_posthocs as sp

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PI_CSV      = os.path.join(PROJECT_DIR, 'output', 'pi_per_gene.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_gtex.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared constants ───────────────────────────────────────────────────────
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

# Population display names and colors
POP_LABELS = {
    'african':    'African',
    'melanesian': 'Melanesian',
    'southasian': 'South Asian',
    'european':   'European',
    'eastasian':  'East Asian',
}

POP_COLORS = {
    'african':    '#C0392B',
    'melanesian': '#2471A3',
    'southasian': '#1E8449',
    'european':   '#7D3C98',
    'eastasian':  '#D4A017',
}

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading data...")
pi_df     = pd.read_csv(PI_CSV)
pbs_df    = pd.read_csv(PBS_CSV)
master_df = pd.read_csv(MASTER_CSV)
master_df['gene'] = master_df['gene'].str.upper()
pi_df['gene']     = pi_df['gene'].str.upper()
pbs_df['gene']    = pbs_df['gene'].str.upper()

# Merge with functional category and LOEUF
pi_df  = pi_df.merge(master_df[['gene', 'functional_category', 'LOEUF']],
                     on='gene', how='left')
pbs_df = pbs_df.merge(master_df[['gene', 'functional_category', 'LOEUF']],
                      on='gene', how='left')

print(f"  π: {len(pi_df)} genes")
print(f"  PBS: {len(pbs_df)} genes")

# Floor PBS at 0: negative values mean no excess differentiation signal
pbs_scan_cols = ['pbs1_african', 'pbs2_african', 'pbs3_melanesian', 'pbs4_melanesian']
for col in pbs_scan_cols:
    if col in pbs_df.columns:
        n_neg = (pbs_df[col] < 0).sum()
        if n_neg > 0:
            print(f"  Flooring {n_neg} negative {col} values to 0")
        pbs_df[col] = pbs_df[col].clip(lower=0)

POPS = ['african', 'melanesian', 'southasian', 'european', 'eastasian']

# ══════════════════════════════════════════════════════════════════════════
# Statistics
# ══════════════════════════════════════════════════════════════════════════
print("\n=== π Statistics ===")
for pop in ['african', 'melanesian']:
    col = f'pi_{pop}'
    sub = pi_df.dropna(subset=[col, 'functional_category'])
    rho, pval = stats.spearmanr(sub[col], sub['LOEUF'].fillna(np.nan))
    print(f"\nSpearman ρ(π_{pop}, LOEUF) = {rho:.3f}, p = {pval:.3e}  (n={sub['LOEUF'].notna().sum()})")

    groups = [sub.loc[sub['functional_category'] == c, col].values
              for c in CATEGORY_ORDER if c in sub['functional_category'].values]
    non_empty = [g for g in groups if len(g) > 0]
    if len(non_empty) >= 2:
        kw_stat, kw_p = stats.kruskal(*non_empty)
        print(f"Kruskal-Wallis π_{pop} by category: H = {kw_stat:.2f}, p = {kw_p:.3e}")

        dunn = sp.posthoc_dunn(sub, val_col=col,
                               group_col='functional_category',
                               p_adjust='bonferroni')
        cats = [c for c in CATEGORY_ORDER if c in sub['functional_category'].values]
        sig_found = False
        for i in range(len(cats)):
            for j in range(i + 1, len(cats)):
                p = dunn.loc[cats[i], cats[j]]
                if p < 0.05:
                    if not sig_found:
                        print(f"  Dunn's posthoc significant pairs:")
                        sig_found = True
                    print(f"    {cats[i]} vs {cats[j]}: p = {p:.4f}")
        if not sig_found:
            print("  No significant pairwise differences (Bonferroni)")

# MC1R case study
print("\n=== MC1R π across populations ===")
mc1r_pi = pi_df[pi_df['gene'] == 'MC1R']
if len(mc1r_pi) > 0:
    row = mc1r_pi.iloc[0]
    for pop in POPS:
        col = f'pi_{pop}'
        if col in row:
            med = pi_df[col].median()
            pct = (pi_df[col] < row[col]).mean() * 100
            print(f"  π_{pop}: {row[col]:.6f}  "
                  f"(dataset median={med:.6f}, percentile={pct:.0f}th)")

print("\n=== MC1R PBS scores ===")
mc1r_pbs = pbs_df[pbs_df['gene'] == 'MC1R']
if len(mc1r_pbs) > 0:
    row = mc1r_pbs.iloc[0]
    for scan, label in [('pbs1_african', 'PBS-1 African (SAS out)'),
                        ('pbs2_african', 'PBS-2 African (EUR out)'),
                        ('pbs3_melanesian', 'PBS-3 Melanesian (SAS out)'),
                        ('pbs4_melanesian', 'PBS-4 Melanesian (EUR out)')]:
        if scan in row:
            med = pbs_df[scan].median()
            pct = (pbs_df[scan] < row[scan]).mean() * 100
            print(f"  {label}: {row[scan]:.4f}  "
                  f"(dataset median={med:.4f}, percentile={pct:.0f}th)")

print("\n=== Top 10 genes by PBS scan ===")
for scan, label in [('pbs1_african',    'PBS-1 African (SAS outgroup)'),
                    ('pbs2_african',    'PBS-2 African (EUR outgroup)'),
                    ('pbs3_melanesian', 'PBS-3 Melanesian (SAS outgroup)'),
                    ('pbs4_melanesian', 'PBS-4 Melanesian (EUR outgroup)')]:
    top = pbs_df.dropna(subset=[scan]).nlargest(10, scan)[
        ['gene', 'functional_category', scan]]
    print(f"\n{label}:")
    print(top.to_string(index=False))


# ══════════════════════════════════════════════════════════════════════════
# Helper: significance bracket
# ══════════════════════════════════════════════════════════════════════════
def draw_bracket(ax, x1, x2, y, p_text, y_range, lw=1.0, fontsize=11):
    h = y_range * 0.025
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h],
            color='black', lw=lw, zorder=5)
    ax.text((x1 + x2) / 2, y + h * 0.3, p_text,
            ha='center', va='bottom', fontsize=fontsize, zorder=5)


def category_boxplot(ax, data, y_col, title, ylabel, panel_letter,
                     highlight_gene=None):
    """LOEUF-style boxplot of y_col by functional category."""
    bp_data, colors, ns = [], [], []
    for cat in CATEGORY_ORDER:
        vals = data.loc[data['functional_category'] == cat, y_col].dropna().values
        bp_data.append(vals)
        colors.append(CATEGORY_COLORS[cat])
        ns.append(len(vals))

    positions = list(range(len(CATEGORY_ORDER)))
    bp = ax.boxplot(bp_data, positions=positions, widths=0.55,
                    patch_artist=True, showfliers=False, zorder=2,
                    medianprops=dict(linewidth=0),
                    whiskerprops=dict(color='#666666', linewidth=1.1),
                    capprops=dict(color='#666666', linewidth=1.1),
                    boxprops=dict(linewidth=0.9))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.80)
        patch.set_edgecolor('#444444')

    rng = np.random.default_rng(42)
    for i, (vals, cat) in enumerate(zip(bp_data, CATEGORY_ORDER)):
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.14, 0.14, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals,
                   c='white', s=18, alpha=0.7,
                   edgecolors=CATEGORY_COLORS[cat], linewidths=0.6, zorder=3)

    # Median as white open circle
    for i, vals in enumerate(bp_data):
        if len(vals) == 0:
            continue
        ax.scatter([i], [np.median(vals)], marker='o', c='white', s=45,
                   zorder=5, edgecolors='black', linewidths=1.1)

    # Highlight a specific gene
    if highlight_gene and highlight_gene in data['gene'].values:
        row = data[data['gene'] == highlight_gene].iloc[0]
        if pd.notna(row[y_col]) and row['functional_category'] in CATEGORY_ORDER:
            cat_idx = CATEGORY_ORDER.index(row['functional_category'])
            ax.scatter([cat_idx], [row[y_col]], marker='*', c='black',
                       s=160, zorder=6, label=highlight_gene)
            ax.annotate(highlight_gene, xy=(cat_idx, row[y_col]),
                        xytext=(cat_idx + 0.35, row[y_col]),
                        fontsize=9, fontweight='bold', style='italic',
                        color='#222222',
                        arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.6))

    all_vals = [v for vals in bp_data for v in vals]
    if not all_vals:
        return
    y_min, y_max = min(all_vals), max(all_vals)
    y_range = y_max - y_min
    y_label_pos = y_min - y_range * 0.09

    for i, n in enumerate(ns):
        ax.text(i, y_label_pos, f'n={n}', ha='center', fontsize=9,
                color='#555555', fontstyle='italic', zorder=4)

    # Significance brackets from Dunn's posthoc
    sub = data.dropna(subset=[y_col, 'functional_category'])
    cats_present = [c for c in CATEGORY_ORDER
                    if c in sub['functional_category'].values and
                    len(sub[sub['functional_category'] == c]) >= 2]
    if len(cats_present) >= 2:
        try:
            dunn = sp.posthoc_dunn(sub, val_col=y_col,
                                   group_col='functional_category',
                                   p_adjust='bonferroni')
            sig_pairs = []
            for i in range(len(cats_present)):
                for j in range(i + 1, len(cats_present)):
                    c1, c2 = cats_present[i], cats_present[j]
                    if c1 in dunn.index and c2 in dunn.columns:
                        p = dunn.loc[c1, c2]
                        if p < 0.05:
                            idx1 = CATEGORY_ORDER.index(c1)
                            idx2 = CATEGORY_ORDER.index(c2)
                            sig_pairs.append((idx1, idx2, p))
            sig_pairs.sort(key=lambda x: x[1] - x[0])

            bracket_y = y_max + y_range * 0.08
            spacing   = y_range * 0.11
            for k, (i, j, p) in enumerate(sig_pairs):
                p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
                draw_bracket(ax, i, j, bracket_y + k * spacing, p_text, y_range)
            top_y = bracket_y + max(len(sig_pairs), 1) * spacing
        except Exception:
            top_y = y_max + y_range * 0.2
    else:
        top_y = y_max + y_range * 0.2

    ax.set_ylim(y_label_pos - abs(y_label_pos) * 0.1,
                top_y + y_range * 0.08)
    ax.set_xticks(positions)
    ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                       fontsize=10, rotation=30, ha='right')
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=10)
    ax.text(-0.10, 1.06, panel_letter, transform=ax.transAxes,
            fontsize=20, fontweight='bold', va='top')


# ══════════════════════════════════════════════════════════════════════════
# Tables
# ══════════════════════════════════════════════════════════════════════════
pi_cols = [f'pi_{p}' for p in POPS]

# Table 1: Pigment-specific gene π values (all populations)
print("\n=== Table 1: Pigment-specific gene π values ===")
pig_tbl = (pi_df[pi_df['functional_category'] == 'Pigment-specific']
           [['gene'] + pi_cols]
           .sort_values('pi_african', ascending=False)
           .reset_index(drop=True))
pig_tbl.columns = ['gene'] + [POP_LABELS[p] for p in POPS]
pig_tbl_print = pig_tbl.copy()
for col in [POP_LABELS[p] for p in POPS]:
    pig_tbl_print[col] = pig_tbl_print[col].map(
        lambda x: f'{x:.6f}' if pd.notna(x) else 'NA')
print(pig_tbl_print.to_string(index=False))
pig_tbl.to_csv(os.path.join(OUT_DIR, 'table_pi_pigment_specific.csv'), index=False)
print(f"  Saved → output/table_pi_pigment_specific.csv")

# Table 2: Top 10 genes by π divergence across populations
# Divergence = max(π) − min(π) across the 5 populations
print("\n=== Table 2: Top 10 genes by π divergence across populations ===")
pi_mat = pi_df[pi_cols].copy()
pi_df['pi_max']       = pi_mat.max(axis=1)
pi_df['pi_min']       = pi_mat.min(axis=1)
pi_df['pi_divergence'] = pi_df['pi_max'] - pi_df['pi_min']
pi_df['pop_max']      = pi_mat.idxmax(axis=1).str.replace('pi_', '')
pi_df['pop_min']      = pi_mat.idxmin(axis=1).str.replace('pi_', '')

div_cols = ['gene', 'functional_category', 'pi_divergence',
            'pop_max', 'pi_max', 'pop_min', 'pi_min']
top_div = (pi_df.dropna(subset=['pi_divergence'])
           .nlargest(10, 'pi_divergence')[div_cols]
           .reset_index(drop=True))
top_div_print = top_div.copy()
for col in ['pi_divergence', 'pi_max', 'pi_min']:
    top_div_print[col] = top_div_print[col].map(lambda x: f'{x:.6f}')
print(top_div_print.to_string(index=False))
top_div.to_csv(os.path.join(OUT_DIR, 'table_pi_top_divergence.csv'), index=False)
print(f"  Saved → output/table_pi_top_divergence.csv")

# ══════════════════════════════════════════════════════════════════════════
# Figure 1: π  —  2 × 3 layout
#   Row 1: African | South Asian | East Asian
#   Row 2: Melanesian | European | MC1R bar chart
# ══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(21, 13))
gs  = fig.add_gridspec(2, 3, wspace=0.35, hspace=0.55,
                       left=0.06, right=0.98, top=0.93, bottom=0.15)

panel_specs = [
    (gs[0, 0], 'pi_african',    'Nucleotide diversity (π)\nAfrican',     'A'),
    (gs[0, 1], 'pi_southasian', 'Nucleotide diversity (π)\nSouth Asian',  'B'),
    (gs[0, 2], 'pi_eastasian',  'Nucleotide diversity (π)\nEast Asian',   'C'),
    (gs[1, 0], 'pi_melanesian', 'Nucleotide diversity (π)\nMelanesian',  'D'),
    (gs[1, 1], 'pi_european',   'Nucleotide diversity (π)\nEuropean',    'E'),
]

for spec, col, title, letter in panel_specs:
    ax = fig.add_subplot(spec)
    category_boxplot(ax, pi_df, col,
                     title=title,
                     ylabel='π (nucleotide diversity)',
                     panel_letter=letter,
                     highlight_gene='MC1R')

# Panel F: MC1R π vs. dataset median across all 5 populations
ax_F = fig.add_subplot(gs[1, 2])
mc1r_row = pi_df[pi_df['gene'] == 'MC1R']
dataset_medians = {pop: pi_df[f'pi_{pop}'].median() for pop in POPS}
mc1r_vals = ([mc1r_row.iloc[0].get(f'pi_{p}', np.nan) for p in POPS]
             if len(mc1r_row) > 0 else [np.nan] * len(POPS))
med_vals  = [dataset_medians[p] for p in POPS]

x = np.arange(len(POPS))
w = 0.35
ax_F.bar(x - w/2, mc1r_vals, width=w, zorder=3,
         color=[POP_COLORS[p] for p in POPS], alpha=0.85,
         edgecolor='white', linewidth=0.5, label='MC1R')
ax_F.bar(x + w/2, med_vals, width=w, zorder=3,
         color=[POP_COLORS[p] for p in POPS], alpha=0.35,
         edgecolor='gray', linewidth=0.5,
         hatch='///', label='Dataset median')

ax_F.set_xticks(x)
ax_F.set_xticklabels([POP_LABELS[p] for p in POPS],
                     fontsize=10, rotation=20, ha='right')
ax_F.set_ylabel('π (nucleotide diversity)', fontsize=12)
ax_F.set_title('MC1R π vs. dataset median\nacross populations',
               fontsize=13, fontweight='bold', loc='left', pad=8)
ax_F.legend(fontsize=10, framealpha=0.9)
ax_F.tick_params(labelsize=10)
ax_F.text(-0.10, 1.06, 'F', transform=ax_F.transAxes,
          fontsize=20, fontweight='bold', va='top')

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pi.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {path}")
plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
# Figure 2: PBS  —  2 × 3 layout
#   Row 1: PBS-1 boxplot | PBS-2 boxplot | lollipop (all 4 scans)
#   Row 2: PBS-3 boxplot | PBS-4 boxplot | (lollipop spans both rows)
# ══════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(21, 13))
gs  = fig.add_gridspec(2, 3, wspace=0.35, hspace=0.55,
                       left=0.06, right=0.98, top=0.93, bottom=0.15)

ax_A  = fig.add_subplot(gs[0, 0])
ax_B  = fig.add_subplot(gs[0, 1])
ax_C  = fig.add_subplot(gs[1, 0])
ax_D  = fig.add_subplot(gs[1, 1])
ax_E  = fig.add_subplot(gs[:, 2])   # lollipop spans both rows

# Panel A: PBS-1 (African, S. Asian outgroup)
category_boxplot(ax_A, pbs_df, 'pbs1_african',
                 title='PBS-1: African selection\n(S. Asian outgroup, Melanesian distant)',
                 ylabel='PBS', panel_letter='A', highlight_gene='MC1R')

# Panel B: PBS-2 (African, European outgroup)
category_boxplot(ax_B, pbs_df, 'pbs2_african',
                 title='PBS-2: African selection\n(European outgroup, Melanesian distant)',
                 ylabel='PBS', panel_letter='B', highlight_gene='MC1R')

# Panel C: PBS-3 (Melanesian, S. Asian outgroup)
category_boxplot(ax_C, pbs_df, 'pbs3_melanesian',
                 title='PBS-3: Melanesian selection\n(S. Asian outgroup, African distant)',
                 ylabel='PBS', panel_letter='C', highlight_gene='MC1R')

# Panel D: PBS-4 (Melanesian, European outgroup)
category_boxplot(ax_D, pbs_df, 'pbs4_melanesian',
                 title='PBS-4: Melanesian selection\n(European outgroup, African distant)',
                 ylabel='PBS', panel_letter='D', highlight_gene='MC1R')

# Panel E: Top-15 genes lollipop — all 4 scans
top_n = 15
top_genes = set()
for col in pbs_scan_cols:
    top_genes |= set(pbs_df.dropna(subset=[col]).nlargest(top_n, col)['gene'])

plot_df = pbs_df[pbs_df['gene'].isin(top_genes)].copy()
plot_df['_sort'] = (plot_df['pbs1_african'].fillna(0) +
                    plot_df['pbs2_african'].fillna(0) +
                    plot_df['pbs3_melanesian'].fillna(0) +
                    plot_df['pbs4_melanesian'].fillna(0))
plot_df = plot_df.sort_values('_sort', ascending=True).reset_index(drop=True)

y_pos  = np.arange(len(plot_df))
colors_lolly = [CATEGORY_COLORS.get(c, '#B0B0B0')
                for c in plot_df['functional_category']]

# 4 bars per gene, offset vertically
SCAN_STYLES = [
    ('pbs1_african',    '#C0392B', 'PBS-1: African / S. Asian outgroup'),
    ('pbs2_african',    '#E8907E', 'PBS-2: African / European outgroup'),
    ('pbs3_melanesian', '#2471A3', 'PBS-3: Melanesian / S. Asian outgroup'),
    ('pbs4_melanesian', '#85C1E9', 'PBS-4: Melanesian / European outgroup'),
]
offsets = [-0.27, -0.09, 0.09, 0.27]
h = 0.16

for (col, color, label), offset in zip(SCAN_STYLES, offsets):
    vals = plot_df[col].fillna(0).values
    ax_E.barh(y_pos + offset, vals, height=h,
              color=color, alpha=0.85, label=label)
    ax_E.scatter(vals, y_pos + offset,
                 c=colors_lolly, s=22, zorder=4,
                 edgecolors='white', linewidths=0.35)

# Highlight MC1R row
if 'MC1R' in plot_df['gene'].values:
    mc1r_y = y_pos[plot_df[plot_df['gene'] == 'MC1R'].index[0]]
    ax_E.axhspan(mc1r_y - 0.4, mc1r_y + 0.4,
                 color='#FDFDE8', zorder=0, alpha=0.8)
    ax_E.annotate('MC1R', xy=(0, mc1r_y), xytext=(-0.005, mc1r_y),
                  fontsize=8, fontweight='bold', color='#555500',
                  ha='right', va='center')

ax_E.set_yticks(y_pos)
ax_E.set_yticklabels(plot_df['gene'], fontsize=8.5, style='italic')
ax_E.set_xlabel('PBS', fontsize=12)
ax_E.set_title('Top genes by PBS\n(all 4 scans)',
               fontsize=13, fontweight='bold', loc='left', pad=8)
ax_E.axvline(0, color='gray', lw=0.8, zorder=0)
ax_E.tick_params(labelsize=9)
ax_E.text(-0.10, 1.03, 'E', transform=ax_E.transAxes,
          fontsize=20, fontweight='bold', va='top')

# Legend 1: bar colors = PBS scan
bar_patches = [mpatches.Patch(facecolor=c, alpha=0.85, label=lbl)
               for _, c, lbl in SCAN_STYLES]
leg1 = ax_E.legend(handles=bar_patches, fontsize=7.5, loc='upper right',
                   framealpha=0.9, title='Bar = PBS scan', title_fontsize=7.5)
ax_E.add_artist(leg1)

# Legend 2: dot colors = functional category
cat_patches = [mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85)
               for c in CATEGORY_ORDER if c in plot_df['functional_category'].values]
ax_E.legend(handles=cat_patches, fontsize=7.5, loc='lower right',
            framealpha=0.9, title='Dot = category', title_fontsize=7.5)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {path}")
plt.close(fig)

print("\nDone!")
