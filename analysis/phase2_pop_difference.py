"""
phase2_pop_difference.py

Figures isolating the **African vs. Melanesian** selection contrast in the
melanogenesis network, with constraint context (LOEUF, τ, KEGG connectivity).

Three figures:

  1. figure_phase2_dpbs_constraint.png
       ΔPBS (Melanesian − African) on y-axis vs. each constraint axis (3 panels).
       Reveals which constraint regime is preferentially targeted by which
       population's selection signal.

  2. figure_phase2_pbs_pbs_scatter.png
       PBS-1 (African) × PBS-3 (Melanesian) scatter, three panels colored by
       LOEUF / τ / KEGG connectivity. Diagonal y=x marks shared selection;
       off-diagonal genes are population-specific.

  3. figure_phase2_pbs_quadrants.png
       Same scatter as (2) partitioned at the 75th-percentile thresholds into
       four quadrants. Lower row: violin plots of LOEUF, τ, KEGG distribution
       per quadrant, with Mann-Whitney p between African-specific and
       Melanesian-specific quadrants.

Inputs:
  data/network_constraint_categorized.csv  — LOEUF, betweenness, tau
  data/kegg_pathway_counts.csv             — KEGG pathway count
  data/pbs_per_gene.csv                    — PBS-1 (African), PBS-3 (Melanesian)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from adjustText import adjust_text

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_counts.csv')
PBS_CSV     = os.path.join(PROJECT_DIR, 'data', 'pbs_per_gene.csv')
if not os.path.exists(PBS_CSV):
    PBS_CSV = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

ALWAYS_LABEL = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R', 'OCA2', 'MLANA',
    'MITF', 'SOX10', 'PAX3', 'TFAP2A', 'KIT', 'KITLG', 'EDNRB',
    'NFKB1', 'AKT1', 'MAPK1', 'BCL2L1', 'BBC3', 'MAPK3', 'ADCY4',
    'PRKACA', 'FASLG', 'FZD3', 'EDN1', 'JUN', 'BAD', 'RHOA',
}

COLOR_AFR = '#6600cc'   # purple — African-specific
COLOR_MEL = '#ccff00'   # yellow-green — Melanesian-specific
COLOR_BOTH = '#00daa7'  # teal — both
COLOR_NONE = '#cccccc'  # light gray — neither

# ── Load ───────────────────────────────────────────────────────────────────
print("Loading data...")
master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()

kegg = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()

pbs = pd.read_csv(PBS_CSV)[['gene', 'pbs1_african', 'pbs3_melanesian']]
pbs['gene'] = pbs['gene'].str.upper()
pbs['pbs1_african']    = pbs['pbs1_african'].clip(lower=0)
pbs['pbs3_melanesian'] = pbs['pbs3_melanesian'].clip(lower=0)

df = (master[['gene', 'functional_category', 'LOEUF', 'tau']]
      .merge(kegg, on='gene', how='left')
      .merge(pbs,  on='gene', how='left'))
df = df.dropna(subset=['pbs1_african', 'pbs3_melanesian'])
df['delta_pbs'] = df['pbs3_melanesian'] - df['pbs1_african']
df['kegg_sqrt'] = np.sqrt(df['kegg_pathway_count'].clip(lower=0))
print(f"  {len(df)} genes with PBS-1 and PBS-3")


def label_top_genes(ax, plot_df, x_col, y_col, label_genes, n_extreme=8,
                    extreme_col=None, fontsize=8):
    """Annotate ALWAYS_LABEL genes plus top-n by |extreme_col|."""
    if extreme_col:
        top = set(plot_df.reindex(plot_df[extreme_col].abs()
                                   .sort_values(ascending=False).index)
                  .head(n_extreme)['gene'])
    else:
        top = set()
    targets = label_genes | top
    texts = []
    for _, r in plot_df.iterrows():
        if r['gene'] in targets and pd.notna(r[x_col]) and pd.notna(r[y_col]):
            texts.append(ax.text(r[x_col], r[y_col], r['gene'],
                                 fontsize=fontsize, fontweight='bold',
                                 style='italic', color='#111111', zorder=5))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.5),
                force_points=(1.6, 1.6), force_text=(1.6, 1.6),
                expand_points=(2.0, 2.0), expand_text=(1.8, 1.8),
                iter_lim=600)


# ════════════════════════════════════════════════════════════════════════════
# Figure 1: ΔPBS vs. constraint axes (3 panels)
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 1: ΔPBS vs. constraint axes...")

PANELS_DPBS = [
    ('LOEUF',     'LOEUF', 'LOEUF (higher = less constrained)'),
    ('tau',       'τ',     'Tissue specificity (τ)\n← uniform · · · tissue-specific →'),
    ('kegg_sqrt', 'KEGG',  'Cross-system connectivity\n(KEGG pathway count)'),
]

fig, axes = plt.subplots(1, 3, figsize=(22, 8))
fig.subplots_adjust(wspace=0.27, left=0.06, right=0.98, top=0.86, bottom=0.16)

for ax, (col, short_label, xlabel) in zip(axes, PANELS_DPBS):
    sub = df.dropna(subset=[col, 'delta_pbs']).copy()

    # Color by sign of ΔPBS
    colors = np.where(sub['delta_pbs'] >= 0, COLOR_MEL, COLOR_AFR)
    sizes  = 30 + 800 * (sub['delta_pbs'].abs() / df['delta_pbs'].abs().max()) ** 2

    ax.axhline(0, color='#444444', lw=0.8, linestyle='--', zorder=1)
    ax.scatter(sub[col], sub['delta_pbs'],
               s=sizes, c=colors, alpha=0.78,
               edgecolors='white', linewidths=0.5, zorder=3)

    # KEGG panel: relabel sqrt-axis ticks back to counts
    if col == 'kegg_sqrt':
        x_ticks = np.linspace(0, sub[col].max(), 6)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f'{v**2:.0f}' for v in x_ticks])

    label_top_genes(ax, sub, col, 'delta_pbs', ALWAYS_LABEL,
                    n_extreme=10, extreme_col='delta_pbs')

    rho, p = stats.spearmanr(sub[col], sub['delta_pbs'])
    n_mel = (sub['delta_pbs'] > 0).sum()
    n_afr = (sub['delta_pbs'] < 0).sum()
    ax.text(0.03, 0.97,
            f'ρ({short_label}, ΔPBS) = {rho:.3f},  p = {p:.2e}\n'
            f'n = {len(sub)}  ({n_mel} Mel ↑, {n_afr} Afr ↓)',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.92))

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel('ΔPBS  (PBS-3 Melanesian − PBS-1 African)', fontsize=12)
    ax.tick_params(labelsize=11)

# Legend
patches = [
    mpatches.Patch(facecolor=COLOR_MEL, edgecolor='white',
                   label='Melanesian-dominant (ΔPBS > 0)'),
    mpatches.Patch(facecolor=COLOR_AFR, edgecolor='white',
                   label='African-dominant (ΔPBS < 0)'),
]
fig.legend(handles=patches, fontsize=10.5, ncol=2,
           loc='lower center', bbox_to_anchor=(0.5, -0.01),
           framealpha=0.92, edgecolor='gray')

fig.suptitle(
    'ΔPBS (Melanesian − African) vs. constraint axes\n'
    'Above zero: Melanesian-specific selection  |  Below zero: African-specific selection',
    fontsize=13, fontweight='bold', y=0.97)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_dpbs_constraint.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 2: PBS-1 × PBS-3 scatter, three panels colored by LOEUF / τ / KEGG
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 2: PBS-1 × PBS-3 scatter colored by constraint...")

PANELS_COLOR = [
    ('LOEUF',     'LOEUF',    'viridis_r', None, None,
     'LOEUF  (low = more constrained)'),
    ('tau',       'τ',        'RdYlBu_r',  0.0, 1.0,
     'Tissue specificity (τ)'),
    ('kegg_sqrt', 'KEGG',     'viridis',   None, None,
     'KEGG pathway count (sqrt-scaled)'),
]

fig, axes = plt.subplots(1, 3, figsize=(22, 8))
fig.subplots_adjust(wspace=0.32, left=0.05, right=0.97, top=0.85, bottom=0.14)

for ax, (col, short_label, cmap, vmin, vmax, cbar_label) in zip(axes, PANELS_COLOR):
    sub = df.dropna(subset=[col]).copy()

    if vmin is None:
        vmin = float(sub[col].quantile(0.05))
        vmax = float(sub[col].quantile(0.95))

    # Diagonal reference line
    lim_max = max(df['pbs1_african'].max(), df['pbs3_melanesian'].max()) * 1.05
    ax.plot([0, lim_max], [0, lim_max],
            color='#888888', lw=0.8, linestyle='--', zorder=1)

    sc = ax.scatter(sub['pbs1_african'], sub['pbs3_melanesian'],
                    c=sub[col], cmap=cmap, vmin=vmin, vmax=vmax,
                    s=110, alpha=0.88,
                    edgecolors='#333333', linewidths=0.5, zorder=3)

    # Label genes far from the diagonal (population-specific)
    sub['off_diag'] = sub['pbs3_melanesian'] - sub['pbs1_african']
    label_top_genes(ax, sub, 'pbs1_african', 'pbs3_melanesian',
                    ALWAYS_LABEL, n_extreme=10, extreme_col='off_diag')

    cb = plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(cbar_label, fontsize=10)
    cb.ax.tick_params(labelsize=8)
    if col == 'kegg_sqrt':
        tick_counts = [0, 5, 19, 43, 77, 120]
        tick_sqrt = [np.sqrt(c) for c in tick_counts if np.sqrt(c) <= sub[col].max() + 0.1]
        cb.set_ticks(tick_sqrt)
        cb.set_ticklabels([f'{int(t**2)}' for t in tick_sqrt])

    ax.set_xlim(0, lim_max)
    ax.set_ylim(0, lim_max)
    ax.set_xlabel('PBS-1 African (S. Asian outgroup)', fontsize=12)
    ax.set_ylabel('PBS-3 Melanesian (S. Asian outgroup)', fontsize=12)
    ax.set_title(f'Color = {short_label}', fontsize=12, fontweight='bold', loc='left', pad=8)
    ax.tick_params(labelsize=10)

    # Annotate the diagonal
    ax.text(0.97, 0.03, 'y = x\n(equal selection)',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.5, color='#666', style='italic')

fig.suptitle(
    'African vs. Melanesian PBS — population-specific selection coloured by constraint\n'
    'Off-diagonal = population-specific  |  Above y=x: Melanesian  |  Below y=x: African',
    fontsize=13, fontweight='bold', y=0.97)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_pbs_scatter.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 3: Quadrant analysis with constraint violins
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 3: Quadrant analysis with constraint violins...")

thr_a = df['pbs1_african'].quantile(0.75)
thr_m = df['pbs3_melanesian'].quantile(0.75)


def assign_quadrant(row):
    hi_a = row['pbs1_african']    >= thr_a
    hi_m = row['pbs3_melanesian'] >= thr_m
    if hi_a and hi_m: return 'Both'
    if hi_a:           return 'African'
    if hi_m:           return 'Melanesian'
    return 'Neither'


df['quadrant'] = df.apply(assign_quadrant, axis=1)
QUAD_ORDER  = ['African', 'Melanesian', 'Both', 'Neither']
QUAD_COLORS = {'African': COLOR_AFR, 'Melanesian': COLOR_MEL,
               'Both': COLOR_BOTH, 'Neither': COLOR_NONE}

print("  Quadrant counts:")
for q in QUAD_ORDER:
    print(f"    {q:11s}: n = {(df['quadrant'] == q).sum():>3d}")

# Mann-Whitney: African vs. Melanesian for each constraint
print("\n  Mann-Whitney (African-specific vs. Melanesian-specific):")
mw_results = {}
for col, label in [('LOEUF', 'LOEUF'), ('tau', 'τ'), ('kegg_pathway_count', 'KEGG')]:
    a = df.loc[df['quadrant'] == 'African',    col].dropna()
    m = df.loc[df['quadrant'] == 'Melanesian', col].dropna()
    if len(a) >= 3 and len(m) >= 3:
        u, p = stats.mannwhitneyu(a, m, alternative='two-sided')
        mw_results[col] = (u, p)
        print(f"    {label:6s}: U = {u:.1f}, p = {p:.3e}  (n_afr = {len(a)}, n_mel = {len(m)})")
    else:
        mw_results[col] = (np.nan, np.nan)

fig = plt.figure(figsize=(20, 13))
gs = fig.add_gridspec(2, 3, height_ratios=[1.4, 1.0],
                      hspace=0.32, wspace=0.30,
                      left=0.06, right=0.98, top=0.92, bottom=0.07)

# Top: PBS-1 × PBS-3 scatter (spans 3 cols)
ax_top = fig.add_subplot(gs[0, :])
lim_max = max(df['pbs1_african'].max(), df['pbs3_melanesian'].max()) * 1.05
ax_top.plot([0, lim_max], [0, lim_max], color='#888888',
            lw=0.7, linestyle='--', zorder=1)
ax_top.axhline(thr_m, color='#bbbbbb', lw=0.6, linestyle=':', zorder=1)
ax_top.axvline(thr_a, color='#bbbbbb', lw=0.6, linestyle=':', zorder=1)

for q in QUAD_ORDER:
    sub = df[df['quadrant'] == q]
    ax_top.scatter(sub['pbs1_african'], sub['pbs3_melanesian'],
                   c=QUAD_COLORS[q], s=130, alpha=0.85,
                   edgecolors='white', linewidths=0.6, zorder=3,
                   label=f'{q} (n={len(sub)})')

# Label outliers
df['rank_sum'] = df['pbs1_african'].rank() + df['pbs3_melanesian'].rank()
top_rs = set(df.nlargest(15, 'rank_sum')['gene'])
label_top_genes(ax_top, df, 'pbs1_african', 'pbs3_melanesian',
                ALWAYS_LABEL | top_rs, n_extreme=0, fontsize=9)

ax_top.set_xlim(0, lim_max)
ax_top.set_ylim(0, lim_max)
ax_top.set_xlabel('PBS-1 African (S. Asian outgroup)', fontsize=13)
ax_top.set_ylabel('PBS-3 Melanesian (S. Asian outgroup)', fontsize=13)
ax_top.set_title('African × Melanesian PBS quadrants  '
                 '(thresholds = 75th percentile of each axis)',
                 fontsize=13, fontweight='bold', loc='left', pad=10)
ax_top.legend(fontsize=10, loc='upper right', framealpha=0.92, edgecolor='gray')
ax_top.tick_params(labelsize=11)
ax_top.text(thr_a, lim_max*0.98, f' x = {thr_a:.3f}', fontsize=8, color='#888', va='top')
ax_top.text(lim_max*0.98, thr_m, f' y = {thr_m:.3f} ', fontsize=8, color='#888',
            ha='right', va='bottom')

# Bottom: 3 violins
VIOLIN_PANELS = [
    (gs[1, 0], 'LOEUF',              'LOEUF (higher = less constrained)'),
    (gs[1, 1], 'tau',                'Tissue specificity (τ)'),
    (gs[1, 2], 'kegg_pathway_count', 'KEGG pathway count'),
]

for slot, col, ylabel in VIOLIN_PANELS:
    ax = fig.add_subplot(slot)
    data = [df.loc[df['quadrant'] == q, col].dropna().values for q in QUAD_ORDER]
    parts = ax.violinplot(data, showmeans=False, showmedians=True, widths=0.8)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(QUAD_COLORS[QUAD_ORDER[i]])
        pc.set_edgecolor('#333')
        pc.set_alpha(0.78)
    for k in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        if k in parts:
            parts[k].set_color('#222')
            parts[k].set_lw(1.0)

    # Overlay individual points
    for i, vals in enumerate(data):
        x_jitter = np.random.normal(i + 1, 0.045, len(vals))
        ax.scatter(x_jitter, vals, s=14, c=QUAD_COLORS[QUAD_ORDER[i]],
                   edgecolors='#222', linewidths=0.3, alpha=0.7, zorder=4)

    ax.set_xticks(np.arange(1, len(QUAD_ORDER) + 1))
    ax.set_xticklabels([f'{q}\n(n={(df["quadrant"]==q).sum()})'
                        for q in QUAD_ORDER], fontsize=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='y', labelsize=10)

    u, p = mw_results.get(col, (np.nan, np.nan))
    if pd.notna(p):
        sig = ' ***' if p < 0.001 else ' **' if p < 0.01 else ' *' if p < 0.05 else ''
        ax.set_title(f'MW (Afr vs. Mel): p = {p:.2e}{sig}',
                     fontsize=11, loc='left', pad=6,
                     fontweight='bold' if p < 0.05 else 'normal')

fig.suptitle(
    'African vs. Melanesian PBS quadrants and their constraint profiles\n'
    'Lower row: distribution of LOEUF / τ / KEGG within each quadrant',
    fontsize=14, fontweight='bold', y=0.985)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_quadrants.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


print("\nDone!")
