"""
phase2_tau_variants.py

Tau-based variants of the phase2 network × constraint × selection figures.

Three figures, each two-panel (PBS-1 African / PBS-3 Melanesian):

  1. figure_phase2_tau_loeuf.png
       x = tissue specificity (τ, Yanai 2005); y = LOEUF
       size = PBS; color = functional category
       Tests: do tissue-specific genes show less LoF constraint?

  2. figure_phase2_tau_centrality.png
       x = tissue specificity (τ); y = within-pathway betweenness centrality
       size = PBS; color = functional category
       Tests: do tissue-specific genes occupy peripheral network positions?

  3. figure_phase2_kegg_loeuf_tau_color.png
       x = cross-system connectivity (KEGG pathway count, sqrt-scaled)
       y = LOEUF; size = PBS; fill COLOR = τ (continuous, RdYlBu_r)
       marker SHAPE = functional category
       Encodes tau as a 3rd visual axis, replacing categorical color.

τ = 0: uniformly expressed across all tissues (housekeeping)
τ = 1: expressed in only one tissue (tissue-specific)

Data sources:
  data/network_constraint_categorized.csv  — LOEUF, betweenness, tau
  data/kegg_pathway_counts.csv             — cross-system connectivity
  output/pbs_per_gene.csv                  — PBS-1, PBS-3 per gene
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy import stats
from adjustText import adjust_text

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_categorized.csv')
KEGG_CSV    = os.path.join(PROJECT_DIR, 'data',   'kegg_pathway_counts.csv')
_PBS_CANDIDATES = [
    os.path.join(PROJECT_DIR, 'data',   'pbs_per_gene.csv'),
    os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv'),
]
PBS_CSV = next((p for p in _PBS_CANDIDATES if os.path.exists(p)), _PBS_CANDIDATES[0])

_PI_CANDIDATES = [
    os.path.join(PROJECT_DIR, 'data',   'pi_per_gene.csv'),
    os.path.join(PROJECT_DIR, 'output', 'pi_per_gene.csv'),
]
PI_CSV = next((p for p in _PI_CANDIDATES if os.path.exists(p)), _PI_CANDIDATES[0])
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

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

# Marker shapes for the tau-as-color figure (category → shape)
MARKER_SHAPES = {
    'Pigment-specific':          'o',   # circle
    'Developmental/NC':          '^',   # triangle up
    'Generic signaling':         's',   # square
    'Cytokines/growth factors':  'D',   # diamond
    'Apoptosis/cell death':      'P',   # filled plus
    'Other':                     'X',   # filled X
}

TAU_CMAP = 'RdYlBu_r'   # blue=τ≈0 (housekeeping), red=τ≈1 (tissue-specific)

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading data...")
master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()

kegg = pd.read_csv(KEGG_CSV)[['gene', 'kegg_pathway_count']]
kegg['gene'] = kegg['gene'].str.upper()

HAS_PBS = os.path.exists(PBS_CSV)
if HAS_PBS:
    pbs = pd.read_csv(PBS_CSV)[['gene', 'pbs1_african', 'pbs3_melanesian']]
    pbs['gene'] = pbs['gene'].str.upper()
    pbs['pbs1_african']    = pbs['pbs1_african'].clip(lower=0)
    pbs['pbs3_melanesian'] = pbs['pbs3_melanesian'].clip(lower=0)
    print("  PBS data loaded")
else:
    print(f"  WARNING: {PBS_CSV} not found — PBS data lives on the cluster.")
    print("  Running with uniform node sizes (NO PBS gradient).")
    print("  scp from greatlakes:.../output/pbs_per_gene.csv to enable size encoding.")
    genes = master['gene'].str.upper().tolist()
    pbs = pd.DataFrame({
        'gene': genes,
        'pbs1_african':    [np.nan] * len(genes),
        'pbs3_melanesian': [np.nan] * len(genes),
    })

if os.path.exists(PI_CSV):
    pi = pd.read_csv(PI_CSV)[['gene', 'pi_african', 'pi_melanesian']]
    pi['gene'] = pi['gene'].str.upper()
    HAS_PI = True
    print("  π data loaded")
else:
    pi = pd.DataFrame({'gene': master['gene'].str.upper(),
                       'pi_african': np.nan, 'pi_melanesian': np.nan})
    HAS_PI = False
    print("  π data not found — Figure 5 will be skipped")

df = (master[['gene', 'functional_category', 'LOEUF', 'betweenness_centrality', 'tau']]
      .merge(kegg, on='gene', how='left')
      .merge(pbs,  on='gene', how='left')
      .merge(pi,   on='gene', how='left'))

print(f"  {len(df)} total genes; {df['tau'].notna().sum()} with tau")

# ── Node size scaling (PBS → marker area) ─────────────────────────────────
S_MIN, S_MAX = 20, 1400
S_FLAT = 180   # used when PBS data is unavailable

if HAS_PBS:
    pbs_global_max = max(df['pbs1_african'].max(), df['pbs3_melanesian'].max())
    def pbs_to_size(pbs_col):
        vals = df[pbs_col].fillna(0)
        return S_MIN + (S_MAX - S_MIN) * (vals / pbs_global_max) ** 2.0
    df['s_afr'] = pbs_to_size('pbs1_african')
    df['s_mel'] = pbs_to_size('pbs3_melanesian')
else:
    pbs_global_max = 0.0
    df['s_afr'] = S_FLAT
    df['s_mel'] = S_FLAT

ALWAYS_LABEL = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R', 'OCA2', 'MLANA',
    'MITF', 'SOX10', 'PAX3', 'TFAP2A', 'KIT', 'KITLG', 'EDNRB',
    'NFKB1', 'AKT1', 'MAPK1',
}

PANEL_SPECS = [
    ('pbs1_african',    's_afr', 'A',
     'PBS-1: African population-specific selection\n(S. Asian outgroup, Melanesian distant)'),
    ('pbs3_melanesian', 's_mel', 'B',
     'PBS-3: Melanesian population-specific selection\n(S. Asian outgroup, African distant)'),
]


def top_pbs_genes(pbs_col, n=8):
    return set(df.nlargest(n, pbs_col)['gene'])


def add_size_legend(ax, loc='lower right'):
    """Add the PBS size legend (or a 'data unavailable' note) directly to ax."""
    if not HAS_PBS:
        ax.text(0.97, 0.03,
                'PBS data unavailable\n(node sizes uniform)',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=8.5, style='italic', color='#666666',
                bbox=dict(boxstyle='round,pad=0.3',
                          fc='#fff8e0', ec='#cca050', alpha=0.95))
        return None
    pbs_vals = [0.0, 0.2, 0.4, 0.6, 0.8]
    handles = []
    for pv in pbs_vals:
        if pv > pbs_global_max + 0.05:
            continue
        s = S_MIN + (S_MAX - S_MIN) * (pv / pbs_global_max) ** 2.0
        h = ax.scatter([], [], s=s, c='#888888', alpha=0.75,
                       edgecolors='white', linewidths=0.5,
                       label=f'PBS = {pv:.1f}')
        handles.append(h)
    leg = ax.legend(handles=handles, title='Node size = PBS',
                    title_fontsize=9, fontsize=8.5, loc=loc,
                    framealpha=0.92, edgecolor='gray',
                    borderpad=0.8, handletextpad=1.0, labelspacing=0.8)
    ax.add_artist(leg)   # preserve when other legends are added
    return leg


def add_cat_legend(fig):
    patches = [
        mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85,
                       edgecolor='white', linewidth=0.5)
        for c in CATEGORY_ORDER
    ]
    fig.legend(handles=patches, fontsize=9.5, ncol=6,
               loc='lower center', bbox_to_anchor=(0.5, -0.02),
               framealpha=0.92, edgecolor='gray',
               title='Functional category (node color)',
               title_fontsize=9.5)


def add_labels(ax, plot_df, label_genes, pbs_col):
    texts = []
    for _, row in plot_df.iterrows():
        if row['gene'] not in label_genes:
            continue
        is_top = row['gene'] in top_pbs_genes(pbs_col, n=8)
        texts.append(ax.text(
            row['x_plot'], row['y_plot'], row['gene'],
            fontsize=8.5 if is_top else 7.5,
            fontweight='bold' if is_top else 'normal',
            style='italic', color='#111111', zorder=5))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#aaaaaa', lw=0.5),
                force_points=(1.8, 1.8), force_text=(1.8, 1.8),
                expand_points=(2.5, 2.5), expand_text=(2.0, 2.0),
                iter_lim=1000)


# ════════════════════════════════════════════════════════════════════════════
# Figure 1: tau × LOEUF
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 1: tau × LOEUF...")
df1 = df.dropna(subset=['tau', 'LOEUF']).copy()
df1['x_plot'] = df1['tau']
df1['y_plot'] = df1['LOEUF']
print(f"  {len(df1)} genes with tau + LOEUF")

r_tl, p_tl = stats.spearmanr(df1['tau'], df1['LOEUF'])
print(f"  ρ(τ, LOEUF) = {r_tl:.3f}, p = {p_tl:.3e}")
for pbs_col, label in [('pbs1_african', 'PBS-1'), ('pbs3_melanesian', 'PBS-3')]:
    sub = df1.dropna(subset=[pbs_col])
    r2, p2 = stats.spearmanr(sub['tau'], sub[pbs_col])
    print(f"  ρ(τ, {label}) = {r2:.3f}, p = {p2:.3e}")

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
fig.subplots_adjust(wspace=0.30, left=0.07, right=0.97, top=0.88, bottom=0.18)

for ax, (pbs_col, size_col, letter, subtitle) in zip(axes, PANEL_SPECS):
    label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col)
    for cat in CATEGORY_ORDER:
        sub = df1[df1['functional_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub['x_plot'], sub['y_plot'],
                   s=sub[size_col], c=CATEGORY_COLORS[cat],
                   alpha=0.82, edgecolors='white', linewidths=0.5,
                   zorder=3, label=f'{cat} (n={len(sub)})')
    add_labels(ax, df1, label_genes, pbs_col)

    ax.text(0.03, 0.97,
            f'ρ(τ, LOEUF) = {r_tl:.3f},  p = {p_tl:.2e}\nn = {len(df1)}',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.9))

    ax.set_xlabel('Tissue specificity (τ)\n← uniform (housekeeping) · · · tissue-specific →',
                  fontsize=13)
    ax.set_ylabel('LOEUF  (higher = less constrained)', fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')
    add_size_legend(ax)

add_cat_legend(fig)
fig.suptitle(
    'Tissue specificity (τ) vs. LoF intolerance and population-specific selection\n'
    'Node size = PBS value  |  τ = 0: uniformly expressed  |  τ = 1: tissue-specific',
    fontsize=12, fontweight='bold', y=0.975)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_tau_loeuf.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 2: tau × betweenness centrality
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 2: tau × betweenness centrality...")
df2 = df.dropna(subset=['tau', 'betweenness_centrality']).copy()
df2['x_plot'] = df2['tau']
df2['y_plot'] = np.sqrt(df2['betweenness_centrality'].clip(lower=0))
print(f"  {len(df2)} genes with tau + betweenness")

r_tb, p_tb = stats.spearmanr(df2['tau'], df2['betweenness_centrality'])
print(f"  ρ(τ, betweenness) = {r_tb:.3f}, p = {p_tb:.3e}")

fig, axes = plt.subplots(1, 2, figsize=(20, 10))
fig.subplots_adjust(wspace=0.30, left=0.07, right=0.97, top=0.88, bottom=0.18)

for ax, (pbs_col, size_col, letter, subtitle) in zip(axes, PANEL_SPECS):
    label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col)
    for cat in CATEGORY_ORDER:
        sub = df2[df2['functional_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub['x_plot'], sub['y_plot'],
                   s=sub[size_col], c=CATEGORY_COLORS[cat],
                   alpha=0.82, edgecolors='white', linewidths=0.5,
                   zorder=3, label=f'{cat} (n={len(sub)})')
    add_labels(ax, df2, label_genes, pbs_col)

    ax.text(0.03, 0.97,
            f'ρ(τ, betweenness) = {r_tb:.3f},  p = {p_tb:.2e}\nn = {len(df2)}',
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.9))

    y_max = df2['y_plot'].max()
    y_ticks = np.linspace(0, y_max, 6)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f'{v**2:.3f}' for v in y_ticks])

    ax.set_xlabel('Tissue specificity (τ)\n← uniform (housekeeping) · · · tissue-specific →',
                  fontsize=13)
    ax.set_ylabel('Within-pathway centrality\n(betweenness centrality)', fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')
    add_size_legend(ax)

add_cat_legend(fig)
fig.suptitle(
    'Tissue specificity (τ) vs. within-pathway centrality and population-specific selection\n'
    'Node size = PBS value  |  τ = 0: uniformly expressed  |  τ = 1: tissue-specific',
    fontsize=12, fontweight='bold', y=0.975)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_tau_centrality.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 3: KEGG × LOEUF — tau as continuous color, shape = category
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 3: KEGG × LOEUF, tau as color, shape=category...")
df3 = df.dropna(subset=['tau', 'LOEUF', 'kegg_pathway_count']).copy()
df3['x_plot'] = np.sqrt(df3['kegg_pathway_count'].clip(lower=0))
df3['y_plot'] = df3['LOEUF']
print(f"  {len(df3)} genes with tau + LOEUF + KEGG")

r_kl, p_kl = stats.spearmanr(df3['kegg_pathway_count'], df3['LOEUF'])
print(f"  ρ(KEGG, LOEUF) = {r_kl:.3f}, p = {p_kl:.3e}")
print(f"  ρ(τ, LOEUF)   = {r_tl:.3f}, p = {p_tl:.3e}")

fig, axes = plt.subplots(1, 2, figsize=(21, 10))
fig.subplots_adjust(wspace=0.32, left=0.07, right=0.91, top=0.88, bottom=0.22)

sc_ref = None
for ax, (pbs_col, size_col, letter, subtitle) in zip(axes, PANEL_SPECS):
    label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col)

    sc = ax.scatter(df3['x_plot'], df3['y_plot'],
                    s=df3[size_col],
                    c=df3['tau'], cmap=TAU_CMAP, vmin=0, vmax=1,
                    alpha=0.87, edgecolors='#444444', linewidths=0.5,
                    zorder=3)
    sc_ref = sc

    add_labels(ax, df3, label_genes, pbs_col)

    x_ticks = np.linspace(0, df3['x_plot'].max(), 6)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{v**2:.0f}' for v in x_ticks])

    ax.set_xlabel('Cross-system connectivity\n(number of KEGG pathways)', fontsize=13)
    ax.set_ylabel('LOEUF  (higher = less constrained)', fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')
    add_size_legend(ax)

# Colorbar for tau (right of panels)
cbar_ax = fig.add_axes([0.925, 0.22, 0.012, 0.56])
cb = fig.colorbar(sc_ref, cax=cbar_ax)
cb.set_label('Tissue specificity (τ)\n← uniform · · · specific →', fontsize=10, labelpad=6)
cb.ax.tick_params(labelsize=9)
cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

fig.suptitle(
    'Cross-system connectivity vs. LoF intolerance — τ as 3rd visual axis (node color)\n'
    'Node size = PBS  |  Node color = tissue specificity (τ)',
    fontsize=12, fontweight='bold', y=0.975)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_kegg_loeuf_tau_color.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 4: τ × LOEUF, color = KEGG pathway count, size = PBS
# ════════════════════════════════════════════════════════════════════════════
print("\nFigure 4: τ × LOEUF with KEGG as color, size = PBS...")
df4 = df.dropna(subset=['tau', 'LOEUF', 'kegg_pathway_count']).copy()
df4['x_plot'] = df4['tau']
df4['y_plot'] = df4['LOEUF']
print(f"  {len(df4)} genes with tau + LOEUF + KEGG")

# sqrt-scale KEGG for color (skewed distribution)
kegg_max_sqrt = float(np.sqrt(df4['kegg_pathway_count'].max()))
df4['kegg_sqrt'] = np.sqrt(df4['kegg_pathway_count'].clip(lower=0))

KEGG_CMAP = 'viridis'

fig, axes = plt.subplots(1, 2, figsize=(21, 10))
fig.subplots_adjust(wspace=0.32, left=0.07, right=0.91, top=0.88, bottom=0.18)

sc_ref = None
for ax, (pbs_col, size_col, letter, subtitle) in zip(axes, PANEL_SPECS):
    label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col)

    sc = ax.scatter(df4['x_plot'], df4['y_plot'],
                    s=df4[size_col],
                    c=df4['kegg_sqrt'], cmap=KEGG_CMAP,
                    vmin=0, vmax=kegg_max_sqrt,
                    alpha=0.88, edgecolors='#333333', linewidths=0.5,
                    zorder=3)
    sc_ref = sc

    add_labels(ax, df4, label_genes, pbs_col)

    ax.set_xlabel('Tissue specificity (τ)\n← uniform (housekeeping) · · · tissue-specific →',
                  fontsize=13)
    ax.set_ylabel('LOEUF  (higher = less constrained)', fontsize=13)
    ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
    ax.tick_params(labelsize=11)
    ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
            fontsize=22, fontweight='bold', va='top')
    add_size_legend(ax)

# Colorbar for KEGG (right of panels) — tick labels back-transformed to count
cbar_ax = fig.add_axes([0.925, 0.22, 0.012, 0.56])
cb = fig.colorbar(sc_ref, cax=cbar_ax)
cb.set_label('Cross-system connectivity\n(KEGG pathway count)', fontsize=10, labelpad=6)
cb.ax.tick_params(labelsize=9)
tick_counts = [0, 5, 19, 43, 77, 120]
tick_sqrt = [np.sqrt(c) for c in tick_counts if np.sqrt(c) <= kegg_max_sqrt + 0.1]
cb.set_ticks(tick_sqrt)
cb.set_ticklabels([f'{int(t**2)}' for t in tick_sqrt])

fig.suptitle(
    'Tissue specificity (τ) vs. LoF intolerance — KEGG connectivity as 3rd visual axis (node color)\n'
    'Node size = PBS  |  Node color = cross-system connectivity (KEGG pathway count)',
    fontsize=12, fontweight='bold', y=0.975)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_tau_loeuf_kegg_color.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved → {path}")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Figure 5: τ × LOEUF, color = π (nucleotide diversity), size = PBS
# Panel A uses π_african (matched to PBS-1 African scan)
# Panel B uses π_melanesian (matched to PBS-3 Melanesian scan)
# ════════════════════════════════════════════════════════════════════════════
if HAS_PI:
    print("\nFigure 5: τ × LOEUF with π as color, size = PBS...")
    df5 = df.dropna(subset=['tau', 'LOEUF']).copy()
    df5['x_plot'] = df5['tau']
    df5['y_plot'] = df5['LOEUF']
    print(f"  {len(df5)} genes; π_afr available for "
          f"{df5['pi_african'].notna().sum()}, π_mel for {df5['pi_melanesian'].notna().sum()}")

    PI_CMAP = 'plasma'

    PI_PANELS = [
        ('pbs1_african',    's_afr', 'pi_african',    'A',
         'PBS-1: African selection — color = π_african\n(S. Asian outgroup, Melanesian distant)'),
        ('pbs3_melanesian', 's_mel', 'pi_melanesian', 'B',
         'PBS-3: Melanesian selection — color = π_melanesian\n(S. Asian outgroup, African distant)'),
    ]

    # Shared color scale across panels (95th-percentile clip for outliers like MC1R)
    pi_all = pd.concat([df5['pi_african'].dropna(),
                        df5['pi_melanesian'].dropna()])
    pi_vmin = float(pi_all.quantile(0.05))
    pi_vmax = float(pi_all.quantile(0.95))

    fig, axes = plt.subplots(1, 2, figsize=(21, 10))
    fig.subplots_adjust(wspace=0.32, left=0.07, right=0.91, top=0.88, bottom=0.18)

    sc_ref = None
    for ax, (pbs_col, size_col, pi_col, letter, subtitle) in zip(axes, PI_PANELS):
        label_genes = ALWAYS_LABEL | top_pbs_genes(pbs_col)
        sub = df5.dropna(subset=[pi_col]).copy()

        sc = ax.scatter(sub['x_plot'], sub['y_plot'],
                        s=sub[size_col],
                        c=sub[pi_col], cmap=PI_CMAP,
                        vmin=pi_vmin, vmax=pi_vmax,
                        alpha=0.88, edgecolors='#333333', linewidths=0.5,
                        zorder=3)
        sc_ref = sc

        add_labels(ax, sub, label_genes, pbs_col)

        ax.set_xlabel('Tissue specificity (τ)\n← uniform (housekeeping) · · · tissue-specific →',
                      fontsize=13)
        ax.set_ylabel('LOEUF  (higher = less constrained)', fontsize=13)
        ax.set_title(subtitle, fontsize=12, fontweight='bold', loc='left', pad=10)
        ax.tick_params(labelsize=11)
        ax.text(-0.08, 1.06, letter, transform=ax.transAxes,
                fontsize=22, fontweight='bold', va='top')
        add_size_legend(ax)

    cbar_ax = fig.add_axes([0.925, 0.22, 0.012, 0.56])
    cb = fig.colorbar(sc_ref, cax=cbar_ax, extend='both')
    cb.set_label('Nucleotide diversity (π)\n(matched population)',
                 fontsize=10, labelpad=6)
    cb.ax.tick_params(labelsize=9)

    fig.suptitle(
        'Tissue specificity (τ) vs. LoF intolerance — π as 3rd visual axis (node color)\n'
        'Node size = PBS  |  Node color = nucleotide diversity (π) in the matched population',
        fontsize=12, fontweight='bold', y=0.975)

    for ext in ('png', 'pdf'):
        path = os.path.join(OUT_DIR, f'figure_phase2_tau_loeuf_pi_color.{ext}')
        fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved → {path}")
    plt.close(fig)


print("\nDone!")
