"""
phase2_pbs_trees.py

Draws PBS population branch-length trees for:
  - Top 3 genes from each of the 4 PBS scans
  - MC1R (case study gene)

Each tree is an unrooted 3-leaf star (additive tree) where the branch
lengths are the individual population PBS values, derived from the
pairwise FST branch-length transforms:

    T_XY = -ln(1 - FST_XY)   (clipped to [0, 0.9999])

    PBS_X = (T_XY + T_XZ - T_YZ) / 2   ← target branch
    PBS_Y = (T_XY + T_YZ - T_XZ) / 2   ← outgroup branch
    PBS_Z = (T_XZ + T_YZ - T_XY) / 2   ← distant outgroup branch

All branch lengths are floored at 0. Branch lengths scaled so the
longest arm in each tree fills ~80% of the panel radius.

Layout: 4 scans × 4 columns (top-1, top-2, top-3, MC1R)
Output: output/figure_phase2_pbs_trees.png/.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PBS_CSV     = os.path.join(PROJECT_DIR, 'output', 'pbs_per_gene.csv')
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data',   'network_constraint_gtex.csv')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared constants ───────────────────────────────────────────────────────
CATEGORY_COLORS = {
    'Pigment-specific':          '#D94040',
    'Developmental/NC':          '#E8907E',
    'Generic signaling':         '#F5C242',
    'Cytokines/growth factors':  '#4878CF',
    'Apoptosis/cell death':      '#6BAD6B',
    'Other':                     '#B0B0B0',
}

# Colors for the 3 arm types in the tree
ARM_COLORS = {
    'target':  '#2C3E50',   # dark slate — target population
    'out':     '#7F8C8D',   # gray — outgroup
    'distant': '#BDC3C7',   # light gray — distant outgroup
}

# ── Load data ──────────────────────────────────────────────────────────────
pbs_df    = pd.read_csv(PBS_CSV)
master_df = pd.read_csv(MASTER_CSV)
master_df['gene'] = master_df['gene'].str.upper()
pbs_df['gene']    = pbs_df['gene'].str.upper()
pbs_df = pbs_df.merge(master_df[['gene', 'functional_category']], on='gene', how='left')

# Floor all PBS scan columns at 0
for col in ['pbs1_african', 'pbs2_african', 'pbs3_melanesian', 'pbs4_melanesian']:
    pbs_df[col] = pbs_df[col].clip(lower=0)

print(f"Loaded {len(pbs_df)} genes")


# ── FST → branch-length transform ─────────────────────────────────────────
def fst_to_t(fst):
    if pd.isna(fst):
        return np.nan
    return float(-np.log(1.0 - max(0.0, min(float(fst), 0.9999))))


# ── Compute all 3 branch lengths for a given gene + scan ──────────────────
# Scan definitions: (scan_col, target_label, outgroup_label, distant_label,
#                    fst_target_out, fst_target_distant, fst_out_distant)
SCAN_DEFS = [
    ('pbs1_african',
     'African', 'S. Asian', 'Melanesian',
     'fst_african_southasian', 'fst_african_melanesian', 'fst_melanesian_southasian'),

    ('pbs2_african',
     'African', 'European', 'Melanesian',
     'fst_african_european', 'fst_african_melanesian', 'fst_european_melanesian'),

    ('pbs3_melanesian',
     'Melanesian', 'S. Asian', 'African',
     'fst_melanesian_southasian', 'fst_african_melanesian', 'fst_african_southasian'),

    ('pbs4_melanesian',
     'Melanesian', 'European', 'African',
     'fst_european_melanesian', 'fst_african_melanesian', 'fst_african_european'),
]

SCAN_TITLES = {
    'pbs1_african':    'PBS-1: African\n(S. Asian out, Melanesian distant)',
    'pbs2_african':    'PBS-2: African\n(European out, Melanesian distant)',
    'pbs3_melanesian': 'PBS-3: Melanesian\n(S. Asian out, African distant)',
    'pbs4_melanesian': 'PBS-4: Melanesian\n(European out, African distant)',
}


def get_branches(row, fst_to_col, fst_td_col, fst_od_col):
    """Return (pbs_target, pbs_out, pbs_distant), all floored at 0."""
    t_to = fst_to_t(row.get(fst_to_col, np.nan))
    t_td = fst_to_t(row.get(fst_td_col, np.nan))
    t_od = fst_to_t(row.get(fst_od_col, np.nan))
    if any(np.isnan(v) for v in [t_to, t_td, t_od]):
        return np.nan, np.nan, np.nan
    pbs_target  = max(0.0, (t_to + t_td - t_od) / 2)
    pbs_out     = max(0.0, (t_to + t_od - t_td) / 2)
    pbs_distant = max(0.0, (t_td + t_od - t_to) / 2)
    return pbs_target, pbs_out, pbs_distant


# ── Tree drawing function ──────────────────────────────────────────────────
def draw_pbs_tree(ax, pbs_target, pbs_out, pbs_distant,
                  label_target, label_out, label_distant,
                  gene, scan_col, cat_color,
                  scan_pbs_value):
    """
    Draw a 3-arm star PBS tree in ax.
    Arms radiate from the origin at fixed angles:
      - Target:  top        (90°)
      - Outgroup: lower-left (210°)
      - Distant:  lower-right (330°)
    Branch lengths are proportional to PBS values, scaled so the longest
    arm reaches ~0.8 of the unit radius.
    """
    ax.set_aspect('equal')
    ax.axis('off')

    # Scale all arms by the same factor (relative to largest)
    max_len = max(pbs_target, pbs_out, pbs_distant, 1e-9)
    scale   = 0.80 / max_len

    angles = {
        'target':  90,
        'out':     210,
        'distant': 330,
    }
    lengths = {
        'target':  pbs_target  * scale,
        'out':     pbs_out     * scale,
        'distant': pbs_distant * scale,
    }
    labels = {
        'target':  label_target,
        'out':     label_out,
        'distant': label_distant,
    }
    pbs_vals = {
        'target':  pbs_target,
        'out':     pbs_out,
        'distant': pbs_distant,
    }

    for arm in ('target', 'out', 'distant'):
        angle_rad = np.radians(angles[arm])
        x_end = lengths[arm] * np.cos(angle_rad)
        y_end = lengths[arm] * np.sin(angle_rad)

        lw     = 3.5 if arm == 'target' else 1.8
        color  = cat_color if arm == 'target' else ARM_COLORS[arm]
        zorder = 4 if arm == 'target' else 3

        ax.annotate('', xy=(x_end, y_end), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='-',
                                   color=color, lw=lw,
                                   connectionstyle='arc3,rad=0'),
                    zorder=zorder)

        # Tip label — position slightly beyond the arm tip
        tip_scale = 1.18
        x_lbl = tip_scale * lengths[arm] * np.cos(angle_rad)
        y_lbl = tip_scale * lengths[arm] * np.sin(angle_rad)

        # Alignment by quadrant
        ha = 'center'
        va = 'center'
        if angles[arm] == 90:
            va = 'bottom'
        elif angles[arm] == 210:
            ha = 'right'
            va = 'top'
        elif angles[arm] == 330:
            ha = 'left'
            va = 'top'

        pop_label = labels[arm]
        pbs_str   = f'{pbs_vals[arm]:.4f}'
        weight    = 'bold' if arm == 'target' else 'normal'
        ax.text(x_lbl, y_lbl,
                f'{pop_label}\n{pbs_str}',
                ha=ha, va=va, fontsize=8.5, fontweight=weight,
                color=color if arm == 'target' else '#333333')

    # Central node dot
    ax.scatter([0], [0], c='black', s=25, zorder=5)

    # Gene name + PBS value as title
    ax.set_title(f'{gene}\nPBS = {scan_pbs_value:.4f}',
                 fontsize=9.5, fontweight='bold', style='italic',
                 color=cat_color, pad=4)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)


# ══════════════════════════════════════════════════════════════════════════
# Build figure: 4 rows (scans) × 4 columns (top-1, top-2, top-3, MC1R)
# ══════════════════════════════════════════════════════════════════════════
N_COLS = 4
N_ROWS = len(SCAN_DEFS)

fig, axes = plt.subplots(N_ROWS, N_COLS,
                         figsize=(N_COLS * 4.5, N_ROWS * 4.5))
fig.subplots_adjust(wspace=0.05, hspace=0.45,
                    left=0.04, right=0.98, top=0.93, bottom=0.04)

ROW_LETTERS = ['A', 'B', 'C', 'D']

for row_idx, (scan_col, lbl_target, lbl_out, lbl_distant,
              fst_to_col, fst_td_col, fst_od_col) in enumerate(SCAN_DEFS):

    # Top 3 genes for this scan (excluding MC1R so it always gets col 3)
    top3 = (pbs_df[pbs_df['gene'] != 'MC1R']
            .dropna(subset=[scan_col])
            .nlargest(3, scan_col)
            .reset_index(drop=True))

    mc1r_row = pbs_df[pbs_df['gene'] == 'MC1R']

    genes_this_row = list(top3['gene']) + ['MC1R']

    # Row label (scan title) on the left
    ax_left = axes[row_idx, 0]
    ax_left.text(-0.18, 0.5, SCAN_TITLES[scan_col],
                 transform=ax_left.transAxes,
                 fontsize=9, fontweight='bold', va='center', ha='right',
                 rotation=90, color='#333333')

    # Panel letter (A–D) top-left of first column
    ax_left.text(-0.18, 1.05, ROW_LETTERS[row_idx],
                 transform=ax_left.transAxes,
                 fontsize=18, fontweight='bold', va='top')

    for col_idx in range(N_COLS):
        ax = axes[row_idx, col_idx]

        if col_idx < 3:
            if col_idx >= len(top3):
                ax.axis('off')
                continue
            gene_row = top3.iloc[col_idx]
        else:
            if len(mc1r_row) == 0:
                ax.axis('off')
                ax.set_title('MC1R\n(not in dataset)', fontsize=9)
                continue
            gene_row = mc1r_row.iloc[0]

        gene      = gene_row['gene']
        cat       = gene_row.get('functional_category', 'Other')
        cat_color = CATEGORY_COLORS.get(cat, '#B0B0B0')
        scan_val  = gene_row.get(scan_col, np.nan)

        pbs_t, pbs_o, pbs_d = get_branches(
            gene_row, fst_to_col, fst_td_col, fst_od_col)

        if any(np.isnan(v) for v in [pbs_t, pbs_o, pbs_d]):
            ax.axis('off')
            ax.set_title(f'{gene}\nInsufficient data', fontsize=9)
            continue

        # Column header (only row 0)
        if row_idx == 0:
            header = (f'Top {col_idx + 1}' if col_idx < 3 else 'MC1R\n(case study)')
            ax.set_title(header, fontsize=10, fontweight='bold',
                         color='#555555', pad=28)

        draw_pbs_tree(ax,
                      pbs_target=pbs_t, pbs_out=pbs_o, pbs_distant=pbs_d,
                      label_target=lbl_target, label_out=lbl_out,
                      label_distant=lbl_distant,
                      gene=gene, scan_col=scan_col,
                      cat_color=cat_color, scan_pbs_value=scan_val)

        # MC1R column: shade background lightly
        if col_idx == 3:
            ax.set_facecolor('#FDFDE8')

# ── Shared legend: functional category colors ──────────────────────────────
cat_order = ['Pigment-specific', 'Developmental/NC', 'Generic signaling',
             'Cytokines/growth factors', 'Apoptosis/cell death', 'Other']
legend_patches = [
    mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85)
    for c in cat_order
]
fig.legend(handles=legend_patches, fontsize=8.5, ncol=6,
           loc='lower center', bbox_to_anchor=(0.5, 0.005),
           framealpha=0.9, edgecolor='gray',
           title='Target branch color = functional category',
           title_fontsize=8.5)

fig.suptitle('PBS population branch trees — top genes per scan + MC1R\n'
             'Branch lengths = PBS values  |  Bold arm = target population',
             fontsize=12, fontweight='bold', y=0.975)

for ext in ('png', 'pdf'):
    path = os.path.join(OUT_DIR, f'figure_phase2_pbs_trees.{ext}')
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {path}")
plt.close(fig)

# ── Print branch lengths table ─────────────────────────────────────────────
print("\n=== PBS branch lengths for plotted genes ===")
header = f"{'Gene':<12} {'Scan':<20} {'PBS_target':>10} {'PBS_out':>10} {'PBS_distant':>11}"
print(header)
print('-' * len(header))

for scan_col, lbl_t, lbl_o, lbl_d, fst_to_col, fst_td_col, fst_od_col in SCAN_DEFS:
    top3 = (pbs_df[pbs_df['gene'] != 'MC1R']
            .dropna(subset=[scan_col])
            .nlargest(3, scan_col))
    mc1r = pbs_df[pbs_df['gene'] == 'MC1R']
    for _, gene_row in pd.concat([top3, mc1r]).iterrows():
        gene = gene_row['gene']
        t, o, d = get_branches(gene_row, fst_to_col, fst_td_col, fst_od_col)
        scan_label = f'{scan_col}({lbl_t})'
        print(f"{gene:<12} {scan_label:<20} {t:>10.4f} {o:>10.4f} {d:>11.4f}")
    print()

print("Done!")
