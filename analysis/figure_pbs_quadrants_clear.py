"""
figure_pbs_quadrants_clear.py

Replacement for the top PBS quadrants scatter in
output/figure_phase2_pbs_quadrants.png with axes that state the biological
point — high PBS on a population's branch = candidate signal of recent
positive selection in that population — and quadrant labels in each corner.

Outputs:  output/figure_pbs_quadrants_clear.{png,pdf}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(PROJECT_DIR, 'data')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

# ── Load ───────────────────────────────────────────────────────────────────
master = pd.read_csv(os.path.join(DATA_DIR, 'network_constraint_categorized.csv'))
master['gene'] = master['gene'].str.upper()
pbs    = pd.read_csv(os.path.join(DATA_DIR, 'pbs_per_gene.csv'))
pbs['gene'] = pbs['gene'].str.upper()
pbs['pbs1_african']    = pbs['pbs1_african'].clip(lower=0)
pbs['pbs3_melanesian'] = pbs['pbs3_melanesian'].clip(lower=0)
df = master.merge(pbs[['gene', 'pbs1_african', 'pbs3_melanesian']],
                  on='gene', how='left')
df = df.dropna(subset=['pbs1_african', 'pbs3_melanesian'])
print(f"{len(df)} genes with PBS values")

# ── Quadrants ──────────────────────────────────────────────────────────────
thr_a = df['pbs1_african'].quantile(0.75)
thr_m = df['pbs3_melanesian'].quantile(0.75)
print(f"75th-percentile thresholds:  African={thr_a:.3f}  Melanesian={thr_m:.3f}")

def assign(row):
    hi_a = row['pbs1_african']    >= thr_a
    hi_m = row['pbs3_melanesian'] >= thr_m
    if hi_a and hi_m: return 'Both'
    if hi_a:           return 'AfricanOnly'
    if hi_m:           return 'MelanesianOnly'
    return 'Neither'
df['quadrant'] = df.apply(assign, axis=1)
counts = df['quadrant'].value_counts().to_dict()
for k in ('AfricanOnly', 'MelanesianOnly', 'Both', 'Neither'):
    print(f"  {k:18s} n={counts.get(k, 0)}")

COLOR = {'AfricanOnly':    '#6600cc',
         'MelanesianOnly': '#9ECC1B',
         'Both':           '#00C09A',
         'Neither':        '#CFCFCF'}

# Genes to always label
ALWAYS_LABEL = {
    'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R', 'OCA2', 'MLANA',
    'MITF', 'SOX10', 'PAX3', 'TFAP2A', 'KIT', 'KITLG', 'EDNRB',
}

# ── Figure ─────────────────────────────────────────────────────────────────
# Square axes so PBS values on the two branches are visually comparable
# and the y=x diagonal is a true 45° reference.
fig, ax = plt.subplots(figsize=(12, 12))

xmax = ymax = 1.0
ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
ax.set_aspect('equal', adjustable='box')

# Quadrant background shading (very subtle)
ax.axhspan(thr_m, ymax, xmin=0, xmax=thr_a/xmax,
           facecolor=COLOR['MelanesianOnly'], alpha=0.06, zorder=0)
ax.axhspan(thr_m, ymax, xmin=thr_a/xmax, xmax=1,
           facecolor=COLOR['Both'], alpha=0.08, zorder=0)
ax.axhspan(0, thr_m, xmin=thr_a/xmax, xmax=1,
           facecolor=COLOR['AfricanOnly'], alpha=0.06, zorder=0)
ax.axhspan(0, thr_m, xmin=0, xmax=thr_a/xmax,
           facecolor='#FFFFFF', alpha=0.0, zorder=0)

# Threshold lines and y=x
ax.plot([0, max(xmax, ymax)], [0, max(xmax, ymax)],
        color='#888', lw=0.7, linestyle='--', zorder=1)
ax.axhline(thr_m, color='#888', lw=0.8, linestyle=':', zorder=1)
ax.axvline(thr_a, color='#888', lw=0.8, linestyle=':', zorder=1)
ax.text(thr_a, ymax * 0.005, f'  top 25 %\n  African',
        fontsize=8.5, color='#555', ha='left', va='bottom')
ax.text(xmax * 0.005, thr_m, f'  top 25 %\n  Melanesian',
        fontsize=8.5, color='#555', ha='left', va='bottom')

# Scatter
order = ['Neither', 'AfricanOnly', 'MelanesianOnly', 'Both']
LABEL_PRETTY = {
    'AfricanOnly':    'African-specific selection',
    'MelanesianOnly': 'Melanesian-specific selection',
    'Both':           'Shared selection (both populations)',
    'Neither':        'No selection signal',
}
for q in order:
    sub = df[df['quadrant'] == q]
    ax.scatter(sub['pbs1_african'], sub['pbs3_melanesian'],
               c=COLOR[q], s=130, alpha=0.88,
               edgecolors='white', linewidths=0.7, zorder=3,
               label=f'{LABEL_PRETTY[q]}  (n = {len(sub)})')

# Label canonical genes + top extreme genes by sum-rank
df['rank_sum'] = df['pbs1_african'].rank() + df['pbs3_melanesian'].rank()
labels = ALWAYS_LABEL | set(df.nlargest(15, 'rank_sum')['gene'])
texts = []
for _, r in df.iterrows():
    if r['gene'] in labels:
        texts.append(ax.text(r['pbs1_african'], r['pbs3_melanesian'], r['gene'],
                             fontsize=9, fontweight='bold', style='italic',
                             color='#111', zorder=5))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999', lw=0.5),
            force_points=(1.6, 1.6), force_text=(1.6, 1.6),
            expand_points=(2.0, 2.0), expand_text=(1.8, 1.8))

# Axis labels — directional, plain-English
ax.set_xlabel(
    'PBS on African branch  →  stronger candidate selection signal in Africans',
    fontsize=13, fontweight='bold', labelpad=10)
ax.set_ylabel(
    'PBS on Melanesian branch  →  stronger candidate selection in Melanesians',
    fontsize=13, fontweight='bold', labelpad=10)

# Quadrant corner labels — anchored to each quadrant's outer corner
# so they don't sit on top of the densest point clusters.
corner_props = dict(fontsize=12.5, fontweight='bold', alpha=0.85,
                    zorder=2,
                    bbox=dict(boxstyle='round,pad=0.40', facecolor='white',
                              edgecolor='#999', alpha=0.92, linewidth=0.8))

ax.text(xmax * 0.98, ymax * 0.95,
        'SHARED selection',
        color='#008572', ha='right', va='top', **corner_props)
ax.text(xmax * 0.02, ymax * 0.95,
        'MELANESIAN-specific\nselection',
        color='#5B7A0E', ha='left', va='top', **corner_props)
ax.text(xmax * 0.98, thr_m * 0.92,
        'AFRICAN-specific\nselection',
        color='#5B009E', ha='right', va='top', **corner_props)

# Title + caption
ax.set_title(
    'Where do melanogenesis-network genes show '
    'population-specific recent-selection signals?',
    fontsize=14.5, fontweight='bold', loc='left', pad=14)

leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2,
                fontsize=10.5, framealpha=0.92, edgecolor='#999',
                title='Quadrant assignment', title_fontsize=11)
ax.tick_params(labelsize=11)

fig.tight_layout(rect=[0, 0.14, 1, 1])

for ext in ('png', 'pdf'):
    out = os.path.join(OUT_DIR, f'figure_pbs_quadrants_clear.{ext}')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"  saved {out}")
plt.close(fig)
