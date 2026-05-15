"""
figure_mc1r_pi_bars.py

Standalone version of the MC1R π vs. dataset median bar chart
(originally panel F of figure_phase2_pi).

Output: output/figure_mc1r_pi_bars.{png,pdf}
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR    = os.path.join(PROJECT_DIR, 'data')
OUT_DIR     = os.path.join(PROJECT_DIR, 'output')

POPS       = ['african', 'melanesian', 'southasian', 'european', 'eastasian']
POP_LABELS = {'african': 'African', 'melanesian': 'Melanesian',
              'southasian': 'South Asian', 'european': 'European',
              'eastasian': 'East Asian'}
POP_COLORS = {'african': '#C0392B', 'melanesian': '#2471A3',
              'southasian': '#1E8449', 'european': '#7D3C98',
              'eastasian': '#D4A017'}

pi_df = pd.read_csv(os.path.join(DATA_DIR, 'pi_per_gene.csv'))
pi_df['gene'] = pi_df['gene'].str.upper()
mc1r = pi_df[pi_df['gene'] == 'MC1R']
assert len(mc1r) == 1, "MC1R not found or duplicated in pi_per_gene.csv"

mc1r_vals = [mc1r.iloc[0][f'pi_{p}'] for p in POPS]
med_vals  = [pi_df[f'pi_{p}'].median()  for p in POPS]

print("MC1R π by population:")
for p, mv, dv in zip(POPS, mc1r_vals, med_vals):
    print(f"  {POP_LABELS[p]:12s}  MC1R π = {mv:.5f}   "
          f"dataset median = {dv:.5f}   "
          f"ratio = {mv/dv:.2f}x")

fig, ax = plt.subplots(figsize=(8.5, 6.5))
x = np.arange(len(POPS))
w = 0.35

ax.bar(x - w/2, mc1r_vals, width=w, zorder=3,
       color=[POP_COLORS[p] for p in POPS], alpha=0.90,
       edgecolor='white', linewidth=0.6, label='MC1R')

ax.bar(x + w/2, med_vals, width=w, zorder=3,
       color=[POP_COLORS[p] for p in POPS], alpha=0.30,
       edgecolor='gray', linewidth=0.6,
       hatch='///', label='Dataset median')

ax.set_xticks(x)
ax.set_xticklabels([POP_LABELS[p] for p in POPS],
                   fontsize=11, rotation=20, ha='right')
ax.set_ylabel('π (nucleotide diversity)', fontsize=12.5)
ax.set_title('MC1R π vs. dataset median across populations',
             fontsize=14, fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=10.5, framealpha=0.92, edgecolor='#999')
ax.tick_params(labelsize=10.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()

for ext in ('png', 'pdf'):
    out = os.path.join(OUT_DIR, f'figure_mc1r_pi_bars.{ext}')
    fig.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"saved {out}")
plt.close(fig)
