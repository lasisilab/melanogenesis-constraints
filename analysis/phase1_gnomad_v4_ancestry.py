"""
phase1_gnomad_v4_ancestry.py

Phase 1.3: Per-ancestry LOEUF for the Raghunath melanogenesis network.

Per-ancestry LOEUF source:
  Han et al. 2025, Nat Commun, Supplementary Data 11
  data/gnomad_v4_ancestry_loeuf.csv  (copied from 41467_2025_57885_MOESM11_ESM.csv)
  Columns: Gene, Maximally.Diverse..n.43k., NFE..n.43k., NFE..n.440k.,
           Full.Dataset..n.460k., AFR, ASJ, EAS, SAS, NFE..n.20k.
  Values: LOEUF (upper bound of LoF o/e 90% CI).
  Gene column carries a leading apostrophe — stripped on load.
  No AMR or FIN column; those panels are omitted.

Global LOEUF baseline:
  data/gnomad_v4_constraint.tsv  (lof.oe_ci.upper)
  Used only for Comparison C (v2.1.1 vs v4 global LOEUF sanity check).

Headline comparisons:
  A. Pigmentation-gene set vs genome-wide background, per ancestry (Mann-Whitney U)
  B. NFE vs AFR LOEUF within the pigmentation set (paired Wilcoxon)
  C. v2.1.1 LOEUF vs v4 global LOEUF (Spearman, sanity check)
  D. NFE vs AFR by functional category
  E. MC1R LOEUF across all ancestries

Outputs:
  data/network_constraint_ancestry_loeuf.csv
  output/figure_phase1_ancestry_loeuf.png / .pdf

Caveat: NFE n=440,000 vs AFR n=8,701 — sample-size asymmetry means AFR LOEUF
estimates carry wider uncertainty, especially for rare-variant–poor genes.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import scikit_posthocs as sp
from adjustText import adjust_text

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR        = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV         = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
V4_TSV             = os.path.join(PROJECT_DIR, 'data', 'gnomad_v4_constraint.tsv')
ANCESTRY_LOEUF_CSV = os.path.join(PROJECT_DIR, 'data', 'gnomad_v4_ancestry_loeuf.csv')
OUT_CSV            = os.path.join(PROJECT_DIR, 'data', 'network_constraint_ancestry_loeuf.csv')
OUT_DIR            = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Column mapping: MOESM11 → internal names ──────────────────────────────
# NFE..n.440k. is the primary NFE column (largest cohort, most power).
# AMR and FIN are not present in MOESM11; those panels are dropped.
ANCESTRY_COL_MAP = {
    'NFE..n.440k.':              'LOEUF_nfe',
    'AFR':                       'LOEUF_afr',
    'EAS':                       'LOEUF_eas',
    'SAS':                       'LOEUF_sas',
    'ASJ':                       'LOEUF_asj',
    'Maximally.Diverse..n.43k.': 'LOEUF_max_diverse',
    'NFE..n.43k.':               'LOEUF_nfe_43k',
    'NFE..n.20k.':               'LOEUF_nfe_20k',
    'Full.Dataset..n.460k.':     'LOEUF_full',
}

# Sample sizes for axis / caption labels
POP_N = {
    'nfe': '440k', 'afr': '8,701', 'eas': '2,150',
    'sas': '9,217', 'asj': '2,671',
}

# ── Shared display constants ───────────────────────────────────────────────
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

LABEL_GENES = {'TYR', 'TYRP1', 'DCT', 'MC1R', 'OCA2', 'MITF', 'SOX10', 'KIT'}


# ═══════════════════════════════════════════════════════════════════════════
# Step 1: Load per-ancestry LOEUF (Han et al. 2025, Supp Data 11 / MOESM11)
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("Step 1: Loading per-ancestry LOEUF (MOESM11)")
print("=" * 70)

if not os.path.exists(ANCESTRY_LOEUF_CSV):
    sys.exit(
        f"ERROR: {ANCESTRY_LOEUF_CSV} not found.\n"
        "Copy ~/Downloads/41467_2025_57885_MOESM11_ESM.csv "
        "to data/gnomad_v4_ancestry_loeuf.csv"
    )

anc_raw = pd.read_csv(ANCESTRY_LOEUF_CSV)
# Strip leading apostrophe from gene symbols (e.g. "'A1BG" → "A1BG")
anc_raw['gene_upper'] = anc_raw['Gene'].str.strip("'").str.strip().str.upper()

# Rename to internal column names; drop any MOESM11 cols not in map
anc = anc_raw.rename(columns=ANCESTRY_COL_MAP)
anc_loeuf_cols = [v for v in ANCESTRY_COL_MAP.values() if v in anc.columns]

print(f"  Loaded {len(anc):,} genes")
print(f"  Available ancestry LOEUF columns: {anc_loeuf_cols}")


# ═══════════════════════════════════════════════════════════════════════════
# Step 2: Load global v4 LOEUF for sanity check (Comparison C)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 2: Loading global v4 LOEUF baseline")
print("=" * 70)

v4_global = None
if os.path.exists(V4_TSV):
    try:
        v4g = pd.read_csv(
            V4_TSV, sep='\t',
            usecols=['gene', 'lof.oe_ci.upper', 'canonical'],
            low_memory=False
        )
        v4g = v4g[v4g['canonical'].astype(str).str.lower() == 'true']
        v4g = (v4g.sort_values('lof.oe_ci.upper', na_position='last')
                  .drop_duplicates('gene', keep='first')
                  .rename(columns={'gene': 'gene_upper',
                                   'lof.oe_ci.upper': 'LOEUF_v4'}))
        v4g['gene_upper'] = v4g['gene_upper'].str.upper()
        v4_global = v4g[['gene_upper', 'LOEUF_v4']]
        print(f"  Loaded {len(v4_global):,} genes (canonical transcripts)")
    except Exception as e:
        print(f"  WARNING: Could not load v4 global file: {e}")
else:
    print("  v4 global file not found — Comparison C will be skipped")


# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Merge with 129-gene master network dataset
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 3: Merging with master 129-gene dataset")
print("=" * 70)

master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
print(f"  Master dataset: {len(master)} genes")

anc_slim = anc[['gene_upper'] + anc_loeuf_cols].rename(columns={'gene_upper': 'gene'})
df = master.merge(anc_slim, on='gene', how='left')

if v4_global is not None:
    df = df.merge(
        v4_global.rename(columns={'gene_upper': 'gene'}),
        on='gene', how='left'
    )

n_matched = df['LOEUF_nfe'].notna().sum() if 'LOEUF_nfe' in df.columns else 0
print(f"  Matched {n_matched}/{len(df)} network genes with per-ancestry LOEUF")
unmatched = df.loc[df['LOEUF_nfe'].isna(), 'gene'].tolist() \
    if 'LOEUF_nfe' in df.columns else []
if unmatched:
    print(f"  Unmatched genes: {unmatched}")

df.to_csv(OUT_CSV, index=False)
print(f"  Saved merged data → {OUT_CSV}")

# ── Compute non-European composite LOEUF (median of AFR, EAS, SAS per gene) ─
# ASJ excluded: Ashkenazi Jewish is a European-origin founder population.
NON_EUR_COLS = [c for c in ['LOEUF_afr', 'LOEUF_eas', 'LOEUF_sas'] if c in df.columns]
NON_EUR_COLS_ANC = [c for c in ['LOEUF_afr', 'LOEUF_eas', 'LOEUF_sas'] if c in anc.columns]
print(f"\n  Non-European composite from: {NON_EUR_COLS}")
df['LOEUF_non_eur']  = df[NON_EUR_COLS].median(axis=1)
anc['LOEUF_non_eur'] = anc[NON_EUR_COLS_ANC].median(axis=1)

# Outlier threshold: genes where |non_eur - nfe| exceeds 1.5 × IQR of the delta
# (computed across all network genes with both values)
_delta = (df['LOEUF_non_eur'] - df['LOEUF_nfe']).dropna()
_q1, _q3 = _delta.quantile(0.25), _delta.quantile(0.75)
OUTLIER_THRESHOLD = 1.5 * (_q3 - _q1)
df['_delta_non_eur'] = df['LOEUF_non_eur'] - df['LOEUF_nfe']
df['_is_outlier']    = df['_delta_non_eur'].abs() > OUTLIER_THRESHOLD
print(f"  Outlier threshold (1.5×IQR of delta): ±{OUTLIER_THRESHOLD:.3f}")
outlier_genes = df.loc[df['_is_outlier'], 'gene'].tolist()
print(f"  Outlier genes (n={len(outlier_genes)}): {outlier_genes}")


# ═══════════════════════════════════════════════════════════════════════════
# Step 4: Statistical comparisons
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 4: Statistical comparisons")
print("=" * 70)

pigment_genes = set(master.loc[master['functional_category'] == 'Pigment-specific', 'gene'])
print(f"  Pigmentation gene set: n={len(pigment_genes)}")

# ── Comparison A: Pigmentation genes vs genome-wide background ────────────
print("\nComparison A — Pigmentation set vs genome-wide background (Mann-Whitney U)")
for pop_label, col in [('European (NFE)', 'LOEUF_nfe'),
                        ('Non-European composite', 'LOEUF_non_eur')]:
    src = anc if col in anc.columns else df
    if col not in src.columns:
        continue
    gene_col = 'gene_upper' if 'gene_upper' in src.columns else 'gene'
    fg = src.loc[src[gene_col].isin(pigment_genes), col].dropna()
    bg = src.loc[~src[gene_col].isin(pigment_genes), col].dropna()
    if len(fg) == 0:
        print(f"  {pop_label}: no pigmentation genes matched")
        continue
    u_stat, p_mwu = stats.mannwhitneyu(fg, bg, alternative='two-sided')
    print(f"  {pop_label}: n_pig={len(fg)}, n_bg={len(bg):,}, "
          f"med_pig={fg.median():.3f}, med_bg={bg.median():.3f}, "
          f"U={u_stat:.0f}, p={p_mwu:.3e}")

# ── Comparison B: Paired European vs non-European within pigmentation set ──
print("\nComparison B — European vs non-European within pigmentation set (Wilcoxon)")
rho_na, p_na = np.nan, np.nan
if 'LOEUF_nfe' in df.columns and 'LOEUF_non_eur' in df.columns:
    pig_paired = df[df['functional_category'] == 'Pigment-specific'].dropna(
        subset=['LOEUF_nfe', 'LOEUF_non_eur'])
    print(f"  n pigmentation genes with both European and non-European LOEUF: {len(pig_paired)}")
    if len(pig_paired) >= 5:
        w_stat, w_p = stats.wilcoxon(pig_paired['LOEUF_nfe'], pig_paired['LOEUF_non_eur'])
        print(f"  Wilcoxon W = {w_stat:.0f}, p = {w_p:.3e}")
        print(f"  Median LOEUF — European: {pig_paired['LOEUF_nfe'].median():.3f}, "
              f"Non-European: {pig_paired['LOEUF_non_eur'].median():.3f}")
    else:
        print(f"  Too few paired genes for Wilcoxon (n={len(pig_paired)})")

    # Spearman across all network genes (used in scatter label)
    both = df.dropna(subset=['LOEUF_nfe', 'LOEUF_non_eur'])
    if len(both) >= 3:
        rho_na, p_na = stats.spearmanr(both['LOEUF_nfe'], both['LOEUF_non_eur'])
        print(f"\n  Spearman ρ(European, Non-European) across {len(both)} network genes: "
              f"ρ = {rho_na:.3f}, p = {p_na:.3e}")
else:
    print("  Required columns missing — skipping")

# ── Comparison C: v2.1.1 LOEUF vs v4 global LOEUF (sanity check) ─────────
print("\nComparison C — v2.1.1 LOEUF vs v4 global LOEUF (sanity check)")
if 'LOEUF_v4' in df.columns and 'LOEUF' in df.columns:
    ab = df.dropna(subset=['LOEUF', 'LOEUF_v4'])
    rho_ab, p_ab = stats.spearmanr(ab['LOEUF'], ab['LOEUF_v4'])
    print(f"  n = {len(ab)}, Spearman ρ = {rho_ab:.3f}, p = {p_ab:.3e}")
    print(f"  Median LOEUF: v2.1.1 = {ab['LOEUF'].median():.3f}, "
          f"v4 = {ab['LOEUF_v4'].median():.3f}")
else:
    print("  Skipping (v4 global LOEUF not available)")

# ── Comparison D: European vs non-European by functional category ──────────
print("\nComparison D — European vs non-European LOEUF by functional category")
if 'LOEUF_nfe' in df.columns and 'LOEUF_non_eur' in df.columns:
    d = df.dropna(subset=['LOEUF_nfe', 'LOEUF_non_eur'])
    print(f"\n  {'Category':<30} {'N':>4}  {'med_EUR':>8}  {'med_nonEUR':>10}  {'diff':>8}")
    print(f"  {'-'*30} {'-'*4}  {'-'*8}  {'-'*10}  {'-'*8}")
    for cat in CATEGORY_ORDER:
        sub = d[d['functional_category'] == cat]
        if len(sub) == 0:
            continue
        mn = sub['LOEUF_nfe'].median()
        ma = sub['LOEUF_non_eur'].median()
        print(f"  {cat:<30} {len(sub):>4}  {mn:>8.3f}  {ma:>10.3f}  {mn-ma:>+8.3f}")
else:
    print("  Required columns missing — skipping")

# ── Comparison E: MC1R across all ancestries ──────────────────────────────
print("\nComparison E — MC1R LOEUF across ancestries")
mc1r = df[df['gene'] == 'MC1R']
if len(mc1r) > 0:
    row = mc1r.iloc[0]
    for col in ['LOEUF', 'LOEUF_v4', 'LOEUF_nfe', 'LOEUF_non_eur',
                'LOEUF_afr', 'LOEUF_eas', 'LOEUF_sas', 'LOEUF_asj']:
        if col in df.columns:
            val = row.get(col, np.nan)
            label = col.replace('LOEUF_', '').upper() if col != 'LOEUF' else 'v2.1.1'
            print(f"  {label:<14}: {val:.3f}" if pd.notna(val) else f"  {label:<14}: N/A")
else:
    print("  MC1R not found in dataset")

# ── Comparison F: Pigmentation vs Generic signaling (EUR and non-EUR) ──────
print("\nComparison F — Pigmentation vs Generic signaling (European vs non-European)")
if 'LOEUF_nfe' in df.columns and 'LOEUF_non_eur' in df.columns:
    for pop_label, col in [('European', 'LOEUF_nfe'), ('Non-European', 'LOEUF_non_eur')]:
        pig = df.loc[df['functional_category'] == 'Pigment-specific', col].dropna()
        sig = df.loc[df['functional_category'] == 'Generic signaling', col].dropna()
        if len(pig) > 0 and len(sig) > 0:
            u, p = stats.mannwhitneyu(pig, sig, alternative='two-sided')
            print(f"  {pop_label} Pigment-specific vs Generic signaling: "
                  f"n_pig={len(pig)}, n_sig={len(sig)}, "
                  f"med_pig={pig.median():.3f}, med_sig={sig.median():.3f}, "
                  f"U={u:.0f}, p={p:.3e}")

# ── Print outlier gene details ─────────────────────────────────────────────
print("\nOutlier genes (|non-European − European| > 1.5×IQR):")
if outlier_genes:
    out_df = df.loc[df['_is_outlier'],
                    ['gene', 'functional_category', 'LOEUF_nfe',
                     'LOEUF_afr', 'LOEUF_eas', 'LOEUF_sas',
                     'LOEUF_non_eur', '_delta_non_eur']
                   ].sort_values('_delta_non_eur')
    print(out_df.to_string(index=False))
else:
    print("  None")


# ═══════════════════════════════════════════════════════════════════════════
# Step 5: Figure
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 5: Generating figure")
print("=" * 70)

if 'LOEUF_nfe' not in df.columns or 'LOEUF_non_eur' not in df.columns:
    print("WARNING: Required columns missing — skipping figure.")
    sys.exit(0)

df_plot = df.dropna(subset=['LOEUF_nfe', 'LOEUF_non_eur'])
print(f"  Plotting {len(df_plot)} genes")

fig, axes = plt.subplots(1, 2, figsize=(18, 8))
fig.subplots_adjust(wspace=0.34, left=0.07, right=0.98, top=0.88, bottom=0.40)

# ── Panel A: European vs non-European scatter with outlier labeling ────────
ax = axes[0]
for cat in CATEGORY_ORDER:
    sub = df_plot[df_plot['functional_category'] == cat]
    if len(sub) == 0:
        continue
    ax.scatter(sub['LOEUF_nfe'], sub['LOEUF_non_eur'],
               c=CATEGORY_COLORS[cat], s=65, alpha=0.85,
               edgecolors='white', linewidths=0.4,
               label=f'{cat} (n={len(sub)})', zorder=3)

xy_min = min(df_plot['LOEUF_nfe'].min(), df_plot['LOEUF_non_eur'].min()) - 0.05
xy_max = max(df_plot['LOEUF_nfe'].max(), df_plot['LOEUF_non_eur'].max()) + 0.05
ax.plot([xy_min, xy_max], [xy_min, xy_max], '--',
        color='#999999', alpha=0.6, lw=1.2, zorder=1, label='y = x')

# Label outliers and key pigmentation genes
label_set = LABEL_GENES | set(outlier_genes)
texts = []
for _, row in df_plot.iterrows():
    if row['gene'] in label_set:
        weight = 'bold' if row['_is_outlier'] else 'normal'
        color  = '#CC0000' if row['_is_outlier'] else '#333333'
        texts.append(ax.text(row['LOEUF_nfe'], row['LOEUF_non_eur'],
                             row['gene'], fontsize=11, fontweight=weight,
                             color=color, style='italic'))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7),
            force_points=(2.5, 2.5), force_text=(2.5, 2.5),
            expand_points=(3.0, 3.0), expand_text=(2.0, 2.0),
            iter_lim=1000)

rho_str = f"ρ = {rho_na:.3f}" if not np.isnan(rho_na) else "ρ = N/A"
p_str   = f"p = {p_na:.2e}" if not np.isnan(p_na) else ""
ax.text(0.98, -0.22, f'Spearman {rho_str},  {p_str}',
        transform=ax.transAxes, fontsize=13, ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='gray', alpha=0.9))
ax.text(0.02, -0.28,
        'Non-European = median(AFR, EAS, SAS); source: Han et al. 2025 Nat Commun, Supp Data 11\n'
        'Red labels = outliers (|Δ| > 1.5×IQR). NFE n=440k; AFR n=8,701; EAS n=2,150; SAS n=9,217',
        transform=ax.transAxes, fontsize=9, ha='left', va='top',
        color='#666666', style='italic')

ax.set_xlabel('LOEUF — European (NFE, n=440k)', fontsize=16)
ax.set_ylabel('LOEUF — Non-European (median AFR/EAS/SAS)', fontsize=16)
ax.set_title('European vs non-European LoF intolerance\n(outlier genes in red)', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.legend(fontsize=11, loc='upper center', bbox_to_anchor=(0.5, -0.34),
          ncol=2, framealpha=0.9, edgecolor='gray',
          handletextpad=0.4, borderpad=0.5, markerscale=1.2)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'A', transform=ax.transAxes, fontsize=24,
        fontweight='bold', va='top')

# ── Panel B: side-by-side boxplots European vs non-European by category ────
ax = axes[1]
n_cats      = len(CATEGORY_ORDER)
box_width   = 0.35
gap         = 0.08
cat_centers = np.arange(n_cats) * 1.3

eur_positions    = cat_centers - (box_width / 2 + gap / 2)
noneur_positions = cat_centers + (box_width / 2 + gap / 2)

eur_data, noneur_data, ns = [], [], []
for cat in CATEGORY_ORDER:
    sub = df_plot[df_plot['functional_category'] == cat]
    eur_data.append(sub['LOEUF_nfe'].dropna().values)
    noneur_data.append(sub['LOEUF_non_eur'].dropna().values)
    ns.append(len(sub))


def _make_bp(ax, data, positions, colors, alpha=0.85, hatch=None):
    bp = ax.boxplot(data, positions=positions, widths=box_width,
                    patch_artist=True, showfliers=False, zorder=2,
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(color='#555555'),
                    capprops=dict(color='#555555'))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(alpha)
        patch.set_edgecolor('gray')
        if hatch:
            patch.set_hatch(hatch)
    return bp


colors_list = [CATEGORY_COLORS[c] for c in CATEGORY_ORDER]
_make_bp(ax, eur_data,    eur_positions,    colors_list, alpha=0.90)
_make_bp(ax, noneur_data, noneur_positions, colors_list, alpha=0.55, hatch='//')

rng = np.random.default_rng(42)
for i, (ev, nv, cat) in enumerate(zip(eur_data, noneur_data, CATEGORY_ORDER)):
    c = CATEGORY_COLORS[cat]
    if len(ev) > 0:
        jitter = rng.uniform(-0.08, 0.08, len(ev))
        ax.scatter(eur_positions[i] + jitter, ev,
                   c=c, s=18, alpha=0.85, edgecolors='white',
                   linewidths=0.3, zorder=3)
    if len(nv) > 0:
        jitter = rng.uniform(-0.08, 0.08, len(nv))
        ax.scatter(noneur_positions[i] + jitter, nv,
                   c=c, s=18, alpha=0.55, edgecolors='white',
                   linewidths=0.3, zorder=3)

all_vals = [v for d in eur_data + noneur_data for v in d]
y_label  = min(all_vals) - 0.08 if all_vals else -0.1
for i, n in enumerate(ns):
    ax.text(cat_centers[i], y_label, f'n={n}', ha='center', fontsize=10,
            color='#555555', fontstyle='italic', zorder=4)

eur_patch    = mpatches.Patch(facecolor='#888888', alpha=0.90,
                               label='European (NFE, n=440k)')
noneur_patch = mpatches.Patch(facecolor='#888888', alpha=0.55, hatch='//',
                               label='Non-European (AFR/EAS/SAS median)')
ax.legend(handles=[eur_patch, noneur_patch], fontsize=12,
          loc='upper right', framealpha=0.9, edgecolor='gray')

ax.set_xticks(cat_centers)
ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                   fontsize=12, rotation=30, ha='right')
ax.set_ylabel('LOEUF (lower = more constrained)', fontsize=16)
ax.set_title('LOEUF by functional category:\nEuropean vs non-European', fontsize=16,
             fontweight='bold', loc='left', pad=10)
ax.tick_params(labelsize=13)
ax.text(-0.08, 1.08, 'B', transform=ax.transAxes, fontsize=24,
        fontweight='bold', va='top')

for ext in ('png', 'pdf'):
    out_path = os.path.join(OUT_DIR, f'figure_phase1_ancestry_loeuf.{ext}')
    fig.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved → {out_path}")
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Step 6: Summary tables
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 6: Summary tables")
print("=" * 70)

for pop_label, col in [('European (NFE)', 'LOEUF_nfe'),
                        ('Non-European composite', 'LOEUF_non_eur')]:
    if col not in df.columns:
        continue
    print(f"\nLOEUF summary by functional category ({pop_label}):")
    s = (df.dropna(subset=[col])
           .groupby('functional_category')[col]
           .agg(['count', 'median', 'mean'])
           .round(3)
           .rename(columns={'count': 'N', 'median': 'Median', 'mean': 'Mean'}))
    s = s.loc[[c for c in CATEGORY_ORDER if c in s.index]]
    print(s.to_string())

print("\nDone!")
