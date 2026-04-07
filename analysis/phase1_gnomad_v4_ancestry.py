"""
phase1_gnomad_v4_ancestry.py

Phase 1.3: gnomAD v4.1 ancestry-stratified LOEUF for the Raghunath melanogenesis network.

Tests whether the constraint architecture (pigment-specific = tolerant, signaling =
constrained) found in v2.1.1 (77% European) is robust to ancestry, using per-ancestry
LOEUF from gnomAD v4.1.

Downloads gnomad.v4.1.constraint_metrics.tsv from Google Cloud Storage, extracts
per-ancestry LOEUF, merges with our 129-gene network, and produces:

  1. data/network_constraint_ancestry_loeuf.csv  — merged dataset with v4 ancestry LOEUF
  2. output/figure_phase1_ancestry_loeuf.png/pdf — two-panel figure:
       Panel A: NFE LOEUF vs AFR LOEUF scatter (per gene, colored by category)
       Panel B: LOEUF by functional category, NFE vs AFR side-by-side boxplots

Scientific motivation:
  v2.1.1 LOEUF is computed on ~77% European individuals. If pigment-specific genes still
  show high LOEUF (tolerance) in African-ancestry data, the biology is robust: heterozygous
  LoF has no fitness cost because these are recessive genes. If MC1R's LOEUF drops in
  African-ancestry data, the v2.1.1 estimate is inflated by European-specific functional
  variants (red hair alleles).

gnomAD v4.1 source:
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/
  gnomad.v4.1.constraint_metrics.tsv
"""

import os
import sys
import urllib.request
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
PROJECT_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV    = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
V4_TSV        = os.path.join(PROJECT_DIR, 'data', 'gnomad_v4_constraint.tsv')
OUT_CSV       = os.path.join(PROJECT_DIR, 'data', 'network_constraint_ancestry_loeuf.csv')
OUT_DIR       = os.path.join(PROJECT_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

V4_URL = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
    "gnomad.v4.1.constraint_metrics.tsv"
)

# ── Shared color/label constants ───────────────────────────────────────────
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

# Population codes to try (gnomAD v4 uses these prefixes)
POP_CODES = ['afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'mid', 'sas']


# ═══════════════════════════════════════════════════════════════════════════
# Step 1: Download gnomAD v4.1 constraint file
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("Step 1: gnomAD v4.1 constraint file")
print("=" * 70)

if os.path.exists(V4_TSV):
    print(f"File already present: {V4_TSV}")
else:
    print(f"Downloading gnomAD v4.1 constraint metrics...")
    print(f"  URL: {V4_URL}")
    print("  (This is a large file — ~500 MB. Please wait...)")
    try:
        def _progress(count, block_size, total_size):
            if total_size > 0 and count % 500 == 0:
                pct = min(count * block_size / total_size * 100, 100)
                print(f"  {pct:.1f}%", end='\r', flush=True)

        urllib.request.urlretrieve(V4_URL, V4_TSV, reporthook=_progress)
        print(f"\n  Saved → {V4_TSV}")
    except Exception as e:
        print(f"\nERROR: Could not download gnomAD v4.1 file: {e}")
        print(
            "\nAlternative URLs to try:\n"
            "  https://gnomad.broadinstitute.org/downloads#v4-constraint\n"
            "  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/\n"
            f"\nIf you download it manually, save to: {V4_TSV}\n"
            "Then re-run this script."
        )
        sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════
# Step 2: Inspect column names
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 2: Inspecting gnomAD v4.1 column names")
print("=" * 70)

# Read only header row first — avoids loading the full file
header_df = pd.read_csv(V4_TSV, sep='\t', nrows=0)
all_cols = list(header_df.columns)

print(f"\nTotal columns: {len(all_cols)}")
print("\nAll columns:")
for c in all_cols:
    print(f"  {c}")

# Find LoF / o/e related columns
lof_cols = [c for c in all_cols if any(k in c.lower() for k in ('lof', 'oe', 'o/e'))]
print(f"\nColumns containing 'lof' or 'oe' ({len(lof_cols)} total):")
for c in lof_cols:
    print(f"  {c}")


# ── Detect LOEUF column naming convention ─────────────────────────────────
# gnomAD v4 typically uses dot-notation: lof.oe_ci.upper
# Per-ancestry: {pop}.lof.oe_ci.upper  OR  lof.oe_ci.upper_{pop}
# We try both patterns and report what we find.

def find_loeuf_col(cols, pop=None):
    """
    Search for the LOEUF column for a given population (or overall if pop=None).
    Returns the column name or None.
    """
    candidates_dot   = []
    candidates_under = []
    if pop is None:
        # Overall — avoid columns with any pop prefix/suffix
        for c in cols:
            cl = c.lower()
            if ('lof' in cl and ('oe_ci' in cl or 'oe_upper' in cl or 'upper' in cl)):
                # Exclude ancestry-specific
                if not any(f'{p}.' in cl or f'.{p}' in cl or f'_{p}' in cl
                           for p in POP_CODES):
                    candidates_dot.append(c)
        # Also try exact known v4 name
        for name in ('lof.oe_ci.upper', 'oe_lof_upper', 'LOEUF'):
            if name in cols:
                return name
    else:
        p = pop.lower()
        for c in cols:
            cl = c.lower()
            if p in cl and ('lof' in cl) and ('upper' in cl or 'oe_ci' in cl):
                candidates_dot.append(c)
        # Try exact patterns
        for pattern in (
            f'{p}.lof.oe_ci.upper',
            f'lof.oe_ci.upper_{p}',
            f'{p}_lof_oe_upper',
            f'oe_lof_upper_{p}',
        ):
            if pattern in cols:
                return pattern

    all_candidates = candidates_dot + candidates_under
    return all_candidates[0] if all_candidates else None


# Identify overall LOEUF column
overall_col = find_loeuf_col(all_cols)
if overall_col is None:
    # Last resort: any column with 'upper' and 'lof'
    for c in all_cols:
        if 'lof' in c.lower() and 'upper' in c.lower():
            overall_col = c
            break

print(f"\nOverall LOEUF column: {overall_col!r}")

# Identify per-ancestry LOEUF columns
pop_loeuf_cols = {}
for pop in POP_CODES:
    col = find_loeuf_col(all_cols, pop=pop)
    if col:
        pop_loeuf_cols[pop] = col
        print(f"  {pop.upper()} LOEUF column: {col!r}")
    else:
        print(f"  {pop.upper()} LOEUF column: NOT FOUND")

# Minimal required: nfe + afr
missing_required = [p for p in ('nfe', 'afr') if p not in pop_loeuf_cols]
if missing_required:
    print(f"\nWARNING: Could not auto-detect LOEUF columns for: {missing_required}")
    print("Check the column listing above and update pop_loeuf_cols manually.")


# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Load gnomAD v4.1 — only the columns we need
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 3: Loading gnomAD v4.1 data (relevant columns only)")
print("=" * 70)

# Identify gene identifier columns
gene_sym_candidates = [c for c in all_cols if c.lower() in ('gene', 'gene_name', 'symbol')]
gene_id_candidates  = [c for c in all_cols if c.lower() in ('gene_id', 'ensembl_gene_id')]
gene_sym_col = gene_sym_candidates[0] if gene_sym_candidates else None
gene_id_col  = gene_id_candidates[0]  if gene_id_candidates  else None

print(f"Gene symbol column: {gene_sym_col!r}")
print(f"Gene ID column:     {gene_id_col!r}")

# Build column list to read
use_cols = []
for c in (gene_sym_col, gene_id_col, overall_col):
    if c and c not in use_cols:
        use_cols.append(c)
for col in pop_loeuf_cols.values():
    if col not in use_cols:
        use_cols.append(col)

# Also grab transcript column if present (for filtering canonical transcripts)
transcript_col = next((c for c in all_cols if 'transcript' in c.lower()), None)
canonical_col  = next(
    (c for c in all_cols if 'canonical' in c.lower() or 'mane' in c.lower()), None
)
if transcript_col and transcript_col not in use_cols:
    use_cols.append(transcript_col)
if canonical_col and canonical_col not in use_cols:
    use_cols.append(canonical_col)

print(f"\nReading {len(use_cols)} columns from gnomAD v4.1...")
print(f"  Columns: {use_cols}")

v4 = pd.read_csv(V4_TSV, sep='\t', usecols=use_cols, low_memory=False)
print(f"  Loaded: {len(v4):,} rows")


# ── Filter to canonical / MANE transcripts ─────────────────────────────────
if canonical_col:
    before = len(v4)
    v4 = v4[v4[canonical_col].astype(str).str.lower().isin(('true', '1', 'yes'))]
    print(f"  After canonical filter ({canonical_col}): {len(v4):,} rows (dropped {before - len(v4):,})")
elif transcript_col:
    # Fall back: keep one row per gene (row with lowest overall LOEUF — most constrained
    # transcript is usually the canonical one in gnomAD)
    if overall_col and overall_col in v4.columns:
        before = len(v4)
        v4 = (v4.sort_values(overall_col, ascending=True, na_position='last')
                .drop_duplicates(subset=gene_sym_col or gene_id_col, keep='first')
                .reset_index(drop=True))
        print(f"  Deduplicated to one row/gene (lowest overall LOEUF): {len(v4):,} rows (was {before:,})")
    else:
        v4 = v4.drop_duplicates(subset=gene_sym_col or gene_id_col).reset_index(drop=True)
        print(f"  Deduplicated: {len(v4):,} rows")


# ── Normalise gene symbols ─────────────────────────────────────────────────
if gene_sym_col:
    v4['gene_upper'] = v4[gene_sym_col].astype(str).str.upper()
elif gene_id_col:
    v4['gene_upper'] = v4[gene_id_col].astype(str).str.upper()
else:
    raise RuntimeError("Could not identify a gene identifier column in v4 file.")


# ── Rename columns to standardised names ─────────────────────────────────
rename_map = {}
if overall_col:
    rename_map[overall_col] = 'LOEUF_v4'
for pop, col in pop_loeuf_cols.items():
    rename_map[col] = f'LOEUF_{pop}'
v4 = v4.rename(columns=rename_map)

print(f"\nSample of v4 data (first 3 rows):")
display_cols = ['gene_upper'] + [c for c in v4.columns if c.startswith('LOEUF')]
print(v4[display_cols].head(3).to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════
# Step 4: Merge with 129-gene master dataset
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 4: Merging with master 129-gene dataset")
print("=" * 70)

master = pd.read_csv(MASTER_CSV)
master['gene'] = master['gene'].str.upper()
print(f"Master dataset: {len(master)} genes")

v4_loeuf_cols = [c for c in v4.columns if c.startswith('LOEUF')]
merge_cols = ['gene_upper'] + v4_loeuf_cols
v4_slim = v4[merge_cols].rename(columns={'gene_upper': 'gene'})

df = master.merge(v4_slim, on='gene', how='left')

n_matched = df[v4_loeuf_cols[0] if v4_loeuf_cols else 'LOEUF_v4'].notna().sum() \
            if v4_loeuf_cols else 0
print(f"Matched: {n_matched} / {len(df)} genes")

unmatched = df.loc[df[v4_loeuf_cols[0]].isna() if v4_loeuf_cols else df.index, 'gene'].tolist()
if unmatched:
    print(f"Unmatched genes: {unmatched}")

df.to_csv(OUT_CSV, index=False)
print(f"Saved merged data → {OUT_CSV}")


# ═══════════════════════════════════════════════════════════════════════════
# Step 5: Statistical comparisons
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 5: Statistical comparisons")
print("=" * 70)

# ── Comparison A: v2.1.1 vs. v4 overall LOEUF ─────────────────────────────
if 'LOEUF_v4' in df.columns:
    print("\nComparison A — v2.1.1 LOEUF vs. v4 overall LOEUF (sanity check)")
    ab = df.dropna(subset=['LOEUF', 'LOEUF_v4'])
    rho_ab, p_ab = stats.spearmanr(ab['LOEUF'], ab['LOEUF_v4'])
    print(f"  n = {len(ab)}")
    print(f"  Spearman ρ = {rho_ab:.3f}, p = {p_ab:.3e}")
    v2_med = ab['LOEUF'].median()
    v4_med = ab['LOEUF_v4'].median()
    print(f"  Median LOEUF: v2.1.1 = {v2_med:.3f}, v4 = {v4_med:.3f}")
else:
    print("\nComparison A: LOEUF_v4 column not available — skipping.")

# ── Comparison B: NFE vs. AFR by category ─────────────────────────────────
print("\nComparison B — NFE vs. AFR LOEUF by functional category")
if 'LOEUF_nfe' in df.columns and 'LOEUF_afr' in df.columns:
    b_cols = ['gene', 'functional_category', 'LOEUF', 'LOEUF_nfe', 'LOEUF_afr']
    b = df[b_cols].dropna(subset=['LOEUF_nfe', 'LOEUF_afr'])
    print(f"  n with both NFE and AFR: {len(b)}")
    print(f"\n  {'Category':<30} {'N':>4}  {'med_NFE':>8}  {'med_AFR':>8}  {'diff':>8}")
    print(f"  {'-'*30} {'-'*4}  {'-'*8}  {'-'*8}  {'-'*8}")
    for cat in CATEGORY_ORDER:
        sub = b[b['functional_category'] == cat]
        if len(sub) == 0:
            continue
        mn = sub['LOEUF_nfe'].median()
        ma = sub['LOEUF_afr'].median()
        print(f"  {cat:<30} {len(sub):>4}  {mn:>8.3f}  {ma:>8.3f}  {mn-ma:>+8.3f}")
else:
    missing = [p for p in ('LOEUF_nfe', 'LOEUF_afr') if p not in df.columns]
    print(f"  Skipping — missing columns: {missing}")
    b = pd.DataFrame()

# ── Comparison C: MC1R across all ancestries ──────────────────────────────
print("\nComparison C — MC1R LOEUF across all ancestries")
mc1r = df[df['gene'] == 'MC1R']
if len(mc1r) > 0:
    mc1r_row = mc1r.iloc[0]
    loeuf_cols_all = ['LOEUF'] + [f'LOEUF_{p}' for p in ['v4'] + POP_CODES]
    loeuf_cols_present = [c for c in loeuf_cols_all if c in df.columns]
    print(f"  {'Column':<18}  {'LOEUF':>8}")
    print(f"  {'-'*18}  {'-'*8}")
    for col in loeuf_cols_present:
        val = mc1r_row.get(col, np.nan)
        val_str = f"{val:.3f}" if pd.notna(val) else "N/A"
        label = col.replace('LOEUF_', '').upper() if col != 'LOEUF' else 'v2.1.1 overall'
        print(f"  {label:<18}  {val_str:>8}")
else:
    print("  MC1R not found in dataset.")

# ── Comparison D: Mann-Whitney U tests ────────────────────────────────────
print("\nComparison D — Mann-Whitney U: Pigment-specific vs. Generic signaling")
for pop_label, loeuf_col in [
    ('v2.1.1 overall', 'LOEUF'),
    ('v4 NFE',         'LOEUF_nfe'),
    ('v4 AFR',         'LOEUF_afr'),
]:
    if loeuf_col not in df.columns:
        continue
    sub = df.dropna(subset=[loeuf_col])
    pig = sub.loc[sub['functional_category'] == 'Pigment-specific', loeuf_col].values
    sig = sub.loc[sub['functional_category'] == 'Generic signaling', loeuf_col].values
    if len(pig) == 0 or len(sig) == 0:
        print(f"  {pop_label}: insufficient data")
        continue
    u_stat, p_mwu = stats.mannwhitneyu(pig, sig, alternative='two-sided')
    print(f"  {pop_label:<18}: n_pig={len(pig)}, n_sig={len(sig)}, "
          f"U={u_stat:.0f}, p={p_mwu:.3e}  "
          f"(med_pig={np.median(pig):.3f}, med_sig={np.median(sig):.3f})")

# NFE vs AFR concordance across all genes
if 'LOEUF_nfe' in df.columns and 'LOEUF_afr' in df.columns:
    both = df.dropna(subset=['LOEUF_nfe', 'LOEUF_afr'])
    rho_na, p_na = stats.spearmanr(both['LOEUF_nfe'], both['LOEUF_afr'])
    print(f"\n  Spearman ρ(LOEUF_nfe, LOEUF_afr) across {len(both)} genes: "
          f"ρ = {rho_na:.3f}, p = {p_na:.3e}")
else:
    rho_na, p_na = np.nan, np.nan


# ═══════════════════════════════════════════════════════════════════════════
# Step 6: Figure
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 6: Generating figure")
print("=" * 70)

# Check we have the required columns for plotting
has_nfe = 'LOEUF_nfe' in df.columns
has_afr = 'LOEUF_afr' in df.columns

if not (has_nfe and has_afr):
    print("WARNING: NFE and/or AFR LOEUF columns missing — skipping figure.")
    print(f"  Available LOEUF columns: {[c for c in df.columns if 'LOEUF' in c]}")
else:
    df_plot = df.dropna(subset=['LOEUF_nfe', 'LOEUF_afr'])
    print(f"  Plotting {len(df_plot)} genes with both NFE and AFR LOEUF")

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    fig.subplots_adjust(wspace=0.34, left=0.07, right=0.98, top=0.88, bottom=0.38)

    # ── Panel A: NFE LOEUF vs AFR LOEUF scatter ────────────────────────────
    ax = axes[0]
    for cat in CATEGORY_ORDER:
        sub = df_plot[df_plot['functional_category'] == cat]
        if len(sub) == 0:
            continue
        ax.scatter(sub['LOEUF_nfe'], sub['LOEUF_afr'],
                   c=CATEGORY_COLORS[cat], s=65, alpha=0.85,
                   edgecolors='white', linewidths=0.4,
                   label=f'{cat} (n={len(sub)})', zorder=3)

    # Diagonal y=x reference line
    xy_min = min(df_plot['LOEUF_nfe'].min(), df_plot['LOEUF_afr'].min()) - 0.05
    xy_max = max(df_plot['LOEUF_nfe'].max(), df_plot['LOEUF_afr'].max()) + 0.05
    ax.plot([xy_min, xy_max], [xy_min, xy_max], '--',
            color='#999999', alpha=0.6, lw=1.2, zorder=1, label='y = x')

    # Gene labels
    texts = []
    for _, row in df_plot.iterrows():
        if row['gene'] in LABEL_GENES:
            texts.append(ax.text(row['LOEUF_nfe'], row['LOEUF_afr'],
                                 row['gene'], fontsize=12, fontweight='bold',
                                 color='#333333', style='italic'))
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.7),
                force_points=(2.5, 2.5), force_text=(2.5, 2.5),
                expand_points=(3.0, 3.0), expand_text=(2.0, 2.0),
                iter_lim=1000)

    rho_str = f"ρ = {rho_na:.3f}" if not np.isnan(rho_na) else "ρ = N/A"
    p_str   = f"p = {p_na:.2e}" if not np.isnan(p_na) else ""
    ax.text(0.98, -0.20, f'Spearman {rho_str},  {p_str}',
            transform=ax.transAxes, fontsize=13, ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='gray', alpha=0.9))

    ax.set_xlabel('LOEUF (Non-Finnish European)', fontsize=16)
    ax.set_ylabel('LOEUF (African)', fontsize=16)
    ax.set_title('Ancestry-stratified LoF intolerance:\nNFE vs. AFR', fontsize=16,
                 fontweight='bold', loc='left', pad=10)
    ax.legend(fontsize=11, loc='upper center', bbox_to_anchor=(0.5, -0.32),
              ncol=2, framealpha=0.9, edgecolor='gray',
              handletextpad=0.4, borderpad=0.5, markerscale=1.2)
    ax.tick_params(labelsize=13)
    ax.text(-0.08, 1.08, 'A', transform=ax.transAxes, fontsize=24,
            fontweight='bold', va='top')

    # ── Panel B: side-by-side boxplots by category (NFE vs AFR) ──────────
    ax = axes[1]
    n_cats = len(CATEGORY_ORDER)
    group_width = 1.0
    box_width   = 0.35
    gap         = 0.08
    cat_centers = np.arange(n_cats) * (group_width + 0.3)

    nfe_positions = cat_centers - (box_width / 2 + gap / 2)
    afr_positions = cat_centers + (box_width / 2 + gap / 2)

    nfe_data, afr_data, ns = [], [], []
    for cat in CATEGORY_ORDER:
        sub = df_plot[df_plot['functional_category'] == cat]
        nfe_data.append(sub['LOEUF_nfe'].dropna().values)
        afr_data.append(sub['LOEUF_afr'].dropna().values)
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
    bp_nfe = _make_bp(ax, nfe_data, nfe_positions, colors_list, alpha=0.90)
    bp_afr = _make_bp(ax, afr_data, afr_positions, colors_list, alpha=0.55, hatch='//')

    # Jitter points
    rng = np.random.default_rng(42)
    for i, (nv, av, cat) in enumerate(zip(nfe_data, afr_data, CATEGORY_ORDER)):
        c = CATEGORY_COLORS[cat]
        if len(nv) > 0:
            jitter = rng.uniform(-0.08, 0.08, len(nv))
            ax.scatter(nfe_positions[i] + jitter, nv,
                       c=c, s=18, alpha=0.85, edgecolors='white',
                       linewidths=0.3, zorder=3)
        if len(av) > 0:
            jitter = rng.uniform(-0.08, 0.08, len(av))
            ax.scatter(afr_positions[i] + jitter, av,
                       c=c, s=18, alpha=0.55, edgecolors='white',
                       linewidths=0.3, zorder=3)

    # Sample size labels
    y_bottom = ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else -0.05
    all_vals = [v for d in nfe_data + afr_data for v in d]
    y_label = min(all_vals) - 0.08 if all_vals else -0.1
    for i, n in enumerate(ns):
        ax.text(cat_centers[i], y_label, f'n={n}', ha='center', fontsize=10,
                color='#555555', fontstyle='italic', zorder=4)

    # Legend for NFE / AFR
    nfe_patch = mpatches.Patch(facecolor='#888888', alpha=0.9, label='NFE (Non-Finnish European)')
    afr_patch = mpatches.Patch(facecolor='#888888', alpha=0.55, hatch='//', label='AFR (African)')
    ax.legend(handles=[nfe_patch, afr_patch], fontsize=12,
              loc='upper right', framealpha=0.9, edgecolor='gray')

    ax.set_xticks(cat_centers)
    ax.set_xticklabels([CATEGORY_SHORT[c] for c in CATEGORY_ORDER],
                       fontsize=12, rotation=30, ha='right')
    ax.set_ylabel('LOEUF (lower = more constrained)', fontsize=16)
    ax.set_title('LOEUF by functional category:\nNFE vs. AFR comparison', fontsize=16,
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
# Step 7: Summary tables for Quarto
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("Step 7: Summary tables")
print("=" * 70)

print("\nLOEUF summary by functional category (NFE):")
if 'LOEUF_nfe' in df.columns:
    s = (df.dropna(subset=['LOEUF_nfe'])
           .groupby('functional_category')['LOEUF_nfe']
           .agg(['count', 'median', 'mean'])
           .round(3)
           .rename(columns={'count': 'N', 'median': 'Median', 'mean': 'Mean'}))
    s = s.loc[[c for c in CATEGORY_ORDER if c in s.index]]
    print(s.to_string())

print("\nLOEUF summary by functional category (AFR):")
if 'LOEUF_afr' in df.columns:
    s = (df.dropna(subset=['LOEUF_afr'])
           .groupby('functional_category')['LOEUF_afr']
           .agg(['count', 'median', 'mean'])
           .round(3)
           .rename(columns={'count': 'N', 'median': 'Median', 'mean': 'Mean'}))
    s = s.loc[[c for c in CATEGORY_ORDER if c in s.index]]
    print(s.to_string())

print("\nDone!")
