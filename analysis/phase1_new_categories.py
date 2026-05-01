"""
phase1_new_categories.py

Build two new gene categorizations alongside the existing hand-curated
`functional_category` from data/LOEUF_by_functional_category.xlsx:

1. `gtex_tissue_category`  — derived from GTEx v8 expression (data-driven)
2. `kegg_primary_pathway`  — derived from KEGG pathway membership (database-driven)

Both columns are added to data/network_constraint_categorized.csv (built on top
of data/network_constraint_gtex.csv) so downstream analyses can compare
hand-curated vs data-driven groupings.

────────────────────────────────────────────────────────────────────────────
Categorization rules (documented here AND mirrored in pages/gtex.qmd)
────────────────────────────────────────────────────────────────────────────

GTEx tissue category (per gene):

  Compute on log2(TPM + 1) across all 54 GTEx v8 tissues:
    Tau (τ) = Σ (1 - x_i / max(x))  /  (n - 1)            Yanai et al. 2005
    n_expr  = # tissues with median TPM > 1
    max_t   = tissue with highest TPM

  Decision tree (first match wins):
    Housekeeping            : tau < 0.4  AND  n_expr ≥ 40
    Skin-restricted         : tau ≥ 0.6  AND  max_t ∈ skin tissues
    Brain-restricted        : tau ≥ 0.6  AND  max_t ∈ CNS/nerve/pituitary
    Reproductive-restricted : tau ≥ 0.6  AND  max_t ∈ gonads/uterus/cervix/prostate
    Immune-restricted       : tau ≥ 0.6  AND  max_t ∈ blood/spleen/lymphocytes
    Liver-restricted        : tau ≥ 0.6  AND  max_t == 'Liver'
    Other-restricted        : tau ≥ 0.6  AND  any other tissue
    Broad                   : everything else (intermediate τ)

  Thresholds (τ ≥ 0.6, n_expr ≥ 40) follow conventions in the tissue-
  specificity literature (Sonawane 2017, Kryuchkova-Mostacci 2017). Adjust
  here if a different stringency is desired.

KEGG primary pathway (per gene):

  KEGG link/pathway/hsa was fetched in fetch_kegg_pathways.py. Each gene's
  full pathway list is in data/kegg_pathway_lists.csv. We assign a single
  "primary pathway" using the first match from this priority list:

    1. Melanogenesis              (hsa04916)
    2. MAPK signaling             (hsa04010)
    3. PI3K-Akt signaling         (hsa04151)
    4. Apoptosis                  (hsa04210)
    5. Cytokine-cytokine receptor (hsa04060)
    6. Wnt signaling              (hsa04310)
    7. JAK-STAT signaling         (hsa04630)
    8. NF-kB signaling            (hsa04064)
    9. Other (gene is in KEGG but none of the above)
   10. Not in KEGG (no entry)

  Priority ordered so pigmentation-specific (Melanogenesis) wins over generic
  signaling cascades. Genes co-membered in many pathways pick the highest-
  priority match.
"""

import os
import gzip
import numpy as np
import pandas as pd

PROJECT_DIR  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GTEX_GCT     = os.path.join(PROJECT_DIR, 'data', 'GTEx_v8_gene_median_tpm.gct.gz')
NETWORK_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
KEGG_LISTS   = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_lists.csv')
OUT_CSV      = os.path.join(PROJECT_DIR, 'data', 'network_constraint_categorized.csv')

TPM_THRESHOLD = 1.0

# ── Tissue groupings for GTEx-based categorization ────────────────────────
SKIN_TISSUES = [
    'Skin - Sun Exposed (Lower leg)',
    'Skin - Not Sun Exposed (Suprapubic)',
]
BRAIN_TISSUES = [
    'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)',
    'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere',
    'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)',
    'Brain - Hippocampus', 'Brain - Hypothalamus',
    'Brain - Nucleus accumbens (basal ganglia)',
    'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
    'Brain - Substantia nigra', 'Pituitary', 'Nerve - Tibial',
]
REPRO_TISSUES = [
    'Ovary', 'Uterus', 'Vagina',
    'Cervix - Ectocervix', 'Cervix - Endocervix',
    'Fallopian Tube', 'Testis', 'Prostate',
]
IMMUNE_TISSUES = ['Whole Blood', 'Spleen', 'Cells - EBV-transformed lymphocytes']
LIVER_TISSUES  = ['Liver']

GTEX_CATEGORY_ORDER = [
    'Housekeeping', 'Skin-restricted', 'Brain-restricted',
    'Reproductive-restricted', 'Immune-restricted', 'Liver-restricted',
    'Other-restricted', 'Broad',
]

# ── KEGG pathway priority for primary-pathway assignment ──────────────────
KEGG_PRIORITY = [
    ('path:hsa04916', 'Melanogenesis'),
    ('path:hsa04010', 'MAPK signaling'),
    ('path:hsa04151', 'PI3K-Akt signaling'),
    ('path:hsa04210', 'Apoptosis'),
    ('path:hsa04060', 'Cytokine-cytokine receptor'),
    ('path:hsa04310', 'Wnt signaling'),
    ('path:hsa04630', 'JAK-STAT signaling'),
    ('path:hsa04064', 'NF-kB signaling'),
]
KEGG_CATEGORY_ORDER = [name for _, name in KEGG_PRIORITY] + ['Other (in KEGG)', 'Not in KEGG']

# ──────────────────────────────────────────────────────────────────────────
# Step 1: Load GTEx, compute per-gene tau / n_expr / max_tissue
# ──────────────────────────────────────────────────────────────────────────
print("Parsing GTEx GCT...")
with gzip.open(GTEX_GCT, 'rt') as f:
    f.readline(); f.readline()
    gtex_df = pd.read_csv(f, sep='\t')

TISSUES = [c for c in gtex_df.columns if c not in ('Name', 'Description')]
gtex_df['gene'] = gtex_df['Description'].str.upper()
tpm = gtex_df.groupby('gene')[TISSUES].max().reset_index()

# IL8 alias
if 'CXCL8' in tpm['gene'].values and 'IL8' not in tpm['gene'].values:
    cxcl8 = tpm[tpm['gene'] == 'CXCL8'].copy()
    cxcl8['gene'] = 'IL8'
    tpm = pd.concat([tpm, cxcl8], ignore_index=True)

log_tpm = np.log2(tpm[TISSUES].values + 1)
row_max = log_tpm.max(axis=1, keepdims=True)
safe = np.where(row_max == 0, 1, row_max)
tau = np.sum(1 - log_tpm / safe, axis=1) / (len(TISSUES) - 1)
tau = np.where(row_max.flatten() == 0, np.nan, tau)
tpm['tau']         = tau
tpm['n_expr']      = (tpm[TISSUES] > TPM_THRESHOLD).sum(axis=1)
tpm['max_tissue']  = tpm[TISSUES].idxmax(axis=1)


def gtex_category(row):
    if pd.isna(row['tau']):
        return 'Broad'  # zero-expression edge case
    if row['tau'] < 0.4 and row['n_expr'] >= 40:
        return 'Housekeeping'
    if row['tau'] >= 0.6:
        m = row['max_tissue']
        if m in SKIN_TISSUES:   return 'Skin-restricted'
        if m in BRAIN_TISSUES:  return 'Brain-restricted'
        if m in REPRO_TISSUES:  return 'Reproductive-restricted'
        if m in IMMUNE_TISSUES: return 'Immune-restricted'
        if m in LIVER_TISSUES:  return 'Liver-restricted'
        return 'Other-restricted'
    return 'Broad'


tpm['gtex_tissue_category'] = tpm.apply(gtex_category, axis=1)

# ──────────────────────────────────────────────────────────────────────────
# Step 2: Load KEGG pathway lists, assign primary pathway per gene
# ──────────────────────────────────────────────────────────────────────────
print("Loading KEGG pathway memberships...")
kegg = pd.read_csv(KEGG_LISTS)
kegg['gene'] = kegg['gene'].str.upper()
gene_pathways = kegg.groupby('gene')['pathway_id'].apply(set).to_dict()
genes_in_kegg = set(kegg[kegg['kegg_id'] != '']['gene'])


def kegg_category(gene):
    g = gene.upper()
    pathways = gene_pathways.get(g, set())
    for pid, name in KEGG_PRIORITY:
        if pid in pathways:
            return name
    if g in genes_in_kegg and pathways:
        return 'Other (in KEGG)'
    return 'Not in KEGG'


# ──────────────────────────────────────────────────────────────────────────
# Step 3: Merge with network, save
# ──────────────────────────────────────────────────────────────────────────
print("Merging into network CSV...")
net = pd.read_csv(NETWORK_CSV)
net['gene'] = net['gene'].str.upper()

merged = net.merge(
    tpm[['gene', 'tau', 'n_expr', 'max_tissue', 'gtex_tissue_category']],
    on='gene', how='left',
)
merged['kegg_primary_pathway'] = merged['gene'].apply(kegg_category)
merged.to_csv(OUT_CSV, index=False)
print(f"Saved → {OUT_CSV}  ({len(merged)} genes)")

# ──────────────────────────────────────────────────────────────────────────
# Step 4: Report breakdown
# ──────────────────────────────────────────────────────────────────────────
print("\nGTEx tissue category — network gene counts:")
counts_g = merged['gtex_tissue_category'].value_counts()
for cat in GTEX_CATEGORY_ORDER:
    if cat in counts_g.index:
        print(f"  {cat:25s}  {counts_g[cat]:>4d}")

print("\nKEGG primary pathway — network gene counts:")
counts_k = merged['kegg_primary_pathway'].value_counts()
for cat in KEGG_CATEGORY_ORDER:
    if cat in counts_k.index:
        print(f"  {cat:30s}  {counts_k[cat]:>4d}")

print("\nCross-tab (rows = GTEx category, cols = KEGG primary pathway):")
ct = pd.crosstab(merged['gtex_tissue_category'], merged['kegg_primary_pathway'])
print(ct.to_string())

print("\nDone!")
