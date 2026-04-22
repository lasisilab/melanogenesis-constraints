"""
fetch_kegg_pathways.py

Fetches the number of KEGG pathways each gene in the 129-gene melanogenesis
network belongs to.

Strategy
--------
1. POST https://mygene.info/v3/querymany
       → batch symbol → NCBI Entrez ID lookup (reliable, designed for this)
2. GET  https://rest.kegg.jp/link/pathway/hsa
       → KEGG ID (hsa:ENTREZID) → pathway set
   KEGG IDs are identical to NCBI gene IDs, so hsa:{ENTREZID} is the key.

Output:  data/kegg_pathway_counts.csv — columns: gene, kegg_id, kegg_pathway_count
"""

import os
import requests
import pandas as pd
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
OUT_CSV     = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_counts.csv')

# ── Load gene list ─────────────────────────────────────────────────────────
df    = pd.read_csv(MASTER_CSV)
genes = df['gene'].str.upper().tolist()
print(f"Loaded {len(genes)} genes from master dataset")

# ── Step 1: gene symbol → NCBI Entrez ID via mygene.info ──────────────────
print("\nStep 1: Looking up NCBI Entrez IDs via mygene.info (batch query)...")
resp = requests.post(
    'http://mygene.info/v3/query',
    json={
        'q':       genes,
        'scopes':  'symbol',
        'species': 'human',
        'fields':  'entrezgene',
    },
    headers={'Content-Type': 'application/json'},
    timeout=60,
)
resp.raise_for_status()

symbol_to_entrez: dict[str, str] = {}
not_found: list[str] = []
for hit in resp.json():
    sym = hit.get('query', '').upper()
    if hit.get('notfound'):
        not_found.append(sym)
    elif 'entrezgene' in hit:
        # entrezgene can be a list if there are multiple hits; take first
        eg = hit['entrezgene']
        if isinstance(eg, list):
            eg = eg[0]
        symbol_to_entrez[sym] = str(int(eg))

print(f"  Mapped {len(symbol_to_entrez)}/{len(genes)} genes to NCBI Entrez IDs")
if not_found:
    print(f"  Not found in mygene.info: {not_found}")

# ── Step 2: KEGG gene-pathway links (hsa:ENTREZID → pathway) ──────────────
print("\nStep 2: Fetching gene-pathway links from KEGG (link/pathway/hsa)...")
resp = requests.get("https://rest.kegg.jp/link/pathway/hsa", timeout=120)
resp.raise_for_status()

kegg_to_pathways: dict[str, set] = defaultdict(set)
n_lines = 0
for line in resp.text.strip().split('\n'):
    parts = line.split('\t')
    if len(parts) == 2:
        kegg_to_pathways[parts[0]].add(parts[1])
        n_lines += 1
print(f"  Loaded {n_lines:,} gene-pathway links for {len(kegg_to_pathways):,} unique genes")

# ── Step 3: Count pathways per melanogenesis gene ──────────────────────────
print("\nStep 3: Counting KEGG pathways per gene...")
rows = []
for gene in genes:
    entrez = symbol_to_entrez.get(gene, '')
    kegg_id = f'hsa:{entrez}' if entrez else ''
    n = len(kegg_to_pathways[kegg_id]) if kegg_id else None
    rows.append({'gene': gene, 'kegg_id': kegg_id, 'kegg_pathway_count': n})

out_df = pd.DataFrame(rows)
out_df.to_csv(OUT_CSV, index=False)
print(f"\nSaved → {OUT_CSV}")

# ── Summary ────────────────────────────────────────────────────────────────
valid = out_df.dropna(subset=['kegg_pathway_count'])
zero  = (valid['kegg_pathway_count'] == 0).sum()
print(f"\nn = {len(valid)} genes with KEGG pathway counts "
      f"({zero} with zero pathways, likely not in KEGG signaling maps)")
if len(valid) > 0:
    print(f"Range: {int(valid['kegg_pathway_count'].min())}–"
          f"{int(valid['kegg_pathway_count'].max())} pathways  "
          f"(median {valid['kegg_pathway_count'].median():.0f})")

    print("\nTop 15 by KEGG pathway count:")
    top = valid.sort_values('kegg_pathway_count', ascending=False).head(15)
    for _, row in top.iterrows():
        print(f"  {row['gene']:12s}  {int(row['kegg_pathway_count'])} pathways")

print("\nDone!")
