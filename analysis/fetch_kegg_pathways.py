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

Output:
  data/kegg_pathway_counts.csv — gene, kegg_id, kegg_pathway_count
  data/kegg_pathway_lists.csv  — gene, kegg_id, pathway_id, pathway_name (long form)
"""

import os
import requests
import pandas as pd
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER_CSV  = os.path.join(PROJECT_DIR, 'data', 'network_constraint_gtex.csv')
OUT_CSV     = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_counts.csv')
LISTS_CSV   = os.path.join(PROJECT_DIR, 'data', 'kegg_pathway_lists.csv')

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

# Pathway ID → name (single REST call, ~340 entries)
print("\nStep 2b: Fetching pathway names (list/pathway/hsa)...")
resp = requests.get("https://rest.kegg.jp/list/pathway/hsa", timeout=60)
resp.raise_for_status()
pathway_names: dict[str, str] = {}
for line in resp.text.strip().split('\n'):
    parts = line.split('\t')
    if len(parts) == 2:
        # Normalize ID to "path:hsa04010" form (matches link output)
        pid = parts[0] if parts[0].startswith('path:') else f'path:{parts[0]}'
        pathway_names[pid] = parts[1]
print(f"  Loaded {len(pathway_names):,} pathway names")

# ── Step 3: Count pathways and emit long-form list per gene ───────────────
print("\nStep 3: Counting KEGG pathways and building long-form pathway list...")
rows_count = []
rows_long  = []
for gene in genes:
    entrez = symbol_to_entrez.get(gene, '')
    kegg_id = f'hsa:{entrez}' if entrez else ''
    pathways = kegg_to_pathways.get(kegg_id, set()) if kegg_id else set()
    rows_count.append({'gene': gene, 'kegg_id': kegg_id,
                       'kegg_pathway_count': len(pathways) if kegg_id else None})
    for pid in pathways:
        rows_long.append({'gene': gene, 'kegg_id': kegg_id,
                          'pathway_id': pid,
                          'pathway_name': pathway_names.get(pid, '')})

out_df  = pd.DataFrame(rows_count)
long_df = pd.DataFrame(rows_long).sort_values(['gene', 'pathway_id'])
out_df.to_csv(OUT_CSV, index=False)
long_df.to_csv(LISTS_CSV, index=False)
print(f"\nSaved counts → {OUT_CSV}")
print(f"Saved per-pathway list → {LISTS_CSV}  ({len(long_df)} rows)")

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
