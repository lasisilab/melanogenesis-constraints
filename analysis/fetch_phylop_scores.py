"""
fetch_phylop_scores.py

Fetches mean PhyloP 100-way vertebrate conservation scores for each gene
in the Raghunath melanogenesis network.

Steps:
  1. Extract gene list from LOEUF_by_functional_category.xlsx (the 129 network genes)
  2. Query Ensembl REST API for hg38 genomic coordinates per gene symbol
  3. Query UCSC REST API for mean phyloP100way score over each gene interval
  4. Write data/phylop_scores.csv

Output columns:
  gene, ensembl_id, chrom, start, end, mean_phylop_100way
"""

import time
import requests
import pandas as pd

LOEUF_FILE = "../LOEUF_by_functional_category.xlsx"
OUTPUT_FILE = "../data/phylop_scores.csv"
ENSEMBL_REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{}"
UCSC_API = "https://api.genome.ucsc.edu/getData/track"
GENOME = "hg38"
TRACK = "phyloP100way"


# ---------------------------------------------------------------------------
# Step 1: Extract gene list from LOEUF Excel file (the 129 network genes)
# ---------------------------------------------------------------------------
def get_network_genes(path):
    df = pd.read_excel(path, sheet_name="All Genes by Category")
    genes = df["Gene"].dropna().str.upper().unique().tolist()
    print(f"  {len(genes)} unique genes from LOEUF_by_functional_category.xlsx")
    return genes


# ---------------------------------------------------------------------------
# Step 2: Get hg38 coordinates from Ensembl REST API
# ---------------------------------------------------------------------------
def get_coords_for_gene(gene):
    url = ENSEMBL_REST.format(gene)
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=15)
        if r.status_code == 200:
            data = r.json()
            return {
                "gene": gene,
                "ensembl_id": data.get("id"),
                "chrom": "chr" + str(data.get("seq_region_name", "")),
                "start": data.get("start"),
                "end": data.get("end"),
            }
    except Exception:
        pass
    return None


def get_gene_coordinates(genes):
    print(f"  Querying Ensembl REST API for {len(genes)} genes...")
    rows = []
    for i, gene in enumerate(genes):
        result = get_coords_for_gene(gene)
        if result:
            rows.append(result)
        else:
            print(f"    WARNING: No coordinates found for {gene}")
        time.sleep(0.1)  # Ensembl rate limit: ~15 req/s
        if (i + 1) % 30 == 0:
            print(f"    {i + 1}/{len(genes)} queried...")

    df = pd.DataFrame(rows)
    # Keep only standard chromosomes
    df = df[df["chrom"].str.match(r"^chr(\d+|[XY])$")]
    print(f"  Coordinates retrieved for {len(df)} / {len(genes)} genes")
    missing = set(genes) - set(df["gene"])
    if missing:
        print(f"  Missing: {sorted(missing)}")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Step 3: Fetch mean PhyloP score for each gene from UCSC REST API
# ---------------------------------------------------------------------------
def fetch_phylop_for_gene(chrom, start, end, retries=3):
    params = {
        "genome": GENOME,
        "track": TRACK,
        "chrom": chrom,
        "start": start,
        "end": end,
    }
    for attempt in range(retries):
        try:
            r = requests.get(UCSC_API, params=params, timeout=60)
            if r.status_code == 200:
                data = r.json()
                # UCSC returns scores under the track name key, not "data"
                scores = data.get(TRACK, data.get("data", []))
                if scores:
                    values = [
                        s.get("value")
                        for s in scores
                        if s.get("value") is not None
                    ]
                    if values:
                        return sum(values) / len(values)
            return None
        except Exception:
            time.sleep(2 ** attempt)
    return None


def fetch_all_phylop(coords_df):
    print(f"  Fetching PhyloP scores for {len(coords_df)} genes (this may take a few minutes)...")
    phylop_scores = []
    for i, row in coords_df.iterrows():
        score = fetch_phylop_for_gene(row["chrom"], row["start"], row["end"])
        phylop_scores.append(score)
        if (i + 1) % 20 == 0:
            print(f"    {i + 1}/{len(coords_df)} done...")
        time.sleep(0.3)  # be polite to UCSC

    coords_df = coords_df.copy()
    coords_df["mean_phylop_100way"] = phylop_scores
    return coords_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Step 1: Extracting gene list from LOEUF_by_functional_category.xlsx...")
    genes = get_network_genes(LOEUF_FILE)

    print("\nStep 2: Getting Ensembl hg38 coordinates...")
    coords = get_gene_coordinates(genes)

    print("\nStep 3: Fetching PhyloP scores from UCSC...")
    scored = fetch_all_phylop(coords)

    print("\nStep 4: Saving output...")
    scored.to_csv(OUTPUT_FILE, index=False)
    print(f"  Saved {len(scored)} rows to {OUTPUT_FILE}")

    # Sanity check
    print("\nSanity check — PhyloP scores for key genes:")
    key_genes = ["TYR", "TYRP1", "DCT", "SOX10", "PAX3", "MC1R", "MITF"]
    check = scored[scored["gene"].isin(key_genes)][
        ["gene", "mean_phylop_100way"]
    ].sort_values("mean_phylop_100way", ascending=False)
    print(check.to_string(index=False))

    n_missing = scored["mean_phylop_100way"].isna().sum()
    print(f"\n  {len(scored) - n_missing} genes with PhyloP scores, {n_missing} missing")
