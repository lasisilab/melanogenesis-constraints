# Melanogenesis Network Constraint Analysis

Investigating the relationship between network topology (centrality) and evolutionary
constraint (LOEUF) in melanogenesis pathway genes.

## data/

| File | Source | Description |
|------|--------|-------------|
| `13104_2015_1128_MOESM2_ESM.xlsx` | Raghunath et al. 2015 (PMID: 26286471) | Supplementary data from melanocyte signaling network paper. Contains `node_properties` sheet with betweenness centrality, in/out-degree for 265 network nodes (208 unique after collapsing `_melan`/`_kerat` suffixes; ~134 are actual genes). |
| `baxter2018_650_pigmentation_genes_tableS7.xlsx` | Baxter et al. 2018 (PMID: 30098108, DOI: 10.1111/pcmr.12743) | Table S7: curated list of 650 cross-species pigmentation genes from OMIM, MGI, ZFIN, GO, and PubMed. Contains 635 unique human gene symbols. Columns include Ensembl ID, human/mouse/zebrafish symbols, phenotype location, source databases, and species with phenotype. |
| `bajpai2023_crispr_screen_tableS1.xlsx` | Bajpai et al. 2023 (PMID: 37535747, DOI: 10.1126/science.ade6289) | Table S1: genome-wide CRISPR screen results in MNT-1 melanoma cells. Contains CasTLE effect sizes, confidence scores, and q-values for ~4,950 genes. 169 melanin-promoting hits at 10% FDR. |
| `gnomad_constraint.txt.bgz` | gnomAD v2.1.1 | Gene-level constraint metrics including LOEUF (loss-of-function observed/expected upper bound fraction) and pLI. |
| `gwas_catalog_pigmentation.csv` | NHGRI-EBI GWAS Catalog | Pigmentation-related GWAS associations downloaded from the catalog. |
| `gwas_gene_summary.csv` | Derived | Summary of GWAS genes mapped from pigmentation associations. |
| `gwas_pigmentation_associations.csv` | Derived | Processed GWAS pigmentation associations. |
| `gwas_pigmentation_with_functional_class.csv` | Derived | GWAS associations annotated with functional class. |
| `pigmentation_gwas_constraint.csv` | Derived | GWAS pigmentation genes merged with gnomAD constraint data. |

## analysis/

| File | Description |
|------|-------------|
| `melanogenesis_network_constraint_v2.ipynb` | Main analysis notebook. Merges Raghunath network with gnomAD constraint. Creates 3-panel figure (centrality scatter, functional category boxplot, disease class boxplot). Defines 6 manually curated functional categories. |
| `gene_list_overlap_analysis.py` | Overlap analysis across Raghunath, Baxter, and Bajpai gene lists. Computes pairwise and three-way overlaps, generates Venn diagram, assesses grant network size claim. |

## output/

| File | Description |
|------|-------------|
| `network_constraint_data.csv` | 130 genes with columns: gene, BetweennessCentrality, LOEUF, pLI, functional_category, disease_class. Core data for all figures. |
| `supplementary_table_genes.csv` | 87 GWAS-associated pigmentation genes. |
| `gene_list_overlap_results.json` | Full overlap analysis results with all gene lists and summary counts. |
| `gene_list_union_annotated.csv` | Union of all 873 genes annotated by source (Raghunath/Baxter/Bajpai). |
| `gene_overlap_venn.png` / `.pdf` | Venn diagram and bar chart of gene overlap across three sources. |

## Key findings

- Raghunath network contains 134 unique gene nodes (plus 74 non-gene nodes)
- Only 3 genes appear in all three sources: **TYR, DCT, OCA2**
- Union of all three sources: **873 unique genes**
- Complete LOEUF separation between syndromic (max 0.658) and isolated pigment (min 0.856) disease genes (p = 6.78e-05)
- Topologically central genes (high betweenness) are generic signaling hubs, not melanogenesis-specific effectors
