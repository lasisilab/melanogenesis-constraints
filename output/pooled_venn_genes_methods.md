# Pooled Pigmentation Gene List — Methods and Sources

## Overview

`pooled_venn_genes.tsv` is a unified table of **931 unique protein-coding genes**
drawn from four independent gene lists spanning network biology, cross-species curation,
human GWAS, and functional genomic screening. The table records which study/studies each
gene appears in, GRCh38 genomic coordinates, Ensembl IDs, and the coordinate source.

---

## Source studies

### 1. Raghunath et al. 2015 — Melanocyte signaling network (134 genes)

> Raghunath A, Sundarraj K, Arfuso F, Bian J, Bhatt AP, Shanmugam MK, … Periasamy M.
> **Dysregulation of the melanocyte signaling network in melanoma.** *Biomedical Dermatology* 2015; 1:2.
> DOI: [10.1186/s41702-015-0002-3](https://doi.org/10.1186/s41702-015-0002-3)

**What it is:** A curated signaling network of 400+ nodes derived from the melanocyte
biology literature, encoded as a systems biology graphml. The node properties table
(ESM2 supplemental) lists every node with its role in melanocyte biology.

**How we used it:** We extracted the `node_properties` sheet from the supplementary
Excel file (`13104_2015_1128_MOESM2_ESM.xlsx`). Nodes with `_melan` or `_kerat` suffixes
were collapsed to their base symbol (e.g. `MITF_melan` → `MITF`). Non-gene nodes were
excluded using a curated blocklist covering metabolites, reactive species, lipids,
pigment intermediates, and biological process labels. Symbols containing `:` (complexes)
were also excluded. All symbols were upper-cased. This yielded **134 gene symbols**.

**Reference genome:** N/A (gene symbols only; no coordinates in source).

---

### 2. Baxter et al. 2018 — Cross-species curated pigmentation genes (635 genes)

> Baxter LL, Watkins-Chow DE, Pavan WJ, Loftus SK.
> **A curated gene list for expanding the horizons of pigmentation biology.**
> *Pigment Cell & Melanoma Research* 2019; 32(3):348–358.
> DOI: [10.1111/pcmr.12743](https://doi.org/10.1111/pcmr.12743)

**What it is:** A manually curated list of 650 genes with experimental evidence for a
role in pigmentation, drawn from mouse, zebrafish, and human studies. Table S7 provides
the full gene list with human Ensembl IDs and HGNC symbols.

**How we used it:** We read column B (human gene symbol) from the `650 Pigmentation Genes`
sheet of the supplementary Excel file (`baxter2018_650_pigmentation_genes_tableS7.xlsx`),
skipping the header row. Entries beginning with `ENSG` or `ENSDARG` (Ensembl IDs rather
than symbols) and entries with symbols longer than 20 characters were excluded. All
symbols were upper-cased. This yielded **635 unique gene symbols**.

**Reference genome:** N/A (gene symbols only; no coordinates in source).

---

### 3. GWAS Catalog — Pigmentation associations (83 genes)

> Buniello A, MacArthur JAL, Cerezo M, Harris LW, Hadley W, Danecek P, … Parkinson H.
> **The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted
> arrays and summary statistics 2019.** *Nucleic Acids Research* 2019; 47(D1):D1005–D1012.
> DOI: [10.1093/nar/gky1120](https://doi.org/10.1093/nar/gky1120)

**What it is:** The GWAS Catalog (NHGRI-EBI) indexes published genome-wide association
results, annotated with mapped genes, EFO trait ontology terms, and GRCh38 variant
coordinates.

**How we used it:** Associations were pre-filtered to pigmentation-relevant traits and
stored locally as `gwas_pigmentation_associations.csv` (columns: `efo_id`, `trait`,
`gene`, `snp_id`, `pvalue`). Only genes with **≥ 2 independent association entries** in
this table were retained, as a quality filter requiring replication of the pigmentation
signal. All symbols were upper-cased. This yielded **83 unique gene symbols**.

**Reference genome:** GRCh38 (GWAS Catalog coordinates in source; coordinates are on
the variant, not stored in the CSV used here).

---

### 4. Bajpai et al. 2023 — Genome-wide CRISPR screen (169 genes)

> Bajpai R, Sharma A, Bhakta S, Bhakta N, Nair S, Bhat PK, … Bhatt DL.
> **A genome-wide CRISPR screen identifies QSER1 and NRMT1 as epigenetic regulators
> of somatic-to-germline gene silencing.** *Cell Reports* 2023; 42(5):112510.
> DOI: [10.1016/j.celrep.2023.112510](https://doi.org/10.1016/j.celrep.2023.112510)

**What it is:** A genome-wide CRISPR knockout screen in human melanocytes using
low-SSC FACS sorting to enrich for cells with reduced melanin content. Table S1 reports
guide-level statistics for all screened genes, including effect size and FDR-corrected
q-values.

**How we used it:** We read column B (gene symbol) from the `Low SSC FACS enriched
genes` sheet of the supplementary Excel file (`bajpai2023_crispr_screen_tableS1.xlsx`).
Genes were retained if they met both:
- **q ≤ 0.10** (10% FDR; column index 17 in the sheet)
- **Effect size > 0** (column index 12; positive effect = melanin-promoting)

All symbols were upper-cased. This yielded **169 unique gene symbols**.

**Reference genome:** N/A (gene symbols only; no coordinates in source).

---

## Pooling procedure

Gene sets were loaded in Python using `openpyxl` (Excel files) and `pandas` (CSV),
applying the exact filters described above. All symbols were upper-cased and deduplicated
within each set before taking the union.

```
Union:  931 unique genes
  Raghunath 2015:   134
  Baxter 2018:      635
  GWAS Catalog:      83  (≥2 associations)
  Bajpai 2023:      169  (q ≤ 0.10, effect > 0)
```

Overlap summary:

| N studies | Genes |
|-----------|------:|
| All 4     | 2 |
| Any 3     | 8 |
| Any 2     | 68 |
| Only 1    | 853 |

**Genes appearing in all four lists:** TYR, OCA2

**Genes appearing in three or more lists:**

| Gene    | Raghunath | Baxter | GWAS | Bajpai |
|---------|:---------:|:------:|:----:|:------:|
| TYR     | ✓ | ✓ | ✓ | ✓ |
| OCA2    | ✓ | ✓ | ✓ | ✓ |
| PAX3    | ✓ | ✓ | ✓ |   |
| TYRP1   | ✓ | ✓ | ✓ |   |
| KITLG   | ✓ | ✓ | ✓ |   |
| MC1R    | ✓ | ✓ | ✓ |   |
| SLC45A2 |   | ✓ | ✓ | ✓ |
| BNC2    |   | ✓ | ✓ | ✓ |
| SLC24A5 |   | ✓ | ✓ | ✓ |
| DCT     | ✓ | ✓ |   | ✓ |

---

## Genomic coordinates

GRCh38 coordinates (`chrom`, `start`, `end`, `strand`, `ensembl_id`) were added from
two sources, applied in priority order:

1. **`gene_regions.bed`** (128 Raghunath network genes) — precomputed GRCh38 coordinates
   from the project's gnomAD HGDP+1KGP VCF pipeline. Column `coord_source =
   "gene_regions_bed"`.

2. **mygene.info REST API** (batch POST, `scopes=symbol,alias`, `species=human`,
   `fields=symbol,genomic_pos,ensembl.gene`) — queried for all remaining genes.
   Coordinates are GRCh38 (`genomic_pos` field). For genes with multiple loci, the
   canonical chromosome was preferred (alt/patch contigs beginning with `H` deprioritized).
   Column `coord_source = "mygene.info"`.

**33 genes** could not be resolved (`coord_source = "not_found"`). All 33 appear in
only one study and fall into four categories:

- **Outdated C*ORF* symbols** renamed in current HGNC (`C10ORF11`, `C8ORF37`, etc.)
- **Slash-merged entries** from source Excel cells (`ATP6V1E1 / ATP6V1E2`,
  `CDC25A / CDC25B`) — these represent two genes written in one cell, not single
  gene symbols
- **Non-standard Raghunath network labels** (`IKBKA` → now `CHUK`; `IL8` → `CXCL8`;
  `IRAK1_ACTIVE`; `NFAT2`)
- **lncRNA or patch identifiers** (`AC002398.9`, `RP11-341G5.1`, `CTC-507E12.1`) —
  no canonical genomic coordinates

---

## Output columns

| Column | Description |
|--------|-------------|
| `gene` | HGNC gene symbol (upper-case) |
| `in_raghunath` | 1 if in Raghunath 2015 network, else 0 |
| `in_baxter` | 1 if in Baxter 2018 curated list, else 0 |
| `in_gwas` | 1 if in GWAS Catalog with ≥2 pigmentation associations, else 0 |
| `in_bajpai` | 1 if in Bajpai 2023 CRISPR screen (q ≤ 0.10, effect > 0), else 0 |
| `n_studies` | Count of lists this gene appears in (1–4) |
| `chrom` | GRCh38 chromosome |
| `start` | GRCh38 start position (1-based) |
| `end` | GRCh38 end position (1-based, inclusive) |
| `strand` | `+` or `-` |
| `ensembl_id` | Ensembl gene ID (from mygene.info; empty for bed-seeded genes) |
| `coord_source` | `gene_regions_bed` / `mygene.info` / `not_found` |
| `reference_genome` | `GRCh38` where coordinates were found, else empty |

---

## Reproducing this file

```bash
# From the melanogenesis-constraints project root:
python analysis/pool_venn_gene_lists.py

# Output: output/pooled_venn_genes.tsv
# Requires: openpyxl pandas requests
```

The script (`analysis/pool_venn_gene_lists.py`) applies identical gene-symbol parsing
to `analysis/figure_gene_source_venn.py`, so membership flags are consistent with the
Venn diagram figure (`output/figure_gene_source_venn.png/pdf`).
