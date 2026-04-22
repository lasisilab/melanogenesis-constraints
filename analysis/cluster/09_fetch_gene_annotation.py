"""
09_fetch_gene_annotation.py

Downloads GENCODE v38 (hg38) GTF and extracts all protein-coding genes
as a BED file with ±10 kb flanking windows.  Used as the gene list for
genome-wide PBS computation in 10_compute_pbs_genomewide.py.

Output:
    data/gencode_v38_protein_coding_10kb.bed
        columns: chrom  start  end  gene_name  gene_id  strand

Run once locally or on the cluster before submitting 10_compute_pbs_genomewide.sh.

Usage:
    python analysis/cluster/09_fetch_gene_annotation.py \
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints
"""

import argparse
import gzip
import os
import urllib.request

GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_38/gencode.v38.annotation.gtf.gz"
)

FLANK = 10_000   # bp added to each side

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base",
                   default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints")
    return p.parse_args()

def main():
    args = parse_args()
    data_dir = os.path.join(args.base, "data")
    os.makedirs(data_dir, exist_ok=True)

    gtf_gz  = os.path.join(data_dir, "gencode.v38.annotation.gtf.gz")
    bed_out = os.path.join(data_dir, "gencode_v38_protein_coding_10kb.bed")

    # ── Download GTF ─────────────────────────────────────────────────────────
    if not os.path.exists(gtf_gz):
        print(f"Downloading GENCODE v38 GTF (~1.5 GB)...")
        urllib.request.urlretrieve(GENCODE_URL, gtf_gz)
        print(f"  Saved → {gtf_gz}")
    else:
        print(f"GTF already present → {gtf_gz}")

    # ── Parse protein-coding genes ────────────────────────────────────────────
    print("Parsing protein-coding genes...")
    genes = {}   # gene_id → (chrom, start, end, name, strand)

    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue

            attrs = {}
            for item in fields[8].rstrip(";").split(";"):
                item = item.strip()
                if " " in item:
                    k, v = item.split(" ", 1)
                    attrs[k] = v.strip('"')

            if attrs.get("gene_type") != "protein_coding":
                continue

            chrom  = fields[0]
            start  = int(fields[3]) - 1   # 0-based BED
            end    = int(fields[4])
            strand = fields[6]
            gene_id   = attrs.get("gene_id",   "").split(".")[0]
            gene_name = attrs.get("gene_name", gene_id)

            # Keep only standard chromosomes (chr1–22, chrX, chrY)
            if not (chrom.startswith("chr") and
                    chrom[3:].lstrip("0") in
                    [str(i) for i in range(1, 23)] + ["X", "Y"]):
                continue

            # If gene appears multiple times (shouldn't for gene lines), keep widest
            if gene_id not in genes:
                genes[gene_id] = (chrom, start, end, gene_name, strand)
            else:
                prev = genes[gene_id]
                genes[gene_id] = (chrom, min(prev[1], start),
                                         max(prev[2], end),
                                  gene_name, strand)

    print(f"  {len(genes):,} protein-coding genes parsed")

    # ── Write BED with flanking windows ──────────────────────────────────────
    # Sort by chrom (numeric), then start
    def sort_key(entry):
        chrom = entry[0]
        c = chrom[3:]
        try:
            return (int(c), entry[1])
        except ValueError:
            return (99, entry[1])

    rows = list(genes.values())
    rows.sort(key=sort_key)

    with open(bed_out, "w") as fh:
        fh.write("# chrom\tstart\tend\tgene_name\tgene_id\tstrand\n")
        for gene_id, (chrom, start, end, name, strand) in genes.items():
            s = max(0, start - FLANK)
            e = end + FLANK
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t{gene_id}\t{strand}\n")

    print(f"  Written → {bed_out}  ({len(genes):,} genes, ±{FLANK//1000} kb flanks)")

if __name__ == "__main__":
    main()
