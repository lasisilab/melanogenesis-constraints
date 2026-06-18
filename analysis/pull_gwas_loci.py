"""
pull_gwas_loci.py

Pull GWAS Catalog loci directly from the EBI REST API and write a tidy table.

Query by any combination of:
  --genes   gene symbols (comma-separated or a file, one per line)  → mapped-gene loci
  --region  chrom:start-end (GRCh38)                                 → loci in an interval
  --trait   EFO trait text (e.g. "skin pigmentation")               → trait associations

Output columns:
  rsid, chrom, pos, mapped_genes, risk_allele, risk_freq,
  p_value, or_beta, ci_text, reported_trait, efo_traits,
  pubmed_id, study_accession

Coordinates are GRCh38 (matches the project VCFs).

Examples:
  python analysis/pull_gwas_loci.py --genes MC1R,OCA2,SLC24A5 --out gwas_pigment_genes.tsv
  python analysis/pull_gwas_loci.py --genes data/network_genes.txt --out gwas_network.tsv
  python analysis/pull_gwas_loci.py --region 15:27719008-28099342 --out gwas_herc2_oca2.tsv
  python analysis/pull_gwas_loci.py --trait "vitamin D measurement" --out gwas_vitd.tsv

Docs: https://www.ebi.ac.uk/gwas/rest/docs/api
"""

import argparse
import os
import sys
import time

import pandas as pd
import requests

BASE = "https://www.ebi.ac.uk/gwas/rest/api"
PAGE_SIZE = 200
SLEEP = 0.1          # be polite to the API
TIMEOUT = 60
SESSION = requests.Session()
SESSION.headers.update({"Accept": "application/json"})


# ── HTTP helpers ─────────────────────────────────────────────────────────────
def get_json(url, params=None, retries=4):
    """GET with simple exponential backoff. Returns parsed JSON or None on 404."""
    for attempt in range(retries):
        try:
            r = SESSION.get(url, params=params, timeout=TIMEOUT)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            time.sleep(SLEEP)
            return r.json()
        except requests.RequestException as e:
            if attempt == retries - 1:
                print(f"  ! request failed: {url} ({e})", file=sys.stderr)
                return None
            time.sleep(2 ** attempt)
    return None


def iter_embedded(url, key, params=None):
    """Yield items from a paginated HAL collection, following _links.next."""
    params = dict(params or {})
    params.setdefault("size", PAGE_SIZE)
    while url:
        data = get_json(url, params=params)
        params = None  # next links already carry query params
        if not data:
            return
        for item in data.get("_embedded", {}).get(key, []):
            yield item
        url = data.get("_links", {}).get("next", {}).get("href")


def follow(obj, rel):
    """Return JSON for a HAL link relation on an object, or None."""
    href = obj.get("_links", {}).get(rel, {}).get("href")
    if not href:
        return None
    return get_json(href.split("{")[0])  # strip templated query suffix


# ── Flattening ───────────────────────────────────────────────────────────────
def snp_location(snp):
    """(chrom, pos, mapped_genes) from a SNP object."""
    locs = snp.get("locations") or [{}]
    chrom = locs[0].get("chromosomeName")
    pos = locs[0].get("chromosomePosition")
    genes = sorted({
        gc.get("gene", {}).get("geneName")
        for gc in (snp.get("genomicContexts") or [])
        if gc.get("gene", {}).get("geneName")
    })
    return chrom, pos, ";".join(genes)


def association_rows(assoc, snp_cache):
    """Flatten one association into one row per linked SNP."""
    pmant = assoc.get("pvalueMantissa")
    pexp = assoc.get("pvalueExponent")
    p_value = (assoc.get("pvalue")
               or (f"{pmant}e{pexp}" if pmant is not None and pexp is not None else None))

    orval = assoc.get("orPerCopyNum")
    beta = assoc.get("betaNum")
    if orval:
        or_beta = f"OR={orval}"
    elif beta is not None:
        or_beta = f"beta={beta} {assoc.get('betaUnit') or ''} {assoc.get('betaDirection') or ''}".strip()
    else:
        or_beta = None

    risk_alleles = []
    reported_genes = set()
    for locus in assoc.get("loci") or []:
        for ra in locus.get("strongestRiskAlleles") or []:
            if ra.get("riskAlleleName"):
                risk_alleles.append(ra["riskAlleleName"])
        for g in locus.get("authorReportedGenes") or []:
            if g.get("geneName"):
                reported_genes.add(g["geneName"])

    efo = follow(assoc, "efoTraits")
    efo_traits = ";".join(
        t.get("trait", "") for t in (efo or {}).get("_embedded", {}).get("efoTraits", [])
    ) if efo else ""

    study = follow(assoc, "study")
    pubmed = accession = None
    if study:
        pubmed = (study.get("publicationInfo") or {}).get("pubmedId")
        accession = study.get("accessionId")

    # SNPs carrying this association
    snps = follow(assoc, "snps")
    snp_list = (snps or {}).get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
    if not snp_list:
        snp_list = [{}]

    rows = []
    for snp in snp_list:
        rsid = snp.get("rsId")
        if rsid and rsid in snp_cache:
            chrom, pos, mapped = snp_cache[rsid]
        else:
            chrom, pos, mapped = snp_location(snp)
            if rsid:
                snp_cache[rsid] = (chrom, pos, mapped)
        rows.append({
            "rsid": rsid,
            "chrom": chrom,
            "pos": pos,
            "mapped_genes": mapped or ";".join(sorted(reported_genes)),
            "risk_allele": ";".join(risk_alleles) or None,
            "risk_freq": assoc.get("riskFrequency"),
            "p_value": p_value,
            "or_beta": or_beta,
            "ci_text": assoc.get("range"),
            "reported_trait": ";".join(sorted(reported_genes)) or None,
            "efo_traits": efo_traits,
            "pubmed_id": pubmed,
            "study_accession": accession,
        })
    return rows


# ── Query modes ──────────────────────────────────────────────────────────────
def loci_by_gene(gene, snp_cache):
    print(f"[gene] {gene}")
    rows = []
    url = f"{BASE}/singleNucleotidePolymorphisms/search/findByGene"
    for snp in iter_embedded(url, "singleNucleotidePolymorphisms", params={"geneName": gene}):
        rsid = snp.get("rsId")
        snp_cache.setdefault(rsid, snp_location(snp)) if rsid else None
        assocs = follow(snp, "associations")
        for assoc in (assocs or {}).get("_embedded", {}).get("associations", []):
            rows.extend(association_rows(assoc, snp_cache))
    return rows


def loci_by_region(chrom, start, end, snp_cache):
    print(f"[region] {chrom}:{start}-{end}")
    rows = []
    url = f"{BASE}/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"
    params = {"chrom": str(chrom), "bpStart": int(start), "bpEnd": int(end)}
    for snp in iter_embedded(url, "singleNucleotidePolymorphisms", params=params):
        rsid = snp.get("rsId")
        snp_cache.setdefault(rsid, snp_location(snp)) if rsid else None
        assocs = follow(snp, "associations")
        for assoc in (assocs or {}).get("_embedded", {}).get("associations", []):
            rows.extend(association_rows(assoc, snp_cache))
    return rows


def loci_by_trait(trait, snp_cache):
    print(f"[trait] {trait}")
    rows = []
    found = get_json(f"{BASE}/efoTraits/search/findByEfoTrait", params={"trait": trait})
    traits = (found or {}).get("_embedded", {}).get("efoTraits", [])
    if not traits:
        print(f"  no EFO trait matched '{trait}'", file=sys.stderr)
        return rows
    for t in traits:
        efo_id = t.get("shortForm")
        print(f"  matched {t.get('trait')} ({efo_id})")
        url = f"{BASE}/efoTraits/{efo_id}/associations"
        for assoc in iter_embedded(url, "associations"):
            rows.extend(association_rows(assoc, snp_cache))
    return rows


# ── Main ─────────────────────────────────────────────────────────────────────
def read_genes(arg):
    if os.path.isfile(arg):
        with open(arg) as f:
            return [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    return [g.strip() for g in arg.split(",") if g.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", help="comma-separated symbols or a file (one per line)")
    ap.add_argument("--region", help="chrom:start-end (GRCh38)")
    ap.add_argument("--trait", help="EFO trait text, e.g. 'skin pigmentation'")
    ap.add_argument("--out", required=True, help="output TSV path")
    args = ap.parse_args()

    if not (args.genes or args.region or args.trait):
        ap.error("provide at least one of --genes / --region / --trait")

    snp_cache = {}
    rows = []
    if args.genes:
        for g in read_genes(args.genes):
            rows.extend(loci_by_gene(g, snp_cache))
    if args.region:
        chrom, span = args.region.split(":")
        start, end = span.split("-")
        rows.extend(loci_by_region(chrom, start, end, snp_cache))
    if args.trait:
        rows.extend(loci_by_trait(args.trait, snp_cache))

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset=["rsid", "p_value", "efo_traits", "study_accession"])
        df = df.sort_values(["chrom", "pos"], na_position="last")
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {len(df)} loci → {args.out}")


if __name__ == "__main__":
    main()
