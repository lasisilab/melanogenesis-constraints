"""
13_make_sgdp_only_sample_lists.py

Builds per-population sample lists for the SGDP-only genome-wide PBS pipeline.
Mirrors 03_make_sample_lists.py but reads sample IDs directly from the SGDP
VCF and assigns populations from the SGDP sample-ID naming convention
(${PANEL}_${POPULATION}-${INDIV_ID}, e.g. "S_Yoruba-1", "B_Papuan-15").

Outputs (in BASE/data/sgdp_only/):
  samples_sgdp_african.txt        — sub-Saharan continental African
  samples_sgdp_melanesian.txt     — Papuan + Bougainville (same 17 as before)
  samples_sgdp_eastasian.txt      — East Asian
  samples_sgdp_southasian.txt     — South Asian
  samples_sgdp_european.txt       — European
  samples_sgdp_all_keep.txt       — union of the above
  sgdp_population_assignments.tsv — full sample → population mapping for audit

Usage:
    module load Bioinformatics bcftools
    python analysis/cluster/13_make_sgdp_only_sample_lists.py \\
        --vcf  /nfs/turbo/lsa-tlasisi1/tlasisi/reference-genomes/sgdp/<file>.vcf.gz \\
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints

If --vcf is omitted, the script searches the SGDP directory for likely
multi-sample VCFs and prints what it finds so you can pick.

Population assignments follow Mallick et al. 2016 (Nature) region groupings.
Native American, Central Asian/Siberian, Oceanian non-Papuan, and admixed
populations are excluded from the 5-pop PBS panel (kept visible in the audit
TSV for transparency).
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


# ─── SGDP populations → super-population for the 5-pop PBS panel ────────────
# All names follow the Reich Lab / SGDP sample-ID convention.
# Source: Mallick et al. 2016, Nature 538:201 (Table S1).
AFRICAN = {
    "Yoruba", "Mandenka", "Mbuti", "BantuKenya", "BantuSouthAfrica",
    "BantuTutsi", "BantuHerero", "Biaka", "Dinka", "Esan", "Gambian",
    "Ju_hoan_North", "Khomani_San", "Luhya", "Luo", "Mende", "Masai",
    "Somali",
    # Note: Mozabite / Saharawi / Tunisian are North African / admixed
    # with non-African ancestry — excluded here to match the AFRICAN_1KGP +
    # AFRICAN_HGDP scheme in 03_make_sample_lists.py.
}
MELANESIAN = {"Papuan", "Bougainville"}
EAST_ASIAN = {
    "Han", "Japanese", "Korean", "Mongolian", "Naxi", "Daur", "Hezhen",
    "Lahu", "Miao", "She", "Tu", "Tujia", "Xibo", "Yi",
    "Ami", "Atayal", "Cambodian", "Dai", "Mlabri", "Thai",
}
SOUTH_ASIAN = {
    "Bengali", "Brahmin", "Burusho", "Hazara", "Iyer", "Iyengar",
    "Kalash", "Kapu", "Khonda_Dora", "Kshatriya", "Kusunda",
    "Madiga", "Mala", "Pathan", "Punjabi", "Relli", "Sindhi",
    "Sherpa", "Tamil", "Vishwabrahmin", "Balochi", "Brahui", "Makrani",
}
EUROPEAN = {
    "French", "Italian_North", "Sardinian", "Tuscan", "Bergamo",
    "Basque", "English", "Estonian", "Finnish", "Greek", "Hungarian",
    "Icelandic", "Norwegian", "Orcadian", "Polish", "Russian",
    "Spanish", "Saami", "Czech", "Albanian", "Bulgarian",
}

POP_MAP = {}
for s, label in [
    (AFRICAN,    "african"),
    (MELANESIAN, "melanesian"),
    (EAST_ASIAN, "eastasian"),
    (SOUTH_ASIAN,"southasian"),
    (EUROPEAN,   "european"),
]:
    for pop in s:
        POP_MAP[pop] = label


def parse_population_from_sample_id(sample_id):
    """SGDP sample IDs are formatted "{PANEL}_{POPULATION}-{IDX}", e.g.
    "S_Yoruba-1". Returns the population name, or None if the format does
    not match."""
    if "_" not in sample_id or "-" not in sample_id:
        return None
    after_panel = sample_id.split("_", 1)[1]   # drop "S_" / "B_" / "A_" prefix
    pop = after_panel.rsplit("-", 1)[0]
    return pop


def get_vcf_samples(vcf_path):
    out = subprocess.run(
        ["bcftools", "query", "-l", vcf_path],
        capture_output=True, text=True, check=True,
    )
    return [s.strip() for s in out.stdout.splitlines() if s.strip()]


def find_likely_vcfs(directory):
    """List multi-sample VCFs in the SGDP reference dir to help the user
    pick the right one."""
    if not os.path.isdir(directory):
        return []
    cands = []
    for fname in sorted(os.listdir(directory)):
        if fname.endswith((".vcf.gz", ".vcf.bgz", ".bcf")):
            cands.append(os.path.join(directory, fname))
    return cands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf",
                    help="Path to an SGDP VCF (any single chromosome works "
                         "— we only need the sample header)")
    ap.add_argument("--sgdp-dir",
                    default="/nfs/turbo/lsa-tlasisi1/tlasisi/reference-genomes/sgdp",
                    help="Directory containing the SGDP VCFs (used only if "
                         "--vcf is omitted, to help locate the file)")
    ap.add_argument("--base",
                    default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints",
                    help="Cluster working directory (sample lists go here)")
    args = ap.parse_args()

    if not args.vcf:
        cands = find_likely_vcfs(args.sgdp_dir)
        if not cands:
            sys.exit(f"ERROR: no VCFs found in {args.sgdp_dir}. "
                     "Pass --vcf with an explicit path.")
        print(f"Found {len(cands)} candidate VCF(s) in {args.sgdp_dir}:")
        for p in cands:
            print(f"  {p}")
        sys.exit("Re-run with  --vcf <path>  to pick one.")

    if not os.path.exists(args.vcf):
        sys.exit(f"ERROR: VCF not found: {args.vcf}")

    out_dir = os.path.join(args.base, "data", "sgdp_only")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Reading sample IDs from {args.vcf} ...")
    samples = get_vcf_samples(args.vcf)
    print(f"  {len(samples)} samples total in VCF")

    # Assign each sample
    by_pop = defaultdict(list)
    audit_rows = []
    unassigned = []
    for s in samples:
        pop_name = parse_population_from_sample_id(s)
        super_pop = POP_MAP.get(pop_name)
        if super_pop is None:
            unassigned.append(s)
            audit_rows.append((s, pop_name or "", "EXCLUDED"))
        else:
            by_pop[super_pop].append(s)
            audit_rows.append((s, pop_name, super_pop))

    print("\nPer-super-population counts (5-pop panel):")
    for sp in ("african", "melanesian", "eastasian", "southasian", "european"):
        print(f"  {sp:11s}: n = {len(by_pop[sp]):>3d}")
    print(f"  excluded   : n = {len(unassigned):>3d}  "
          "(Native American, Siberian, North African, etc.)")

    # Write per-pop lists + union
    def write_list(name, items):
        path = os.path.join(out_dir, f"samples_sgdp_{name}.txt")
        with open(path, "w") as f:
            f.write("\n".join(items) + ("\n" if items else ""))
        print(f"  → {path}  ({len(items)} samples)")

    print("\nWriting sample lists...")
    for sp in ("african", "melanesian", "eastasian", "southasian", "european"):
        write_list(sp, by_pop[sp])
    union = []
    for sp in ("african", "melanesian", "eastasian", "southasian", "european"):
        union.extend(by_pop[sp])
    write_list("all_keep", union)

    # Audit TSV
    audit_path = os.path.join(out_dir, "sgdp_population_assignments.tsv")
    with open(audit_path, "w") as f:
        f.write("sample_id\tsgdp_population\tsuper_population\n")
        for row in audit_rows:
            f.write("\t".join(row) + "\n")
    print(f"  → {audit_path}  ({len(audit_rows)} rows)")

    print("\nDone.")


if __name__ == "__main__":
    main()
