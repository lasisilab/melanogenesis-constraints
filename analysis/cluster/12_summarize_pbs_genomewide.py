"""
12_summarize_pbs_genomewide.py

Concatenates per-chromosome genome-wide PBS tables (from 11_compute_pbs_genomewide.py),
computes genome-wide percentile ranks for each scan, and generates a comparison
figure showing where the 129-gene melanogenesis network falls in the genome-wide
distribution.

Run after all 22 chromosome PBS jobs complete:
    python analysis/cluster/12_summarize_pbs_genomewide.py \
        --base /nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints

Outputs:
    output/pbs_genomewide_all.csv           — full genome-wide gene-level PBS table
    output/pbs_genomewide_percentiles.csv   — network genes with genome-wide percentiles
    output/figure_pbs_genomewide_compare.png/pdf
"""

import argparse
import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CATEGORY_ORDER = [
    "Pigment-specific",
    "Developmental/NC",
    "Generic signaling",
    "Cytokines/growth factors",
    "Apoptosis/cell death",
    "Other",
]
CATEGORY_COLORS = {
    "Pigment-specific":          "#D94040",
    "Developmental/NC":          "#E8907E",
    "Generic signaling":         "#F5C242",
    "Cytokines/growth factors":  "#4878CF",
    "Apoptosis/cell death":      "#6BAD6B",
    "Other":                     "#B0B0B0",
}

PBS_SCANS = [
    ("pbs1_african",    "PBS-1\nAfrican / S.Asian out"),
    ("pbs2_african",    "PBS-2\nAfrican / Eur out"),
    ("pbs3_melanesian", "PBS-3\nMelanesian / S.Asian out"),
    ("pbs4_melanesian", "PBS-4\nMelanesian / Eur out"),
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base",
                   default="/nfs/turbo/lsa-tlasisi1/tlasisi/melanogenesis-constraints")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir  = os.path.join(args.base, "output")
    data_dir = os.path.join(args.base, "data")

    # ── Concatenate per-chromosome PBS tables ─────────────────────────────────
    chr_files = sorted(glob.glob(os.path.join(out_dir, "pbs_genomewide_chr*.csv")))
    if not chr_files:
        raise FileNotFoundError(
            "No pbs_genomewide_chr*.csv files found in output/. "
            "Run 11_compute_pbs_genomewide.sh first."
        )
    missing = [c for c in range(1, 23)
               if not os.path.exists(os.path.join(out_dir, f"pbs_genomewide_chr{c}.csv"))]
    if missing:
        print(f"WARNING: Missing chromosomes: {missing}. Proceeding with available data.")

    print(f"Concatenating {len(chr_files)} chromosome files...")
    gw = pd.concat([pd.read_csv(f) for f in chr_files], ignore_index=True)

    for col in ["pbs1_african", "pbs2_african", "pbs3_melanesian", "pbs4_melanesian"]:
        gw[col] = gw[col].clip(lower=0)

    print(f"  {len(gw):,} genes in genome-wide dataset")

    gw_csv = os.path.join(out_dir, "pbs_genomewide_all.csv")
    gw.to_csv(gw_csv, index=False)
    print(f"  Saved → {gw_csv}")

    # ── Genome-wide percentile ranks ──────────────────────────────────────────
    for col, _ in PBS_SCANS:
        gw[f"pct_{col}"] = gw[col].rank(pct=True, na_option="keep")

    # ── Merge with network gene list ──────────────────────────────────────────
    master_csv = os.path.join(data_dir, "network_constraint_gtex.csv")
    master = pd.read_csv(master_csv)
    master["gene_name"] = master["gene"].str.upper()

    # Also load targeted PBS for comparison
    targeted_csv = os.path.join(out_dir, "pbs_per_gene.csv")
    if os.path.exists(targeted_csv):
        targeted = pd.read_csv(targeted_csv)
        targeted["gene_name"] = targeted["gene"].str.upper()
        for col in ["pbs1_african", "pbs2_african", "pbs3_melanesian", "pbs4_melanesian"]:
            targeted[col] = targeted[col].clip(lower=0)
        targeted = targeted.rename(columns={c: f"targeted_{c}" for c, _ in PBS_SCANS})
    else:
        targeted = None

    network = master[["gene_name", "functional_category", "LOEUF"]].copy()
    network = network.merge(
        gw[["gene_name"] + [c for c, _ in PBS_SCANS] + [f"pct_{c}" for c, _ in PBS_SCANS]],
        on="gene_name", how="left"
    )
    if targeted is not None:
        network = network.merge(
            targeted[["gene_name"] + [f"targeted_{c}" for c, _ in PBS_SCANS]],
            on="gene_name", how="left"
        )

    pct_csv = os.path.join(out_dir, "pbs_genomewide_percentiles.csv")
    network.to_csv(pct_csv, index=False)
    print(f"  Saved → {pct_csv}")

    print("\n=== Melanogenesis network — genome-wide PBS percentiles ===")
    for col, label in PBS_SCANS:
        sub = network.dropna(subset=[f"pct_{col}"])
        med = sub[f"pct_{col}"].median()
        top10 = (sub[f"pct_{col}"] >= 0.90).sum()
        print(f"  {label.replace(chr(10),' ')}: median percentile = {med:.2f}, "
              f"genes ≥ 90th pct = {top10}/{len(sub)}")

    # ── Figure ────────────────────────────────────────────────────────────────
    # 4 panels (one per scan): genome-wide distribution (histogram/density) with
    # network genes overlaid as coloured rug ticks.
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    fig.subplots_adjust(hspace=0.40, wspace=0.30,
                        left=0.08, right=0.97, top=0.88, bottom=0.10)

    for ax, (col, label) in zip(axes, PBS_SCANS):
        gw_vals = gw[col].dropna()
        net_sub = network.dropna(subset=[col])

        # Genome-wide histogram
        ax.hist(gw_vals, bins=80, color="#CCCCCC", edgecolor="none",
                density=True, alpha=0.85, label=f"All genes (n={len(gw_vals):,})")

        # Vertical lines for 90th and 95th percentile thresholds
        p90 = gw_vals.quantile(0.90)
        p95 = gw_vals.quantile(0.95)
        ax.axvline(p90, color="#888888", lw=1.2, linestyle="--",
                   label=f"90th pct ({p90:.3f})")
        ax.axvline(p95, color="#555555", lw=1.2, linestyle=":",
                   label=f"95th pct ({p95:.3f})")

        # Rug ticks for network genes, coloured by functional category
        y_rug = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1.0
        for cat in CATEGORY_ORDER:
            sub = net_sub[net_sub["functional_category"] == cat]
            if len(sub) == 0:
                continue
            ax.scatter(sub[col],
                       np.full(len(sub), y_rug * 0.97),
                       color=CATEGORY_COLORS[cat],
                       s=40, marker="|", lw=1.5, zorder=5,
                       label=f"{cat} (n={len(sub)})")

        ax.set_xlabel("PBS value", fontsize=11)
        ax.set_ylabel("Density", fontsize=11)
        ax.set_title(label.replace("\n", " — "), fontsize=11,
                     fontweight="bold", loc="left")
        ax.tick_params(labelsize=9)

        # Inset: median percentile rank of network genes
        med_pct = net_sub[f"pct_{col}"].median() * 100
        n_top90 = (net_sub[f"pct_{col}"] >= 0.90).sum()
        ax.text(0.97, 0.97,
                f"Network median: {med_pct:.0f}th pct\n"
                f"Genes ≥ 90th pct: {n_top90}/{len(net_sub)}",
                transform=ax.transAxes, fontsize=8.5,
                va="top", ha="right",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="gray", alpha=0.9))

    # Shared category legend
    cat_patches = [
        mpatches.Patch(facecolor=CATEGORY_COLORS[c], label=c, alpha=0.85,
                       edgecolor="white")
        for c in CATEGORY_ORDER
    ]
    fig.legend(handles=cat_patches, fontsize=9, ncol=6,
               loc="lower center", bbox_to_anchor=(0.5, 0.002),
               framealpha=0.92, edgecolor="gray",
               title="Functional category (rug tick color)",
               title_fontsize=9)

    fig.suptitle(
        "Melanogenesis network PBS values in genome-wide context\n"
        "Grey histogram = all protein-coding genes  |  "
        "Coloured ticks = 129-gene network",
        fontsize=12, fontweight="bold", y=0.975)

    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, f"figure_pbs_genomewide_compare.{ext}")
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"\nSaved → {path}")
    plt.close(fig)

    print("\nDone.")


if __name__ == "__main__":
    main()
