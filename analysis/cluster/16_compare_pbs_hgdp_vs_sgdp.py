"""
16_compare_pbs_hgdp_vs_sgdp.py

Compares the HGDP-only and SGDP-only genome-wide PBS results gene-by-gene.

Inputs (in BASE/output/):
  pbs_genomewide_chr{1..22}.csv         — HGDP-only (from 11_compute_pbs_genomewide.py)
  pbs_sgdp_genomewide_chr{1..22}.csv    — SGDP-only (from 15_compute_pbs_sgdp_genomewide.py)
  network_constraint_categorized.csv    — 129 network genes for highlighting

Outputs (in BASE/output/):
  pbs_hgdp_vs_sgdp_per_gene.csv         — merged gene-level PBS table
  figure_pbs_hgdp_vs_sgdp_scatter.png   — 2×2 scatter: pbs1/pbs3 each, network
                                          genes highlighted, Spearman ρ + p
                                          per panel
  pbs_hgdp_vs_sgdp_summary.txt          — text summary (correlations,
                                          top-percentile concordance, top-10
                                          network genes in each dataset)

Headline questions answered:
  1. Genome-wide Spearman ρ between HGDP and SGDP PBS at the gene level
     → "Do the two datasets generally agree about which genes are selected?"
  2. Network-only Spearman ρ
     → "Do the two datasets agree on the melanogenesis network?"
  3. Top-1% concordance (Jaccard overlap of top-1% gene sets)
     → "Do the strongest signals replicate across datasets?"

Run locally after both pipelines have completed and the per-chrom CSVs have
been pulled back from the cluster.

Usage:
    python analysis/cluster/16_compare_pbs_hgdp_vs_sgdp.py \\
        --base /path/to/repo
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=os.getcwd(),
                   help="Repo root (contains output/ and data/).")
    return p.parse_args()


def load_concat(prefix, out_dir):
    """Concatenate pbs_*_chr{1..22}.csv into one frame; missing chrs are
    skipped with a warning so partial pipelines still work."""
    pieces = []
    for chrom in range(1, 23):
        path = os.path.join(out_dir, f"{prefix}{chrom}.csv")
        if not os.path.exists(path):
            print(f"  WARNING: missing {path}", file=sys.stderr)
            continue
        pieces.append(pd.read_csv(path))
    if not pieces:
        raise FileNotFoundError(f"No per-chrom CSVs found for prefix '{prefix}' "
                                f"in {out_dir}.")
    return pd.concat(pieces, ignore_index=True)


def top_pct_jaccard(a_rank, b_rank, pct):
    """Jaccard overlap of the top-`pct` of two ranked sets (indexed identically)."""
    n = len(a_rank)
    k = max(1, int(round(n * pct / 100)))
    top_a = set(np.argsort(-a_rank)[:k])
    top_b = set(np.argsort(-b_rank)[:k])
    j = len(top_a & top_b) / len(top_a | top_b)
    return k, j


def main():
    args = parse_args()
    out_dir  = os.path.join(args.base, "output")
    data_dir = os.path.join(args.base, "data")

    print("Loading HGDP-only PBS results...")
    hgdp = load_concat("pbs_genomewide_chr",       out_dir)
    print(f"  {len(hgdp):,} HGDP gene rows")

    print("Loading SGDP-only PBS results...")
    sgdp = load_concat("pbs_sgdp_genomewide_chr",  out_dir)
    print(f"  {len(sgdp):,} SGDP gene rows")

    pbs_cols = ["pbs1_african", "pbs2_african",
                "pbs3_melanesian", "pbs4_melanesian"]
    merged = hgdp[["gene_name", "gene_id", "chrom", "start", "end",
                   "n_snps"] + pbs_cols].rename(
        columns={c: c + "_hgdp" for c in pbs_cols + ["n_snps"]}
    ).merge(
        sgdp[["gene_name", "gene_id"] + pbs_cols + ["n_snps"]].rename(
            columns={c: c + "_sgdp" for c in pbs_cols + ["n_snps"]}
        ),
        on=["gene_name", "gene_id"], how="inner",
    )
    print(f"  {len(merged):,} genes present in both pipelines")

    # Network gene set (for highlight + per-network correlation)
    network_csv = os.path.join(data_dir, "network_constraint_categorized.csv")
    if os.path.exists(network_csv):
        network = set(pd.read_csv(network_csv)["gene"].str.upper().tolist())
        merged["is_network"] = merged["gene_name"].str.upper().isin(network)
    else:
        merged["is_network"] = False
        print(f"  WARNING: {network_csv} not found; network highlight disabled.")

    out_csv = os.path.join(out_dir, "pbs_hgdp_vs_sgdp_per_gene.csv")
    merged.to_csv(out_csv, index=False)
    print(f"Saved → {out_csv}")

    # ── Stats ────────────────────────────────────────────────────────────
    summary = ["HGDP vs. SGDP genome-wide PBS — comparison summary",
               "=" * 60,
               f"Genes in both: {len(merged):,}",
               f"Network genes in both: {int(merged['is_network'].sum())}",
               ""]

    panels = [
        ("pbs1_african",    "PBS-1  African (target) / S.Asian (out) / Melanesian (dist)"),
        ("pbs2_african",    "PBS-2  African (target) / European (out) / Melanesian (dist)"),
        ("pbs3_melanesian", "PBS-3  Melanesian (target) / S.Asian (out) / African (dist)"),
        ("pbs4_melanesian", "PBS-4  Melanesian (target) / European (out) / African (dist)"),
    ]

    summary.append("Spearman ρ (HGDP vs SGDP) per PBS scan:")
    summary.append(f"  {'scan':<25} {'all genes':>12} {'network':>12}")
    panel_stats = {}
    for col, label in panels:
        a = merged[col + "_hgdp"].values
        b = merged[col + "_sgdp"].values
        mask = np.isfinite(a) & np.isfinite(b)
        rho_all, p_all = stats.spearmanr(a[mask], b[mask])
        net = merged["is_network"] & mask
        if net.sum() >= 5:
            rho_net, p_net = stats.spearmanr(a[net], b[net])
        else:
            rho_net, p_net = (np.nan, np.nan)
        panel_stats[col] = dict(rho_all=rho_all, p_all=p_all,
                                rho_net=rho_net, p_net=p_net,
                                n_all=int(mask.sum()),
                                n_net=int(net.sum()))
        summary.append(f"  {label[:25]:<25}  ρ = {rho_all:+.3f} (p={p_all:.1e})"
                       f"   ρ_net = {rho_net:+.3f} (p={p_net:.1e})")
    summary.append("")

    # Top-percentile concordance for pbs1 + pbs3 (the two used in poster)
    summary.append("Top-percentile Jaccard concordance (all genes):")
    for col, label in [("pbs1_african",    "PBS-1 African"),
                       ("pbs3_melanesian", "PBS-3 Melanesian")]:
        a = merged[col + "_hgdp"].fillna(-np.inf).values
        b = merged[col + "_sgdp"].fillna(-np.inf).values
        for pct in (0.5, 1.0, 5.0, 10.0):
            k, j = top_pct_jaccard(a, b, pct)
            summary.append(f"  {label}  top {pct:>4.1f}%  (n={k})   "
                           f"Jaccard = {j:.3f}")
    summary.append("")

    # Top-10 network genes per dataset (per axis)
    if merged["is_network"].any():
        net_df = merged.loc[merged["is_network"]].copy()
        for col, label in [("pbs1_african",    "PBS-1 African"),
                           ("pbs3_melanesian", "PBS-3 Melanesian")]:
            top_hgdp = net_df.nlargest(10, col + "_hgdp")[
                ["gene_name", col + "_hgdp", col + "_sgdp"]]
            top_sgdp = net_df.nlargest(10, col + "_sgdp")[
                ["gene_name", col + "_hgdp", col + "_sgdp"]]
            summary.append(f"Top 10 network genes by {label} (HGDP):")
            summary.append(top_hgdp.to_string(index=False))
            summary.append("")
            summary.append(f"Top 10 network genes by {label} (SGDP):")
            summary.append(top_sgdp.to_string(index=False))
            summary.append("")

    summary_path = os.path.join(out_dir, "pbs_hgdp_vs_sgdp_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary) + "\n")
    print(f"Saved → {summary_path}")

    # ── 2×2 scatter figure ──────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    plot_panels = [
        ("pbs1_african",    "PBS-1  African (SAS out, MEL dist)", axes[0, 0]),
        ("pbs2_african",    "PBS-2  African (EUR out, MEL dist)", axes[0, 1]),
        ("pbs3_melanesian", "PBS-3  Melanesian (SAS out, AFR dist)", axes[1, 0]),
        ("pbs4_melanesian", "PBS-4  Melanesian (EUR out, AFR dist)", axes[1, 1]),
    ]
    for col, label, ax in plot_panels:
        a = merged[col + "_hgdp"].values
        b = merged[col + "_sgdp"].values
        ax.scatter(a, b, s=5, alpha=0.15, color="#888",
                   edgecolors="none", zorder=2)
        if merged["is_network"].any():
            net = merged["is_network"].values
            ax.scatter(a[net], b[net], s=42, color="#D94040",
                       edgecolors="black", linewidths=0.4, alpha=0.95,
                       zorder=3, label=f"Network (n={int(net.sum())})")
        lim_max = float(np.nanmax([a, b]))
        ax.plot([0, lim_max], [0, lim_max], "--", color="#555", lw=0.8,
                zorder=1)
        s = panel_stats[col]
        ax.text(0.04, 0.96,
                f"ρ_all = {s['rho_all']:+.3f}  (n = {s['n_all']:,})\n"
                f"ρ_net = {s['rho_net']:+.3f}  (n = {s['n_net']})",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3",
                          facecolor="white", edgecolor="#aaa", lw=0.6))
        ax.set_xlabel("HGDP PBS", fontsize=11)
        ax.set_ylabel("SGDP PBS", fontsize=11)
        ax.set_title(label, fontsize=11, fontweight="bold", loc="left")
        if merged["is_network"].any():
            ax.legend(loc="lower right", fontsize=9, frameon=True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.suptitle("Genome-wide PBS — HGDP-only vs. SGDP-only "
                 "(per-gene concordance)",
                 fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout()
    fig_path = os.path.join(out_dir, "figure_pbs_hgdp_vs_sgdp_scatter.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(fig_path.replace(".png", ".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"Saved → {fig_path}  (+ .pdf)")


if __name__ == "__main__":
    main()
