"""
17_network_enrichment_genomewide.py

Tests whether the melanogenesis network genes carry unusually strong PBS
signal relative to the rest of the genome-wide scan — using ONLY a single
genome-wide scan (no targeted/HGDP-restricted scan needed).

For each PBS scan it computes:
  1. Per-network-gene genome-wide percentile
       "MC1R is at the 99.4th percentile of genome-wide PBS-3."
       → replaces the within-network 75th-percentile threshold with a
         genome-wide-calibrated one.
  2. Mann-Whitney U  (network vs. all non-network genes, one-sided 'greater')
       Fast, standard, but does NOT control for gene length / SNP count.
  3. SNP-count-matched permutation test
       Draws n_perm random gene sets from the non-network background, matched
       to the network's n_snps distribution (decile bins), builds a null of
       the median PBS, and reports an empirical p-value. Controls for the
       "longer genes have higher/steadier PBS" confound.

Inputs (in BASE/output/):
  {prefix}{1..22}.csv                    — per-chrom genome-wide PBS
                                           (default prefix: pbs_sgdp_genomewide_chr)
  data/network_constraint_categorized.csv — melanogenesis network gene list

Outputs (in BASE/output/), tagged with --label:
  network_enrichment_{label}_summary.txt
  network_gene_genomewide_percentiles_{label}.csv
  figure_network_enrichment_{label}.png / .pdf

Usage:
  # SGDP (data already present locally)
  python analysis/cluster/17_network_enrichment_genomewide.py \\
      --base . --prefix pbs_sgdp_genomewide_chr --label SGDP

  # HGDP (once pbs_genomewide_chr*.csv exist)
  python analysis/cluster/17_network_enrichment_genomewide.py \\
      --base . --prefix pbs_genomewide_chr --label HGDP
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


SCANS = [
    ("pbs1_african",    "PBS-1  African (S.Asian out, Melanesian dist)"),
    ("pbs2_african",    "PBS-2  African (European out, Melanesian dist)"),
    ("pbs3_melanesian", "PBS-3  Melanesian (S.Asian out, African dist)"),
    ("pbs4_melanesian", "PBS-4  Melanesian (European out, African dist)"),
]
PRIMARY_SCANS = ["pbs1_african", "pbs3_melanesian"]  # plotted; all 4 in the text summary


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=os.getcwd())
    p.add_argument("--prefix", default="pbs_sgdp_genomewide_chr",
                   help="Per-chrom CSV prefix (without the chromosome number).")
    p.add_argument("--label", default="SGDP",
                   help="Short tag for output filenames / titles.")
    p.add_argument("--n-perm", type=int, default=10000,
                   help="Permutations for the SNP-matched null.")
    p.add_argument("--n-bins", type=int, default=10,
                   help="n_snps quantile bins for matching (deciles by default).")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def load_genomewide(prefix, out_dir):
    pieces = []
    for chrom in range(1, 23):
        path = os.path.join(out_dir, f"{prefix}{chrom}.csv")
        if not os.path.exists(path):
            print(f"  WARNING: missing {path}", file=sys.stderr)
            continue
        pieces.append(pd.read_csv(path))
    if not pieces:
        raise FileNotFoundError(
            f"No per-chrom CSVs found for prefix '{prefix}' in {out_dir}."
        )
    return pd.concat(pieces, ignore_index=True)


def snp_matched_permutation(df, scan_col, is_net, n_perm, n_bins, rng):
    """Empirical p for 'network median PBS > matched-background median PBS'.

    Background gene sets are drawn from non-network genes, matched to the
    network's n_snps distribution via quantile bins. Returns
    (observed_median, p_value, null_medians)."""
    valid = df[scan_col].notna() & (df["n_snps"] > 0)
    sub = df[valid].copy()
    net_mask = is_net[valid].values

    net_pbs = sub.loc[net_mask, scan_col].values
    observed = np.median(net_pbs)
    n_net = len(net_pbs)

    # Quantile bins on n_snps over ALL valid genes
    # (so network and background share the same bin edges).
    snps = sub["n_snps"].values
    edges = np.quantile(snps, np.linspace(0, 1, n_bins + 1))
    edges[-1] += 1  # make the top edge inclusive
    bin_idx = np.clip(np.digitize(snps, edges[1:-1]), 0, n_bins - 1)

    bg_mask = ~net_mask
    # How many network genes fall in each bin
    net_bin_counts = np.bincount(bin_idx[net_mask], minlength=n_bins)
    # Background gene row-indices per bin
    bg_by_bin = {b: np.where(bg_mask & (bin_idx == b))[0] for b in range(n_bins)}
    pbs_all = sub[scan_col].values

    null = np.empty(n_perm)
    for i in range(n_perm):
        drawn = []
        for b in range(n_bins):
            k = net_bin_counts[b]
            if k == 0:
                continue
            pool = bg_by_bin[b]
            if len(pool) == 0:
                continue
            replace = len(pool) < k
            drawn.append(rng.choice(pool, size=k, replace=replace))
        idx = np.concatenate(drawn) if drawn else np.array([], dtype=int)
        null[i] = np.median(pbs_all[idx]) if len(idx) else np.nan

    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return observed, p, null, n_net


def main():
    args = parse_args()
    out_dir  = os.path.join(args.base, "output")
    data_dir = os.path.join(args.base, "data")
    rng = np.random.default_rng(args.seed)

    print(f"Loading genome-wide scan ({args.label}, prefix '{args.prefix}')...")
    gw = load_genomewide(args.prefix, out_dir)
    print(f"  {len(gw):,} gene rows")

    # Flag network genes
    network_csv = os.path.join(data_dir, "network_constraint_categorized.csv")
    if not os.path.exists(network_csv):
        sys.exit(f"ERROR: {network_csv} not found — need the network gene list.")
    network = set(pd.read_csv(network_csv)["gene"].str.upper().tolist())
    gw["is_network"] = gw["gene_name"].str.upper().isin(network)
    n_net_total = int(gw["is_network"].sum())
    print(f"  {n_net_total} network genes matched in the genome-wide scan "
          f"(of {len(network)} in the network list)")

    is_net = gw["is_network"]

    summary = [
        f"Melanogenesis network PBS enrichment vs. genome-wide background "
        f"[{args.label}]",
        "=" * 70,
        f"Genome-wide genes: {len(gw):,}",
        f"Network genes matched: {n_net_total}",
        f"Permutations: {args.n_perm:,}   (SNP-matched, {args.n_bins} bins)",
        "",
    ]

    percentile_rows = {}  # gene_name → {scan: percentile}

    for scan_col, scan_label in SCANS:
        valid = gw[scan_col].notna() & (gw["n_snps"] > 0)
        sub = gw[valid]
        net_vals = sub.loc[sub["is_network"], scan_col].values
        bg_vals  = sub.loc[~sub["is_network"], scan_col].values
        if len(net_vals) < 5 or len(bg_vals) < 50:
            summary.append(f"{scan_label}: insufficient data, skipped.")
            continue

        # 1. Per-network-gene genome-wide percentile (vs ALL valid genes)
        all_vals = sub[scan_col].values
        order = np.argsort(all_vals)
        ranks = np.empty(len(all_vals))
        ranks[order] = np.arange(len(all_vals))
        pct = 100.0 * ranks / (len(all_vals) - 1)
        net_names = sub.loc[sub["is_network"], "gene_name"].values
        net_pct = pct[sub["is_network"].values]
        for g, pc, v in zip(net_names, net_pct, net_vals):
            percentile_rows.setdefault(g, {})[scan_col] = pc

        # 2. Mann-Whitney U (one-sided: network > background)
        u, p_mwu = stats.mannwhitneyu(net_vals, bg_vals, alternative="greater")

        # 3. SNP-matched permutation
        observed, p_perm, null, n_net = snp_matched_permutation(
            gw, scan_col, is_net, args.n_perm, args.n_bins, rng)

        net_med, bg_med = np.median(net_vals), np.median(bg_vals)
        null_med = np.nanmedian(null)
        n_top1 = int((net_pct >= 99).sum())
        n_top5 = int((net_pct >= 95).sum())

        summary += [
            scan_label,
            f"  network median PBS   = {net_med:.4f}",
            f"  background median    = {bg_med:.4f}",
            f"  matched-null median  = {null_med:.4f}",
            f"  Mann-Whitney U       p = {p_mwu:.3e}  (one-sided, network>bg)",
            f"  SNP-matched perm     p = {p_perm:.4f}  "
            f"(network median vs {args.n_perm:,} matched draws)",
            f"  network genes in genome-wide top 1%: {n_top1} / {n_net} "
            f"(expected ~{n_net*0.01:.1f})",
            f"  network genes in genome-wide top 5%: {n_top5} / {n_net} "
            f"(expected ~{n_net*0.05:.1f})",
            "",
        ]

    # Per-gene percentile table
    pct_df = pd.DataFrame.from_dict(percentile_rows, orient="index")
    pct_df = pct_df.rename(columns={c: f"{c}_pctile" for c in pct_df.columns})
    pct_df.index.name = "gene_name"
    pct_df = pct_df.sort_values("pbs3_melanesian_pctile", ascending=False) \
        if "pbs3_melanesian_pctile" in pct_df.columns else pct_df
    pct_csv = os.path.join(out_dir,
                           f"network_gene_genomewide_percentiles_{args.label}.csv")
    pct_df.to_csv(pct_csv)
    print(f"Saved → {pct_csv}")

    # Top network genes by genome-wide percentile (PBS-1 and PBS-3)
    for scan_col in PRIMARY_SCANS:
        col = f"{scan_col}_pctile"
        if col in pct_df.columns:
            top = pct_df.sort_values(col, ascending=False).head(15)
            summary.append(f"Top 15 network genes by genome-wide {scan_col} percentile:")
            for g, r in top.iterrows():
                summary.append(f"  {g:<12} {r[col]:6.2f} pctile")
            summary.append("")

    summary_path = os.path.join(out_dir,
                                f"network_enrichment_{args.label}_summary.txt")
    with open(summary_path, "w") as f:
        f.write("\n".join(summary) + "\n")
    print(f"Saved → {summary_path}")
    print("\n" + "\n".join(summary))

    # ── Figure: ECDF network vs background, PBS-1 and PBS-3 ─────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, scan_col in zip(axes, PRIMARY_SCANS):
        scan_label = dict(SCANS)[scan_col]
        valid = gw[scan_col].notna() & (gw["n_snps"] > 0)
        sub = gw[valid]
        net_vals = np.sort(sub.loc[sub["is_network"], scan_col].values)
        bg_vals  = np.sort(sub.loc[~sub["is_network"], scan_col].values)

        def ecdf(v):
            return v, np.arange(1, len(v) + 1) / len(v)

        bx, by = ecdf(bg_vals)
        nx, ny = ecdf(net_vals)
        ax.plot(bx, by, color="#888", lw=1.6, label=f"Genome-wide (n={len(bg_vals):,})")
        ax.plot(nx, ny, color="#D94040", lw=2.2,
                label=f"Network (n={len(net_vals)})")
        ax.axvline(np.median(bg_vals),  color="#888", ls=":", lw=1)
        ax.axvline(np.median(net_vals), color="#D94040", ls=":", lw=1)

        u, p_mwu = stats.mannwhitneyu(net_vals, bg_vals, alternative="greater")
        # x-limit: clip to 99.5th pct of background so the bulk is visible
        ax.set_xlim(0, np.quantile(bg_vals, 0.995) * 1.5 + 1e-6)
        ax.set_xlabel(f"{scan_col}", fontsize=11)
        ax.set_ylabel("Cumulative fraction of genes", fontsize=11)
        ax.set_title(scan_label, fontsize=11, fontweight="bold", loc="left")
        ax.text(0.97, 0.05,
                f"MWU p = {p_mwu:.1e}\n"
                f"network median = {np.median(net_vals):.3f}\n"
                f"background median = {np.median(bg_vals):.3f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#aaa", lw=0.6))
        ax.legend(loc="lower right", fontsize=9, frameon=True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.suptitle(f"Melanogenesis network vs. genome-wide PBS  [{args.label}]\n"
                 "network ECDF shifted right of background = enrichment for "
                 "high PBS",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig_path = os.path.join(out_dir, f"figure_network_enrichment_{args.label}.png")
    fig.savefig(fig_path, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(fig_path.replace(".png", ".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"Saved → {fig_path}  (+ .pdf)")


if __name__ == "__main__":
    main()
