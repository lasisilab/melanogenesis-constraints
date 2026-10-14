# Provenance note — within-network PBS τ results

**Status:** working note, 2026-10-14. Written to resolve the question of why the poster's
τ coefficients did not replicate on refitting. Every number below was regenerated from the
repo as committed; the scripts and output files are named inline so anyone can re-run them.

**Short answer:** the poster's numbers are correct and reproduce exactly from the committed
data. The refit that came back null was run against a *different dataset* that lives in this
repo under a nearly identical filename, and that dataset is already documented as null. There
is no data-integrity problem. There are, however, three separate statistical caveats that the
poster materials currently overstate or omit, and two of them are more serious than the
provenance question that started this.

---

## 1. The poster numbers reproduce exactly

Re-running the committed model — outcome = within-network PBS percentile, combined
HGDP+1KGP+SGDP callset, covariates `z(log1p n_sites) + z(log10 body length)` — gives:

| Branch | Predictor | adj β | parametric p | perm p |
|---|---|---|---|---|
| AFR | τ | −6.89 | 0.0088 | 0.011 |
| AFR | tissue breadth | +7.88 | 0.0027 | 0.003 |
| MEL | τ | +5.28 | 0.0473 | 0.048 |
| MEL | tissue breadth | −6.18 | 0.0203 | 0.021 |
| EUR | τ | −0.98 | 0.72 | ns |
| SAS | τ | −0.16 | 0.95 | ns |

These match [`output/phase2_pbs_within_network_combined.txt`](../output/phase2_pbs_within_network_combined.txt)
and the poster to the reported precision. Source:
[`analysis/phase2_pbs_within_network_combined.py`](../analysis/phase2_pbs_within_network_combined.py).
Nothing is missing, stale, or fabricated.

## 2. Why the refit came back null — two analyses, nearly identical names

The repo contains two parallel within-network PBS analyses:

| Script | Dataset | Outcome | AFR τ | MEL τ |
|---|---|---|---|---|
| `phase2_pbs_within_network_combined.py` | HGDP+1KGP+SGDP (AFR n=747, MEL n=47) | within-network percentile | **−6.89, p=0.011 \*** | **+5.28, p=0.048 \*** |
| `phase2_pbs_within_network.py` | SGDP only (AFR n≈18) | genome-wide percentile | +2.37, p=0.38 ns | −0.67, p=0.79 ns |

The SGDP-only version is null in all four branches with the τ sign flipped in both AFR and
MEL. That is the signature of the refit. The combined script's own docstring says this
outright:

> the SGDP-only genome-wide version of this analysis is NULL; this signal does not replicate
> at SGDP sample sizes. Report as dataset-dependent.

[`pages/popgen_selection.qmd`](../pages/popgen_selection.qmd) also carries the SGDP null as an
explicit robustness panel and attributes the flip to sample size rather than to the percentile
transform (correctly — Spearman ρ is rank-invariant, so the transform cannot produce it). The
disclosure exists in the Quarto site. It does **not** exist in
[`POSTER_METHODS_GUIDE.md`](../POSTER_METHODS_GUIDE.md), `PEQG_cheat_card.md`,
`PEQG_stats_deepdive.md`, or `CHANGELOG.md`, all of which state the combined numbers flat.

**`tau_pct` is not the cause.** That column only exists in
[`phase2_gene_profile_clustering.py`](../analysis/phase2_gene_profile_clustering.py) as a rank
score for the composite profile figure; it is never a regression predictor. Substituted into
the combined model it barely moves anything: AFR −6.09 (p=0.024), MEL +5.14 (p=0.053).
Rank-transforming τ did not kill the result.

**Recommended wording for the poster**, if the combined result stays as primary: *"Combined
callset (AFR n=747, MEL n=47). The SGDP-only subset (AFR n≈18) is null in all branches; we
report the better-powered dataset and treat this association as dataset-dependent."*

---

## 3. Caveat one — MEL τ p = 0.048 is seed-dependent (resolved: p = 0.0489)

The permutation p-value depends on the RNG seed, and the script shares one
`default_rng(0)` stream across all 20 branch×predictor loops, so each test's permutation draw
depends on loop order. Re-running each of the four headline tests across 40 *independent*
seeds × 10,000 permutations:

| Test | β | perm p min | max | mean | sd | fraction of seeds < 0.05 |
|---|---|---|---|---|---|---|
| **MEL τ** | +5.28 | 0.0460 | 0.0526 | **0.0490** | 0.0017 | **62%** |
| AFR τ | −6.89 | 0.0080 | 0.0117 | 0.0100 | 0.0009 | 100% |
| MEL breadth | −6.18 | 0.0188 | 0.0248 | 0.0213 | 0.0013 | 100% |
| AFR breadth | +7.88 | 0.0018 | 0.0047 | 0.0032 | 0.0006 | 100% |

MEL τ sits on 0.05 — at 10,000 permutations the asterisk on that panel is decided by the seed,
and 15 of 40 seeds take it away. The other three are stable. The stored value of 0.0475 is not
wrong; it is one draw from a distribution centred just under 0.05.

**Resolved at higher resolution.** Re-running MEL τ once at 200,000 permutations (seed
20261014) gives:

> **MEL τ, β = +5.28, permutation p = 0.04888.**

So it does clear 0.05 — but by 0.0011, with the Monte Carlo standard error on that estimate
itself around 0.0005. The asterisk survives and is defensible as printed. It is also the
thinnest result on the poster, and the margin is smaller than the width of the seed-to-seed
spread above. If anyone asks about that panel at PEQG, the honest line is that MEL τ is
significant at p = 0.049 on 200,000 permutations and should be read as suggestive: the MEL
branch rests on n = 47 samples (30 HGDP + 17 SGDP), and a p that lands one thousandth under
threshold is not a result to lean on. The AFR/MEL *contrast* in particular should not be
described as if both halves were equally supported.

(An earlier draft of this note reported the seed sweep as mean 0.0500 with 45% of seeds under
0.05, from a throwaway script. The committed script above gives mean 0.0490 and 62% — same
conclusion, and the 200,000-permutation value supersedes both.)

## 4. Caveat two — τ and tissue breadth are one finding, not two

τ and `tissue_breadth` correlate at **r = −0.890** (Pearson; Spearman −0.888, n=127). They
are two parameterisations of the same quantity. Fitting both in one model:

| Branch | Model | β(τ) | β(breadth) | R² |
|---|---|---|---|---|
| AFR | τ alone | −6.89 | — | 0.095 |
| AFR | breadth alone | — | +7.88 | 0.111 |
| AFR | both (VIF ≈ 5.1 each) | +0.59 (p=0.92) | +8.42 (p=0.14) | 0.112 |
| MEL | τ alone | +5.28 | — | 0.065 |
| MEL | breadth alone | — | −6.18 | 0.076 |
| MEL | both (VIF ≈ 5.1 each) | −1.04 (p=0.86) | −7.11 (p=0.22) | 0.076 |

Adding τ on top of breadth buys R² = 0.001 in AFR and 0.000 in MEL. Neither predictor
survives when both are present, which is the expected behaviour for two near-collinear
measures of the same axis — it does not mean both are false, it means they cannot be
separated at n=123.

The poster currently presents the τ and breadth coefficients as mutually corroborating
("β_τ = −6.89 ... β_breadth = +7.88"). They cannot corroborate each other. This matters for
the framing more than for the result: the finding is *one* association between the
tissue-specificity axis and branch-specific PBS, measured twice. It is defensible to report
both as a consistency check — but then say that is what it is.

`POSTER_METHODS_GUIDE.md` pre-empts the joint-model question by arguing the analysis is about
marginal association, which is a legitimate choice. The problem is pairing that choice with
language that treats the two coefficients as independent support.

## 5. Caveat three — the reported VIF reassures about the wrong thing

`POSTER_METHODS_GUIDE.md` states: *"all VIFs < 1.5, confirming the predictor and covariates
aren't badly collinear."* This is not what the script computes. It calls
`variance_inflation_factor(Xv, 1)` — position 1 only, i.e. the predictor. The covariates are
never checked.

On the actual analysis subset (n=123, after the `n_sites_shared > 0` filter),
`log1p(n_sites_shared)` and `log10(body length)` correlate at **r = 0.961**, giving
**VIF ≈ 13.6** for each. The signature shows in the fitted coefficients: in the AFR both-in
model the two covariates come out at −16.79 and +17.43, large and offsetting.

The predictor's own VIF really is ~1.04, so the τ estimate is not directly destabilised by
this. But the covariate adjustment itself is unstable, and the sentence in the methods guide
is not supportable as written. Either compute VIF across all positions and report what it
says, or drop one covariate — gene length and SNP count are largely redundant here.

---

## What I would change before the poster goes up

1. Add the dataset-dependence caveat to the poster materials; the Quarto site already has it.
2. ~~Re-run MEL τ with more permutations~~ — done (section 3): p = 0.04888 at 200,000
   permutations. The asterisk stays, but describe MEL τ as suggestive and stop presenting
   the AFR/MEL contrast as two equally supported halves.
3. Stop describing τ and breadth as two findings.
4. Correct or delete the VIF sentence in `POSTER_METHODS_GUIDE.md`.

None of this touches the AFR result, which is the strongest of the set: β = −6.89, perm
p = 0.011, stable across every seed tested, and backed by the breadth parameterisation at
p = 0.003. The AFR/MEL *contrast* is the part that now rests on a borderline MEL coefficient.

## Reproducing the numbers in this note

Sections 3–5 are generated by a committed script:

```bash
python analysis/phase2_pbs_within_network_robustness.py
```

It rebuilds the model from the same committed inputs as
[`phase2_pbs_within_network_combined.py`](../analysis/phase2_pbs_within_network_combined.py)
(same outcome, predictors, covariates and filter) and writes:

- [`output/phase2_pbs_within_network_robustness.txt`](../output/phase2_pbs_within_network_robustness.txt) — full report
- `output/phase2_pbs_within_network_robustness_seeds.csv` — the seed-sensitivity table
- `output/phase2_pbs_within_network_robustness_collinearity.csv` — the τ/breadth model comparison

It reads only; it does not modify the headline analysis or any poster material. Runtime is
about 10 seconds, including the 200,000-permutation MEL τ run — permutation coefficients are
computed in closed form (`pinv(X) @ y_perm`) rather than by refitting OLS, which is
algebraically identical (verified to 1e-15) and is what makes that resolution cheap.
