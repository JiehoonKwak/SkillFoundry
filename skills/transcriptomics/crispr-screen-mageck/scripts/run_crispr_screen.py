#!/usr/bin/env python3
"""Analyze a CRISPR screen count table using MAGeCK-style statistics.

Takes guide-level counts across treatment and control samples, performs
median-ratio normalization, computes guide-level fold changes, and
aggregates to gene-level scores using robust rank aggregation (RRA).

Outputs a ranked gene list with p-values, FDR, and summary statistics.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_count_table(path: Path) -> tuple[list[str], list[str], list[str], list[list[float]]]:
    """Load a MAGeCK-format count table.

    Expected format: sgRNA \\t gene \\t sample1 \\t sample2 ...
    Returns (samples, guides, genes, counts_matrix).
    counts_matrix[guide_idx][sample_idx].
    """
    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)

    if len(header) < 3:
        raise ValueError(
            f"Count table must have >= 3 columns (sgRNA, gene, samples), got {len(header)}"
        )

    samples = header[2:]
    guides: list[str] = []
    genes: list[str] = []
    counts: list[list[float]] = []

    with open(path, encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Row for guide {row[0]!r} has {len(row)} fields, expected {len(header)}"
                )
            guides.append(row[0])
            genes.append(row[1])
            counts.append([float(x) for x in row[2:]])

    return samples, guides, genes, counts


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------


def median_ratio_normalize(
    counts: list[list[float]],
    pseudocount: float = 0.5,
) -> tuple[list[list[float]], list[float]]:
    """Median-ratio normalization (DESeq2 / MAGeCK style).

    1. Compute geometric mean per guide across samples.
    2. Ratio = (count + pseudo) / geometric mean.
    3. Size factor per sample = median of ratios.
    4. Normalized count = (count + pseudo) / size factor.
    """
    n_guides = len(counts)
    n_samples = len(counts[0]) if counts else 0

    # Geometric mean per guide
    geo_means: list[float] = []
    for g in range(n_guides):
        log_sum = 0.0
        for s in range(n_samples):
            log_sum += math.log(counts[g][s] + pseudocount)
        geo_means.append(math.exp(log_sum / n_samples))

    # Size factors
    size_factors: list[float] = []
    for s in range(n_samples):
        ratios = []
        for g in range(n_guides):
            if geo_means[g] > 0:
                ratios.append((counts[g][s] + pseudocount) / geo_means[g])
        ratios.sort()
        median = ratios[len(ratios) // 2] if ratios else 1.0
        size_factors.append(median)

    # Normalize
    normalized: list[list[float]] = []
    for g in range(n_guides):
        row = []
        for s in range(n_samples):
            sf = size_factors[s] if size_factors[s] > 0 else 1.0
            row.append((counts[g][s] + pseudocount) / sf)
        normalized.append(row)

    return normalized, size_factors


# ---------------------------------------------------------------------------
# Guide-level statistics
# ---------------------------------------------------------------------------


@dataclass
class GuideResult:
    guide: str
    gene: str
    control_mean: float
    treatment_mean: float
    log2fc: float
    p_value: float  # one-sided, testing for depletion


def compute_guide_stats(
    normalized: list[list[float]],
    guides: list[str],
    genes: list[str],
    control_idx: list[int],
    treatment_idx: list[int],
) -> list[GuideResult]:
    """Compute per-guide log2 fold change and p-value.

    Uses a modified Z-test based on the mean-variance model:
    variance is estimated from the control samples with a shrinkage prior.
    """
    results: list[GuideResult] = []
    n_ctrl = len(control_idx)
    n_treat = len(treatment_idx)

    # Estimate global mean-variance relationship from controls
    # (simplified: use pooled variance estimate with empirical Bayes shrinkage)
    all_ctrl_means: list[float] = []
    all_ctrl_vars: list[float] = []
    for g in range(len(normalized)):
        vals = [normalized[g][s] for s in control_idx]
        mu = sum(vals) / n_ctrl
        all_ctrl_means.append(mu)
        if n_ctrl > 1:
            var = sum((v - mu) ** 2 for v in vals) / (n_ctrl - 1)
        else:
            var = mu  # Poisson assumption when only 1 replicate
        all_ctrl_vars.append(var)

    # Prior variance: median of all variances (empirical Bayes shrinkage)
    sorted_vars = sorted(all_ctrl_vars)
    prior_var = sorted_vars[len(sorted_vars) // 2] if sorted_vars else 1.0
    # Shrinkage weight
    shrink_weight = 0.5

    for g in range(len(normalized)):
        ctrl_vals = [normalized[g][s] for s in control_idx]
        treat_vals = [normalized[g][s] for s in treatment_idx]

        ctrl_mean = sum(ctrl_vals) / n_ctrl
        treat_mean = sum(treat_vals) / n_treat

        # Shrunken variance
        raw_var = all_ctrl_vars[g]
        shrunk_var = shrink_weight * prior_var + (1 - shrink_weight) * raw_var
        shrunk_var = max(shrunk_var, 1.0)  # floor

        # Standard error of the difference in means
        se = math.sqrt(shrunk_var / n_ctrl + shrunk_var / n_treat)

        # Log2 fold change
        if ctrl_mean > 0 and treat_mean > 0:
            lfc = math.log2(treat_mean / ctrl_mean)
        elif treat_mean <= 0:
            lfc = -10.0  # cap
        else:
            lfc = 10.0

        # Z-score (negative = depleted in treatment)
        z = (treat_mean - ctrl_mean) / se if se > 0 else 0.0

        # One-sided p-value for depletion (treatment < control)
        p_depletion = _normal_cdf(z)

        results.append(
            GuideResult(
                guide=guides[g],
                gene=genes[g],
                control_mean=round(ctrl_mean, 4),
                treatment_mean=round(treat_mean, 4),
                log2fc=round(lfc, 4),
                p_value=p_depletion,
            )
        )

    return results


def _normal_cdf(x: float) -> float:
    """Standard normal CDF approximation (Abramowitz & Stegun 26.2.17)."""
    if x > 8:
        return 1.0
    if x < -8:
        return 0.0
    sign = 1.0 if x >= 0 else -1.0
    x = abs(x)
    t = 1.0 / (1.0 + 0.2316419 * x)
    poly = t * (
        0.319381530
        + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
    pdf = math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)
    cdf = 1.0 - pdf * poly
    return 0.5 + sign * (cdf - 0.5)


# ---------------------------------------------------------------------------
# Gene-level aggregation: Robust Rank Aggregation (RRA)
# ---------------------------------------------------------------------------


@dataclass
class GeneResult:
    gene: str
    n_guides: int
    neg_score: float  # RRA score for negative selection (depletion)
    neg_p_value: float
    neg_fdr: float = 0.0
    neg_rank: int = 0
    neg_lfc: float = 0.0  # median log2FC across guides
    pos_score: float = 0.0  # RRA score for positive selection (enrichment)
    pos_p_value: float = 0.0
    pos_fdr: float = 0.0
    pos_rank: int = 0
    pos_lfc: float = 0.0
    guide_details: list[dict[str, Any]] = field(default_factory=list)


def robust_rank_aggregation(
    guide_results: list[GuideResult],
) -> list[GeneResult]:
    """Aggregate guide-level p-values to gene-level using RRA.

    For each gene:
    1. Rank the guide p-values.
    2. Compute the minimum of beta(k, n) p-values for the k-th ranked guide.
    3. Use the minimum as the gene-level score.
    """
    # Group guides by gene
    gene_guides: dict[str, list[GuideResult]] = {}
    for gr in guide_results:
        gene_guides.setdefault(gr.gene, []).append(gr)

    gene_results: list[GeneResult] = []
    n_total_guides = len(guide_results)

    for gene, guides in gene_guides.items():
        n = len(guides)

        # --- Negative selection (depletion) ---
        sorted_neg = sorted(guides, key=lambda g: g.p_value)
        neg_scores: list[float] = []
        for k, gr in enumerate(sorted_neg, 1):
            # Beta(k, n-k+1) CDF at the guide p-value = regularized incomplete beta
            beta_p = _beta_cdf(gr.p_value, k, n - k + 1)
            neg_scores.append(beta_p)
        neg_rra = min(neg_scores) * n if neg_scores else 1.0
        neg_rra = min(neg_rra, 1.0)

        # Approximate p-value: Bonferroni-corrected minimum
        neg_p = min(neg_rra, 1.0)

        # Median log2FC
        lfcs = sorted(g.log2fc for g in guides)
        neg_lfc = lfcs[len(lfcs) // 2]

        # --- Positive selection (enrichment) ---
        pos_pvals = [1.0 - g.p_value for g in guides]
        sorted_pos_idx = sorted(range(n), key=lambda i: pos_pvals[i])
        pos_scores: list[float] = []
        for k_idx, i in enumerate(sorted_pos_idx, 1):
            beta_p = _beta_cdf(pos_pvals[i], k_idx, n - k_idx + 1)
            pos_scores.append(beta_p)
        pos_rra = min(pos_scores) * n if pos_scores else 1.0
        pos_rra = min(pos_rra, 1.0)
        pos_p = min(pos_rra, 1.0)
        pos_lfc = lfcs[len(lfcs) // 2]

        guide_details = [
            {
                "guide": g.guide,
                "control_mean": g.control_mean,
                "treatment_mean": g.treatment_mean,
                "log2fc": g.log2fc,
                "p_value": round(g.p_value, 6),
            }
            for g in sorted_neg
        ]

        gene_results.append(
            GeneResult(
                gene=gene,
                n_guides=n,
                neg_score=neg_rra,
                neg_p_value=neg_p,
                neg_lfc=neg_lfc,
                pos_score=pos_rra,
                pos_p_value=pos_p,
                pos_lfc=pos_lfc,
                guide_details=guide_details,
            )
        )

    # FDR correction (Benjamini-Hochberg) for both directions
    _bh_correction(gene_results, direction="neg")
    _bh_correction(gene_results, direction="pos")

    # Rank by negative selection score
    gene_results.sort(key=lambda g: g.neg_score)
    for i, gr in enumerate(gene_results, 1):
        gr.neg_rank = i

    # Also assign positive ranks
    by_pos = sorted(gene_results, key=lambda g: g.pos_score)
    for i, gr in enumerate(by_pos, 1):
        gr.pos_rank = i

    # Re-sort by negative rank for output
    gene_results.sort(key=lambda g: g.neg_rank)

    return gene_results


def _bh_correction(results: list[GeneResult], direction: str = "neg") -> None:
    """Benjamini-Hochberg FDR correction in place."""
    if direction == "neg":
        pvals = [(r.neg_p_value, i) for i, r in enumerate(results)]
    else:
        pvals = [(r.pos_p_value, i) for i, r in enumerate(results)]

    pvals.sort()
    n = len(pvals)
    fdr_values = [0.0] * n

    prev_fdr = 1.0
    for rank_idx in range(n - 1, -1, -1):
        p, orig_idx = pvals[rank_idx]
        fdr = p * n / (rank_idx + 1)
        fdr = min(fdr, prev_fdr, 1.0)
        fdr_values[orig_idx] = fdr
        prev_fdr = fdr

    for i, r in enumerate(results):
        if direction == "neg":
            r.neg_fdr = round(fdr_values[i], 6)
        else:
            r.pos_fdr = round(fdr_values[i], 6)


def _beta_cdf(x: float, a: float, b: float) -> float:
    """Regularized incomplete beta function I_x(a, b).

    Uses a continued fraction expansion (Lentz's method).
    """
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0

    # Use symmetry relation when appropriate
    if x > (a + 1) / (a + b + 2):
        return 1.0 - _beta_cdf(1.0 - x, b, a)

    # Log of the prefactor: x^a * (1-x)^b / (a * B(a,b))
    ln_prefactor = (
        a * math.log(x) + b * math.log(1.0 - x) - math.log(a) - _log_beta(a, b)
    )

    # Continued fraction (Lentz's method)
    cf = _beta_cf(x, a, b)
    return math.exp(ln_prefactor) * cf


def _log_beta(a: float, b: float) -> float:
    """Log of beta function B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)."""
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _beta_cf(x: float, a: float, b: float, max_iter: int = 200, eps: float = 1e-14) -> float:
    """Continued fraction for incomplete beta (modified Lentz)."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d

    for m in range(1, max_iter + 1):
        m2 = 2 * m
        # Even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c

        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta

        if abs(delta - 1.0) < eps:
            break

    return h


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def format_output(
    gene_results: list[GeneResult],
    size_factors: list[float],
    samples: list[str],
    control_samples: list[str],
    treatment_samples: list[str],
) -> dict[str, Any]:
    """Format results as JSON payload."""
    return {
        "analysis": "crispr-screen-mageck",
        "normalization": "median-ratio",
        "control_samples": control_samples,
        "treatment_samples": treatment_samples,
        "size_factors": {s: round(sf, 4) for s, sf in zip(samples, size_factors)},
        "n_genes": len(gene_results),
        "n_guides": sum(g.n_guides for g in gene_results),
        "gene_summary": [
            {
                "gene": g.gene,
                "n_guides": g.n_guides,
                "neg": {
                    "rank": g.neg_rank,
                    "score": round(g.neg_score, 6),
                    "p_value": round(g.neg_p_value, 6),
                    "fdr": g.neg_fdr,
                    "lfc": g.neg_lfc,
                },
                "pos": {
                    "rank": g.pos_rank,
                    "score": round(g.pos_score, 6),
                    "p_value": round(g.pos_p_value, 6),
                    "fdr": g.pos_fdr,
                    "lfc": g.pos_lfc,
                },
                "guides": g.guide_details,
            }
            for g in gene_results
        ],
    }


def write_gene_summary_tsv(gene_results: list[GeneResult], path: Path) -> None:
    """Write a MAGeCK-compatible gene summary TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "id", "num", "neg|score", "neg|p-value", "neg|fdr", "neg|rank",
            "neg|lfc", "pos|score", "pos|p-value", "pos|fdr", "pos|rank", "pos|lfc",
        ])
        for g in gene_results:
            writer.writerow([
                g.gene, g.n_guides,
                f"{g.neg_score:.6f}", f"{g.neg_p_value:.6f}", f"{g.neg_fdr:.6f}", g.neg_rank,
                f"{g.neg_lfc:.4f}",
                f"{g.pos_score:.6f}", f"{g.pos_p_value:.6f}", f"{g.pos_fdr:.6f}", g.pos_rank,
                f"{g.pos_lfc:.4f}",
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_sample_labels(
    label_str: str, all_samples: list[str],
) -> list[int]:
    """Parse comma-separated sample names into column indices."""
    names = [n.strip() for n in label_str.split(",")]
    indices = []
    for name in names:
        if name not in all_samples:
            raise ValueError(
                f"Sample {name!r} not found. Available: {all_samples}"
            )
        indices.append(all_samples.index(name))
    return indices


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts", type=Path, required=True,
        help="Count table (TSV): sgRNA, gene, sample1, sample2, ...",
    )
    parser.add_argument(
        "--control", required=True,
        help="Comma-separated control sample names",
    )
    parser.add_argument(
        "--treatment", required=True,
        help="Comma-separated treatment sample names",
    )
    parser.add_argument(
        "--json-out", type=Path, default=None,
        help="JSON output path for full results",
    )
    parser.add_argument(
        "--tsv-out", type=Path, default=None,
        help="TSV output path for gene summary (MAGeCK-compatible format)",
    )
    args = parser.parse_args()

    # Load
    samples, guides, genes, counts = load_count_table(args.counts)
    control_idx = parse_sample_labels(args.control, samples)
    treatment_idx = parse_sample_labels(args.treatment, samples)

    control_names = [samples[i] for i in control_idx]
    treatment_names = [samples[i] for i in treatment_idx]

    # Normalize
    normalized, size_factors = median_ratio_normalize(counts)

    # Guide-level stats
    guide_results = compute_guide_stats(
        normalized, guides, genes, control_idx, treatment_idx,
    )

    # Gene-level aggregation
    gene_results = robust_rank_aggregation(guide_results)

    # Output
    payload = format_output(
        gene_results, size_factors, samples, control_names, treatment_names,
    )

    json_text = json.dumps(payload, indent=2)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json_text + "\n", encoding="utf-8")
    else:
        print(json_text)

    if args.tsv_out:
        write_gene_summary_tsv(gene_results, args.tsv_out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
