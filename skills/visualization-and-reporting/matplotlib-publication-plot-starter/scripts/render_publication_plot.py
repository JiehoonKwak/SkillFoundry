#!/usr/bin/env python3
"""Render a deterministic publication-style Matplotlib figure from toy measurements."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_measurements.tsv"


def load_plotting_stack():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated statistics prefix at slurm/envs/statistics. "
            "Run it with slurm/envs/statistics/bin/python."
        ) from exc
    return plt, np


def load_points(path: Path) -> tuple[list[float], list[float]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    x_values = [float(row["x"]) for row in rows]
    y_values = [float(row["y"]) for row in rows]
    return x_values, y_values


def fit_line(x_values: list[float], y_values: list[float]) -> dict:
    n = len(x_values)
    if n != len(y_values) or n < 2:
        raise ValueError("Need at least two paired x/y observations.")
    mean_x = sum(x_values) / n
    mean_y = sum(y_values) / n
    ss_xx = sum((x - mean_x) ** 2 for x in x_values)
    ss_xy = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_values, y_values, strict=True))
    if ss_xx == 0:
        raise ValueError("x values must not all be identical.")
    slope = ss_xy / ss_xx
    intercept = mean_y - slope * mean_x
    fitted = [intercept + slope * x for x in x_values]
    residuals = [y - f for y, f in zip(y_values, fitted, strict=True)]
    ss_res = sum((y - f) ** 2 for y, f in zip(y_values, fitted, strict=True))
    ss_tot = sum((y - mean_y) ** 2 for y in y_values)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    return {
        "slope": round(float(slope), 6),
        "intercept": round(float(intercept), 6),
        "r_squared": round(float(r_squared), 6),
        "residual_mae": round(float(sum(abs(value) for value in residuals) / n), 6),
        "fitted": [round(float(value), 6) for value in fitted],
        "residuals": [round(float(value), 6) for value in residuals],
    }


def format_x_tick(value: float) -> str:
    return f"{value:g}"


def render_figure(x_values: list[float], y_values: list[float], fit: dict, png_out: Path) -> None:
    plt, np = load_plotting_stack()
    x_array = np.asarray(x_values, dtype=float)
    y_array = np.asarray(y_values, dtype=float)
    fitted = np.asarray(fit["fitted"], dtype=float)
    residuals = np.asarray(fit["residuals"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    axes[0].scatter(x_array, y_array, color="#1b4965", label="observed", s=42)
    axes[0].plot(x_array, fitted, color="#ca6702", linewidth=2.2, label="fitted")
    axes[0].set_title("Trend")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].legend(frameon=False)

    axes[1].axhline(0.0, color="#6c757d", linewidth=1.2)
    axes[1].bar([format_x_tick(float(x)) for x in x_array], residuals, color="#5b8e7d")
    axes[1].set_title("Residuals")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("observed - fitted")

    png_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_out, dpi=180)
    plt.close(fig)


def write_summary(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Tab-separated x/y measurements.")
    parser.add_argument("--png-out", type=Path, required=True, help="Output PNG path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    x_values, y_values = load_points(args.input)
    fit = fit_line(x_values, y_values)
    render_figure(x_values, y_values, fit, args.png_out)
    summary = {
        "point_count": len(x_values),
        "x_range": [min(x_values), max(x_values)],
        "y_range": [min(y_values), max(y_values)],
        "figure_path": str(args.png_out),
        **fit,
    }
    write_summary(summary, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
