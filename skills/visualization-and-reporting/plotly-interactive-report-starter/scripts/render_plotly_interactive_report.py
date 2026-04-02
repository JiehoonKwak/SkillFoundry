#!/usr/bin/env python3
"""Render a small deterministic interactive Plotly report from toy measurements."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_measurements.tsv"


def load_plotting_stack():
    try:
        import plotly.graph_objects as go
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the statistics prefix with plotly installed at slurm/envs/statistics."
        ) from exc
    return go


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
    ss_res = sum((y - f) ** 2 for y, f in zip(y_values, fitted, strict=True))
    ss_tot = sum((y - mean_y) ** 2 for y in y_values)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    return {
        "slope": round(float(slope), 6),
        "intercept": round(float(intercept), 6),
        "r_squared": round(float(r_squared), 6),
        "fitted": [round(float(value), 6) for value in fitted],
    }


def render_report(x_values: list[float], y_values: list[float], fit: dict, html_out: Path) -> None:
    go = load_plotting_stack()
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            mode="markers",
            marker={"size": 10, "color": "#1b4965"},
            name="observed",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=x_values,
            y=fit["fitted"],
            mode="lines",
            line={"width": 3, "color": "#ca6702"},
            name="fitted",
        )
    )
    figure.update_layout(
        title="Toy Interactive Trend Report",
        xaxis_title="x",
        yaxis_title="y",
        template="plotly_white",
    )
    html_out.parent.mkdir(parents=True, exist_ok=True)
    # Keep the generated report self-contained so committed assets remain viewable offline.
    figure.write_html(html_out, include_plotlyjs=True, full_html=True)


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
    parser.add_argument("--html-out", type=Path, required=True, help="Output HTML report path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    x_values, y_values = load_points(args.input)
    fit = fit_line(x_values, y_values)
    render_report(x_values, y_values, fit, args.html_out)
    summary = {
        "point_count": len(x_values),
        "trace_count": 2,
        "title": "Toy Interactive Trend Report",
        "html_path": str(args.html_out),
        "html_size_bytes": args.html_out.stat().st_size,
        **fit,
    }
    write_summary(summary, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
