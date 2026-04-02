#!/usr/bin/env python3
"""Build a deterministic Dash scientific-dashboard starter from toy measurements."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_measurements.tsv"


def load_dash_stack():
    try:
        from dash import Dash, Input, Output, dcc, html
        import plotly.graph_objects as go
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires dash and plotly in slurm/envs/reporting.") from exc
    return Dash, Input, Output, dcc, html, go


def load_measurements(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"Measurement table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Measurement table is empty.")
    required = {"time_h", "signal", "control"}
    if missing := sorted(required - set(rows[0].keys())):
        raise ValueError(f"Measurement table missing required columns: {', '.join(missing)}")
    measurements = []
    for row in rows:
        measurements.append({key: round(float(value), 6) for key, value in row.items()})
    return measurements


def make_metric_options(rows: list[dict[str, float]]) -> list[str]:
    return [key for key in rows[0].keys() if key != "time_h"]


def create_figure(rows: list[dict[str, float]], metric: str):
    _, _, _, _, _, go = load_dash_stack()
    x_values = [row["time_h"] for row in rows]
    y_values = [row[metric] for row in rows]
    mean_value = sum(y_values) / len(y_values)
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=x_values,
            y=y_values,
            mode="lines+markers",
            name=metric,
            line={"width": 3, "color": "#185adb"},
        )
    )
    figure.add_trace(
        go.Scatter(
            x=x_values,
            y=[mean_value] * len(x_values),
            mode="lines",
            name="mean",
            line={"dash": "dash", "width": 2, "color": "#b23a48"},
        )
    )
    figure.update_layout(
        title=f"Toy {metric.title()} Dashboard",
        xaxis_title="Time (h)",
        yaxis_title=metric.title(),
        template="plotly_white",
    )
    return figure


def collect_component_types(component, collected: list[str]) -> None:
    collected.append(type(component).__name__)
    children = getattr(component, "children", None)
    if children is None:
        return
    if isinstance(children, (list, tuple)):
        for child in children:
            if hasattr(child, "to_plotly_json"):
                collect_component_types(child, collected)
    elif hasattr(children, "to_plotly_json"):
        collect_component_types(children, collected)


def build_app(rows: list[dict[str, float]]):
    Dash, Input, Output, dcc, html, _ = load_dash_stack()
    metric_options = make_metric_options(rows)
    default_metric = metric_options[0]
    app = Dash(__name__)
    app.layout = html.Div(
        [
            html.H1("Scientific Dashboard Demo"),
            html.Div(
                [
                    html.Div(f"Rows: {len(rows)}", id="row-count"),
                    html.Div(f"Metrics: {len(metric_options)}", id="metric-count"),
                ],
                id="summary-cards",
            ),
            dcc.Dropdown(
                id="metric",
                options=[{"label": metric.title(), "value": metric} for metric in metric_options],
                value=default_metric,
                clearable=False,
            ),
            dcc.Graph(id="trend"),
        ]
    )

    @app.callback(Output("trend", "figure"), Input("metric", "value"))
    def update_figure(selected_metric: str):
        return create_figure(rows, selected_metric)

    return app, default_metric, metric_options, update_figure


def default_shell_path(html_out: Path) -> Path:
    return html_out.with_name(f"{html_out.stem}_shell.html")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_json(path: Path, payload: dict) -> None:
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Tab-separated measurement table.")
    parser.add_argument("--html-out", type=Path, required=True, help="Standalone preview HTML path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    parser.add_argument("--shell-out", type=Path, default=None, help="Optional Dash shell HTML path.")
    args = parser.parse_args()

    try:
        rows = load_measurements(args.input)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    app, default_metric, metric_options, update_figure = build_app(rows)
    figure = update_figure(default_metric)
    shell_out = args.shell_out or default_shell_path(args.html_out)
    args.html_out.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.html_out, include_plotlyjs=True, full_html=True)
    write_text(shell_out, app.index())

    component_types: list[str] = []
    collect_component_types(app.layout, component_types)
    summary = {
        "measurement_count": len(rows),
        "metric_options": metric_options,
        "default_metric": default_metric,
        "trace_count": len(figure.data),
        "callback_count": len(app.callback_map),
        "layout_component_types": sorted(set(component_types)),
        "html_path": str(args.html_out),
        "html_size_bytes": args.html_out.stat().st_size,
        "shell_path": str(shell_out),
        "shell_size_bytes": shell_out.stat().st_size,
    }
    if args.summary_out is None:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        write_json(args.summary_out, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
