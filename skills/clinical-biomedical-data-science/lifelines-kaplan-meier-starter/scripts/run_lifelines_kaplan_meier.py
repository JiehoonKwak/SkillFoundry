#!/usr/bin/env python3
"""Fit deterministic Kaplan-Meier curves for a toy cohort."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_survival_cohort.tsv"


def load_stack():
    try:
        import pandas as pd
        from lifelines import KaplanMeierFitter
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires the statistics prefix at slurm/envs/statistics.") from exc
    return pd, KaplanMeierFitter, plt


def load_table(path: Path):
    pd, _, _ = load_stack()
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"No cohort rows found in {path}")
    frame = pd.DataFrame(rows)
    frame["time"] = frame["time"].astype(float)
    frame["event"] = frame["event"].astype(int)
    return frame


def run_kaplan_meier(input_path: Path, png_out: Path | None) -> dict:
    _, KaplanMeierFitter, plt = load_stack()
    cohort = load_table(input_path)
    groups = sorted(cohort["group"].unique().tolist())
    times = [0, 3, 4, 6]
    medians = {}
    survival_at_times = {}

    if png_out is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = ax = None

    for group in groups:
        subset = cohort[cohort["group"] == group]
        kmf = KaplanMeierFitter()
        kmf.fit(subset["time"], subset["event"], label=group)
        medians[group] = round(float(kmf.median_survival_time_), 6)
        survival_at_times[group] = {str(time): round(float(kmf.predict(time)), 6) for time in times}
        if ax is not None:
            kmf.plot(ax=ax, ci_show=False)

    if ax is not None and fig is not None:
        ax.set_title("Toy Kaplan-Meier Survival Curves")
        ax.set_xlabel("Time")
        ax.set_ylabel("Survival probability")
        ax.grid(alpha=0.2)
        png_out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(png_out, dpi=150)
        plt.close(fig)

    return {
        "input_path": str(input_path),
        "sample_count": int(len(cohort)),
        "event_rate": round(float(cohort["event"].mean()), 6),
        "groups": groups,
        "median_survival_by_group": medians,
        "survival_at_times": survival_at_times,
        "time_grid": times,
        "png_output": str(png_out) if png_out is not None else None,
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Cohort TSV path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--png-out", type=Path, default=None, help="Optional Kaplan-Meier PNG path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input cohort not found: {args.input}")

    payload = run_kaplan_meier(input_path=args.input, png_out=args.png_out)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
