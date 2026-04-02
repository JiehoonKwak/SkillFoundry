#!/usr/bin/env python3
"""Run a deterministic fairness audit over a toy clinical cohort."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_fairness_cohort.tsv"
DEFAULT_LABEL_COLUMN = "readmission_30d"
DEFAULT_PREDICTION_COLUMN = "predicted_high_risk"
DEFAULT_GROUP_COLUMN = "sensitive_group"


def load_stack():
    try:
        import pandas as pd
        from fairlearn.metrics import (
            MetricFrame,
            demographic_parity_difference,
            demographic_parity_ratio,
            equalized_odds_difference,
            equalized_odds_ratio,
            false_positive_rate,
            selection_rate,
            true_positive_rate,
        )
        from sklearn.metrics import accuracy_score
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires the statistics prefix at slurm/envs/statistics.") from exc
    return (
        pd,
        MetricFrame,
        demographic_parity_difference,
        demographic_parity_ratio,
        equalized_odds_difference,
        equalized_odds_ratio,
        false_positive_rate,
        selection_rate,
        true_positive_rate,
        accuracy_score,
    )


def format_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def round_float(value: float) -> float:
    return round(float(value), 6)


def load_table(path: Path):
    pd, *_ = load_stack()
    if not path.exists():
        raise SystemExit(f"Input cohort not found: {path}")
    frame = pd.read_csv(path, sep="\t")
    if frame.empty:
        raise SystemExit(f"No cohort rows found in {path}")
    return frame


def validate_columns(frame, label_column: str, prediction_column: str, group_column: str) -> None:
    required = [label_column, prediction_column, group_column]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        missing_text = ", ".join(missing)
        raise SystemExit(f"Missing required columns: {missing_text}")


def validate_binary_column(frame, column: str) -> None:
    pd, *_ = load_stack()
    try:
        numeric = pd.to_numeric(frame[column], errors="raise")
    except Exception as exc:
        raise SystemExit(f"Column '{column}' must contain only binary 0/1 values.") from exc
    values = set(numeric.dropna().tolist())
    if not values.issubset({0, 1, 0.0, 1.0}):
        raise SystemExit(f"Column '{column}' must contain only binary 0/1 values.")
    frame[column] = numeric.astype(int)


def build_flags(
    fairness_metrics: dict,
    demographic_parity_threshold: float,
    equalized_odds_threshold: float,
    minimum_ratio_threshold: float,
) -> list[str]:
    flags: list[str] = []
    if fairness_metrics["demographic_parity_difference"] > demographic_parity_threshold:
        flags.append(
            f"demographic_parity_difference>{demographic_parity_threshold:.3f}"
        )
    if fairness_metrics["equalized_odds_difference"] > equalized_odds_threshold:
        flags.append(f"equalized_odds_difference>{equalized_odds_threshold:.3f}")
    if fairness_metrics["demographic_parity_ratio"] < minimum_ratio_threshold:
        flags.append(f"demographic_parity_ratio<{minimum_ratio_threshold:.3f}")
    if fairness_metrics["equalized_odds_ratio"] < minimum_ratio_threshold:
        flags.append(f"equalized_odds_ratio<{minimum_ratio_threshold:.3f}")
    return flags


def run_audit(
    input_path: Path,
    label_column: str,
    prediction_column: str,
    group_column: str,
    demographic_parity_threshold: float,
    equalized_odds_threshold: float,
    minimum_ratio_threshold: float,
) -> dict:
    (
        _,
        MetricFrame,
        demographic_parity_difference,
        demographic_parity_ratio,
        equalized_odds_difference,
        equalized_odds_ratio,
        false_positive_rate,
        selection_rate,
        true_positive_rate,
        accuracy_score,
    ) = load_stack()

    frame = load_table(input_path)
    validate_columns(frame, label_column, prediction_column, group_column)
    validate_binary_column(frame, label_column)
    validate_binary_column(frame, prediction_column)
    if frame[group_column].isna().any():
        raise SystemExit(f"Column '{group_column}' must not contain missing values.")
    frame[group_column] = frame[group_column].astype(str)
    if (frame[group_column].str.strip() == "").any():
        raise SystemExit(f"Column '{group_column}' must not contain blank values.")

    metrics = {
        "accuracy": accuracy_score,
        "selection_rate": selection_rate,
        "true_positive_rate": true_positive_rate,
        "false_positive_rate": false_positive_rate,
    }
    metric_frame = MetricFrame(
        metrics=metrics,
        y_true=frame[label_column],
        y_pred=frame[prediction_column],
        sensitive_features=frame[group_column],
    )
    by_group = metric_frame.by_group
    groups = sorted(str(group) for group in by_group.index.tolist())
    metrics_by_group = {
        group: {
            metric_name: round_float(by_group.loc[group, metric_name])
            for metric_name in metrics
        }
        for group in groups
    }

    fairness_metrics = {
        "demographic_parity_difference": round_float(
            demographic_parity_difference(
                frame[label_column],
                frame[prediction_column],
                sensitive_features=frame[group_column],
            )
        ),
        "demographic_parity_ratio": round_float(
            demographic_parity_ratio(
                frame[label_column],
                frame[prediction_column],
                sensitive_features=frame[group_column],
            )
        ),
        "equalized_odds_difference": round_float(
            equalized_odds_difference(
                frame[label_column],
                frame[prediction_column],
                sensitive_features=frame[group_column],
            )
        ),
        "equalized_odds_ratio": round_float(
            equalized_odds_ratio(
                frame[label_column],
                frame[prediction_column],
                sensitive_features=frame[group_column],
            )
        ),
    }

    selection_rates = {group: values["selection_rate"] for group, values in metrics_by_group.items()}
    group_counts = {
        group: int((frame[group_column] == group).sum())
        for group in groups
    }
    flags = build_flags(
        fairness_metrics=fairness_metrics,
        demographic_parity_threshold=demographic_parity_threshold,
        equalized_odds_threshold=equalized_odds_threshold,
        minimum_ratio_threshold=minimum_ratio_threshold,
    )

    return {
        "input_path": format_path(input_path),
        "sample_count": int(len(frame)),
        "columns": {
            "label": label_column,
            "prediction": prediction_column,
            "sensitive_group": group_column,
        },
        "group_order": groups,
        "group_counts": group_counts,
        "overall_accuracy": round_float(accuracy_score(frame[label_column], frame[prediction_column])),
        "overall_positive_label_rate": round_float(frame[label_column].mean()),
        "overall_selection_rate": round_float(frame[prediction_column].mean()),
        "metrics_by_group": metrics_by_group,
        "selection_rate_extremes": {
            "highest_group": max(selection_rates, key=selection_rates.get),
            "highest_value": selection_rates[max(selection_rates, key=selection_rates.get)],
            "lowest_group": min(selection_rates, key=selection_rates.get),
            "lowest_value": selection_rates[min(selection_rates, key=selection_rates.get)],
        },
        "fairness_metrics": fairness_metrics,
        "thresholds": {
            "demographic_parity_difference": round_float(demographic_parity_threshold),
            "equalized_odds_difference": round_float(equalized_odds_threshold),
            "minimum_ratio": round_float(minimum_ratio_threshold),
        },
        "flags": flags,
        "interpretation": (
            "Selection and error-rate gaps remain visible between groups; "
            "inspect thresholding, label quality, and cohort balance before deployment."
            if flags
            else "No configured fairness thresholds were violated for this audit."
        ),
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="TSV cohort input path.")
    parser.add_argument("--label-column", default=DEFAULT_LABEL_COLUMN, help="Binary label column.")
    parser.add_argument("--prediction-column", default=DEFAULT_PREDICTION_COLUMN, help="Binary prediction column.")
    parser.add_argument("--group-column", default=DEFAULT_GROUP_COLUMN, help="Sensitive group column.")
    parser.add_argument(
        "--demographic-parity-threshold",
        type=float,
        default=0.1,
        help="Flag when demographic parity difference exceeds this value.",
    )
    parser.add_argument(
        "--equalized-odds-threshold",
        type=float,
        default=0.1,
        help="Flag when equalized odds difference exceeds this value.",
    )
    parser.add_argument(
        "--minimum-ratio-threshold",
        type=float,
        default=0.8,
        help="Flag when fairness ratios fall below this value.",
    )
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    payload = run_audit(
        input_path=args.input,
        label_column=args.label_column,
        prediction_column=args.prediction_column,
        group_column=args.group_column,
        demographic_parity_threshold=args.demographic_parity_threshold,
        equalized_odds_threshold=args.equalized_odds_threshold,
        minimum_ratio_threshold=args.minimum_ratio_threshold,
    )
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
