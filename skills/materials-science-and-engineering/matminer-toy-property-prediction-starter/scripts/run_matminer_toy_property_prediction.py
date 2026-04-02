#!/usr/bin/env python3
"""Fit a tiny deterministic property predictor on toy material compositions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_materials.tsv"


def load_modeling_stack():
    try:
        import pandas as pd
        from matminer.featurizers.conversions import StrToComposition
        from matminer.featurizers.composition import ElementProperty
        from sklearn.tree import DecisionTreeRegressor
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated materials prefix at slurm/envs/materials. "
            "Run it with slurm/envs/materials/bin/python."
        ) from exc
    return pd, StrToComposition, ElementProperty, DecisionTreeRegressor


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def rank_feature_importances(feature_names: list[str], importances: list[float], limit: int = 5) -> list[dict[str, float]]:
    pairs = sorted(zip(feature_names, importances, strict=True), key=lambda item: (-item[1], item[0]))
    return [
        {"feature": name, "importance": round(float(value), 6)}
        for name, value in pairs[:limit]
        if float(value) > 0.0
    ]


def train_predictor(path: Path) -> dict:
    pd, StrToComposition, ElementProperty, DecisionTreeRegressor = load_modeling_stack()
    dataframe = pd.read_csv(path, sep="\t")
    dataframe = StrToComposition(target_col_id="composition").featurize_dataframe(dataframe, "formula")
    featurizer = ElementProperty.from_preset("magpie")
    dataframe = featurizer.featurize_dataframe(dataframe, "composition", ignore_errors=False)
    feature_columns = [column for column in dataframe.columns if column not in {"formula", "toy_property", "composition"}]
    model = DecisionTreeRegressor(random_state=20260312)
    model.fit(dataframe[feature_columns], dataframe["toy_property"])
    predictions = model.predict(dataframe[feature_columns])
    residuals = [abs(float(prediction) - float(target)) for prediction, target in zip(predictions, dataframe["toy_property"], strict=True)]
    return {
        "sample_count": int(len(dataframe)),
        "feature_count": int(len(feature_columns)),
        "training_mae": round(sum(residuals) / len(residuals), 6),
        "predictions": [
            {
                "formula": str(formula),
                "target": round(float(target), 6),
                "predicted": round(float(predicted), 6),
            }
            for formula, target, predicted in zip(
                dataframe["formula"],
                dataframe["toy_property"],
                predictions,
                strict=True,
            )
        ],
        "top_feature_importances": rank_feature_importances(
            feature_columns,
            list(model.feature_importances_),
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Tab-separated formula/property table.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = train_predictor(args.input)
    payload["input_path"] = str(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
