#!/usr/bin/env python3
"""Featurize formulas with Matminer and fit a deterministic toy regression model."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path


TOY_ROWS = [
    {"formula": "Fe2O3", "target": 1.2},
    {"formula": "LiFePO4", "target": 3.4},
    {"formula": "SiO2", "target": 6.1},
    {"formula": "NaCl", "target": 8.5},
    {"formula": "MgO", "target": 7.8},
    {"formula": "Al2O3", "target": 5.9},
]


def load_stack():
    warnings.filterwarnings("ignore", message="MagpieData\\(impute_nan=False\\).*", category=UserWarning)
    try:
        import numpy as np
        from matminer.featurizers.composition.composite import ElementProperty
        from pymatgen.core import Composition
        from sklearn.linear_model import Ridge
        from sklearn.metrics import r2_score
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the dedicated materials prefix at slurm/envs/materials. "
            "Run it with slurm/envs/materials/bin/python."
        ) from exc
    return np, ElementProperty, Composition, Ridge, r2_score


def featurize_rows(rows: list[dict]) -> tuple[list[str], list[list[float]], list[float], list[str]]:
    _, ElementProperty, Composition, _, _ = load_stack()
    featurizer = ElementProperty.from_preset("magpie")
    labels = featurizer.feature_labels()
    features = []
    targets = []
    formulas = []
    for row in rows:
        formulas.append(row["formula"])
        targets.append(float(row["target"]))
        features.append(featurizer.featurize(Composition(row["formula"])))
    return labels, features, targets, formulas


def run_regression(rows: list[dict]) -> dict:
    np, _, _, Ridge, r2_score = load_stack()
    labels, features, targets, formulas = featurize_rows(rows)
    x_array = np.asarray(features, dtype=float)
    y_array = np.asarray(targets, dtype=float)
    train_x = x_array[:4]
    train_y = y_array[:4]
    test_x = x_array[4:]
    test_y = y_array[4:]
    model = Ridge(alpha=1.0)
    model.fit(train_x, train_y)
    train_predictions = model.predict(train_x)
    holdout_predictions = model.predict(test_x)
    return {
        "feature_count": len(labels),
        "train_count": int(train_x.shape[0]),
        "test_count": int(test_x.shape[0]),
        "train_r2": round(float(r2_score(train_y, train_predictions)), 6),
        "holdout_predictions": [
            {
                "formula": formulas[idx + 4],
                "target": round(float(test_y[idx]), 6),
                "predicted": round(float(holdout_predictions[idx]), 6),
            }
            for idx in range(len(holdout_predictions))
        ],
        "top_feature_names": labels[:5],
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = run_regression(TOY_ROWS)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
