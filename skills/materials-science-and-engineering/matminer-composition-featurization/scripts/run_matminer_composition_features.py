#!/usr/bin/env python3
"""Compute compact matminer stoichiometry features for one or more formulas."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any


DEFAULT_FORMULAS = ["Fe2O3", "LiFePO4", "SiO2"]


def load_toolkit():
    try:
        from matminer.featurizers.composition import Stoichiometry
        from pymatgen.core import Composition
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires the materials prefix at slurm/envs/materials.") from exc
    return Stoichiometry, Composition


def normalize_formulas(formulas: list[str] | None) -> list[str]:
    values = DEFAULT_FORMULAS if not formulas else formulas
    cleaned = [formula.strip() for formula in values if formula.strip()]
    if not cleaned:
        raise ValueError("At least one non-empty formula is required.")
    return cleaned


def featurize_formulas(formulas: list[str]) -> dict[str, Any]:
    Stoichiometry, Composition = load_toolkit()
    featurizer = Stoichiometry()
    labels = featurizer.feature_labels()
    summaries = []
    for formula in formulas:
        composition = Composition(formula)
        values = featurizer.featurize(composition)
        amount_dict = composition.get_el_amt_dict()
        total = float(sum(amount_dict.values()))
        top_fractions = sorted(
            (
                {
                    "element": element,
                    "fraction": round(float(amount / total), 6),
                }
                for element, amount in amount_dict.items()
            ),
            key=lambda item: (-item["fraction"], item["element"]),
        )
        summaries.append(
            {
                "input_formula": formula,
                "reduced_formula": composition.reduced_formula,
                "stoichiometry_features": {
                    label: round(float(value), 6)
                    for label, value in zip(labels, values, strict=True)
                },
                "top_element_fractions": top_fractions[:3],
            }
        )
    return {
        "matminer_version": version("matminer"),
        "pymatgen_version": version("pymatgen"),
        "formula_count": len(summaries),
        "summaries": summaries,
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--formula",
        action="append",
        default=None,
        help="Formula to featurize. Repeat for multiple formulas. Defaults to a built-in example set.",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        formulas = normalize_formulas(args.formula)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    payload = featurize_formulas(formulas)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
