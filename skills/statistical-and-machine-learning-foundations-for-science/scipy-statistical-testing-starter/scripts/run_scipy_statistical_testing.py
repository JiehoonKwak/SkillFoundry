#!/usr/bin/env python3
"""Run deterministic statistical tests over two toy groups."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from scipy import stats


def load_groups(path: Path) -> dict[str, list[float]]:
    groups: dict[str, list[float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            groups.setdefault(row["group"], []).append(float(row["value"]))
    if len(groups) != 2:
        raise ValueError("Expected exactly two groups in the input.")
    return groups


def cohen_d(a: list[float], b: list[float]) -> float:
    mean_a = sum(a) / len(a)
    mean_b = sum(b) / len(b)
    var_a = sum((value - mean_a) ** 2 for value in a) / (len(a) - 1)
    var_b = sum((value - mean_b) ** 2 for value in b) / (len(b) - 1)
    pooled = (((len(a) - 1) * var_a) + ((len(b) - 1) * var_b)) / (len(a) + len(b) - 2)
    return (mean_b - mean_a) / (pooled ** 0.5)


def run_tests(input_path: Path) -> dict:
    groups = load_groups(input_path)
    group_names = sorted(groups)
    a, b = groups[group_names[0]], groups[group_names[1]]
    ttest = stats.ttest_ind(a, b, equal_var=False)
    mann = stats.mannwhitneyu(a, b, alternative="two-sided")
    mean_a = sum(a) / len(a)
    mean_b = sum(b) / len(b)
    return {
        "input_path": str(input_path.resolve()),
        "groups": group_names,
        "group_sizes": {group_names[0]: len(a), group_names[1]: len(b)},
        "group_means": {group_names[0]: round(mean_a, 6), group_names[1]: round(mean_b, 6)},
        "higher_mean_group": group_names[0] if mean_a > mean_b else group_names[1],
        "welch_ttest_pvalue": round(float(ttest.pvalue), 8),
        "mannwhitney_pvalue": round(float(mann.pvalue), 8),
        "cohen_d": round(float(cohen_d(a, b)), 6),
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_tests(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
