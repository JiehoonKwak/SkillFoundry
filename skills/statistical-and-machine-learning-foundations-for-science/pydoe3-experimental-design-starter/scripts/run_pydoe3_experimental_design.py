#!/usr/bin/env python3
"""Generate a tiny full-factorial experimental design with pyDOE3."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyDOE3 import ff2n


def load_factor_spec(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    factor_names = payload["factor_names"]
    low = payload["low"]
    high = payload["high"]
    if not (len(factor_names) == len(low) == len(high)):
        raise ValueError("factor_names, low, and high must have the same length.")
    return payload


def run_design(input_path: Path) -> dict:
    spec = load_factor_spec(input_path)
    factor_names = spec["factor_names"]
    design = ff2n(len(factor_names))
    rows = []
    for coded_row in design:
        mapped = {}
        for index, factor_name in enumerate(factor_names):
            mapped[factor_name] = float(spec["low"][index] if coded_row[index] < 0 else spec["high"][index])
        rows.append(mapped)
    return {
        "input_path": str(input_path.resolve()),
        "factor_names": factor_names,
        "row_count": len(rows),
        "rows": rows,
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
    payload = run_design(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
