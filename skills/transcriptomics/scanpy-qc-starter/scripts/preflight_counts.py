#!/usr/bin/env python3
"""Compute simple preflight statistics for a toy count matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_matrix(path: Path) -> tuple[list[str], list[list[int]]]:
    rows = path.read_text(encoding="utf-8").strip().splitlines()
    header = rows[0].split("\t")
    matrix: list[list[int]] = []
    for raw_row in rows[1:]:
        fields = raw_row.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"Ragged row in {path}: {raw_row}")
        matrix.append([int(value) for value in fields[1:]])
    return header[1:], matrix


def summarize(matrix: list[list[int]]) -> dict:
    n_rows = len(matrix)
    n_cols = len(matrix[0]) if matrix else 0
    flat = [value for row in matrix for value in row]
    zero_fraction = sum(1 for value in flat if value == 0) / len(flat) if flat else 0.0
    return {
        "genes": n_rows,
        "cells": n_cols,
        "total_counts": sum(flat),
        "zero_fraction": round(zero_fraction, 4),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Tab-separated count matrix")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    _, matrix = parse_matrix(args.input)
    payload = summarize(matrix)
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
