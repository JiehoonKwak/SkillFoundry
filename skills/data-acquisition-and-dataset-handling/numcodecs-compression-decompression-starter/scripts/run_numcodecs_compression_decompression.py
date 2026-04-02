#!/usr/bin/env python3
"""Round-trip a tiny integer matrix through numcodecs and summarize the result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numcodecs import Blosc


CODEC = Blosc(cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE)


def load_matrix(path: Path) -> np.ndarray:
    rows: list[list[int]] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        rows.append([int(item) for item in line.split("\t")])
    matrix = np.array(rows, dtype=np.int32)
    if matrix.ndim != 2:
        raise ValueError("Expected a 2D integer matrix.")
    return matrix


def summarize_roundtrip(matrix: np.ndarray, encoded: bytes, decoded: np.ndarray, input_path: Path) -> dict:
    return {
        "input_path": str(input_path.resolve()),
        "shape": list(matrix.shape),
        "dtype": str(matrix.dtype),
        "codec_name": "blosc-zstd",
        "encoded_nbytes": len(encoded),
        "input_nbytes": int(matrix.nbytes),
        "compression_ratio": round(len(encoded) / float(matrix.nbytes), 6),
        "decoded_equal": bool(np.array_equal(matrix, decoded)),
        "row_sums": [int(value) for value in np.sum(matrix, axis=1)],
        "column_maxima": [int(value) for value in np.max(matrix, axis=0)],
    }


def run_roundtrip(input_path: Path) -> dict:
    matrix = load_matrix(input_path)
    encoded = CODEC.encode(matrix.tobytes())
    decoded = np.frombuffer(CODEC.decode(encoded), dtype=matrix.dtype).reshape(matrix.shape)
    return summarize_roundtrip(matrix, encoded, decoded, input_path)


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Input TSV matrix.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary output.")
    args = parser.parse_args()

    payload = run_roundtrip(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
