#!/usr/bin/env python3
"""Write a tiny numeric matrix to a chunked Zarr store and summarize the result."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy
from numcodecs import Blosc


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_matrix.tsv"
DEFAULT_CHUNKS = (2, 2)


def read_matrix(path: Path) -> numpy.ndarray:
    rows: list[list[float]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            rows.append([float(value) for value in row])
    if not rows:
        raise ValueError("Input matrix is empty.")
    width = len(rows[0])
    if width == 0:
        raise ValueError("Input matrix must contain at least one column.")
    for index, row in enumerate(rows, start=1):
        if len(row) != width:
            raise ValueError(f"Row {index} has {len(row)} columns; expected {width}.")
    return numpy.asarray(rows, dtype="float32")


def build_summary(matrix: numpy.ndarray, store_out: Path, chunks: tuple[int, int]) -> dict:
    compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    if store_out.exists():
        shutil.rmtree(store_out)
    store_out.parent.mkdir(parents=True, exist_ok=True)
    store_out.mkdir(parents=True, exist_ok=True)
    metadata = {
        "chunks": [int(chunks[0]), int(chunks[1])],
        "compressor": compressor.get_config(),
        "dtype": "<f4",
        "fill_value": 0.0,
        "filters": None,
        "order": "C",
        "shape": [int(dimension) for dimension in matrix.shape],
        "zarr_format": 2,
    }
    (store_out / ".zarray").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (store_out / ".zattrs").write_text("{}\n", encoding="utf-8")
    chunk_file_count = 0
    row_step, col_step = chunks
    for row_offset in range(0, matrix.shape[0], row_step):
        for col_offset in range(0, matrix.shape[1], col_step):
            chunk = matrix[row_offset : row_offset + row_step, col_offset : col_offset + col_step]
            if chunk.shape != chunks:
                padded = numpy.zeros(chunks, dtype=matrix.dtype)
                padded[: chunk.shape[0], : chunk.shape[1]] = chunk
                chunk = padded
            encoded = compressor.encode(numpy.ascontiguousarray(chunk).tobytes())
            (store_out / f"{row_offset // row_step}.{col_offset // col_step}").write_bytes(encoded)
            chunk_file_count += 1
    first_chunk = matrix[: row_step, : col_step]
    return {
        "input_path": "",
        "store_path": str(store_out),
        "shape": [int(dimension) for dimension in matrix.shape],
        "chunk_shape": [int(dimension) for dimension in chunks],
        "dtype": str(matrix.dtype),
        "zarr_format": 2,
        "compressor": {
            "codec": "blosc",
            "cname": "zstd",
            "clevel": 3,
            "shuffle": "bitshuffle",
        },
        "matrix_sum": round(float(matrix.sum()), 6),
        "row_means": [round(float(value), 6) for value in matrix.mean(axis=1)],
        "first_chunk_sum": round(float(first_chunk.sum()), 6),
        "stored_entry_count": int(matrix.size),
        "store_file_count": chunk_file_count + 2,
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Tab-delimited numeric matrix.")
    parser.add_argument("--store-out", type=Path, required=True, help="Output Zarr store directory.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--chunk-rows", type=int, default=DEFAULT_CHUNKS[0], help="Chunk height.")
    parser.add_argument("--chunk-cols", type=int, default=DEFAULT_CHUNKS[1], help="Chunk width.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input matrix not found: {args.input}")
    if args.chunk_rows < 1 or args.chunk_cols < 1:
        raise SystemExit("Chunk dimensions must be positive.")

    matrix = read_matrix(args.input)
    payload = build_summary(matrix, args.store_out, (args.chunk_rows, args.chunk_cols))
    payload["input_path"] = str(args.input)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
