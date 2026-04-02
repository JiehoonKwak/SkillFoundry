#!/usr/bin/env python3
"""Convert a tiny tabular file to Parquet with PyArrow and summarize the round trip."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def coerce_value(value: str) -> object:
    value = value.strip()
    if not value:
        return value
    if value.isdigit():
        return int(value)
    try:
        return float(value)
    except ValueError:
        return value


def load_rows(path: Path) -> list[dict[str, object]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [{key: coerce_value(value) for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError("No tabular rows were found in the input.")
    return rows


def summarize_table(input_path: Path, parquet_out: Path, table: pa.Table, roundtrip_table: pa.Table) -> dict:
    return {
        "input_path": str(input_path.resolve()),
        "parquet_path": str(parquet_out.resolve()),
        "row_count": int(table.num_rows),
        "column_count": int(table.num_columns),
        "columns": [field.name for field in table.schema],
        "schema": [{"name": field.name, "type": str(field.type)} for field in table.schema],
        "roundtrip_rows_equal": bool(table.to_pylist() == roundtrip_table.to_pylist()),
        "first_row": roundtrip_table.slice(0, 1).to_pylist()[0],
    }


def run_conversion(input_path: Path, parquet_out: Path) -> dict:
    rows = load_rows(input_path)
    table = pa.Table.from_pylist(rows)
    parquet_out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, parquet_out)
    roundtrip_table = pq.read_table(parquet_out)
    return summarize_table(input_path, parquet_out, table, roundtrip_table)


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Input TSV table.")
    parser.add_argument("--parquet-out", type=Path, required=True, help="Output Parquet file path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    payload = run_conversion(args.input, args.parquet_out)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
