#!/usr/bin/env python3
"""Validate a small tabular file against a Frictionless schema."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from frictionless import Resource, Schema


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "toy_people_valid.csv"
DEFAULT_SCHEMA = ROOT / "examples" / "toy_people_schema.json"


def load_schema(path: Path) -> Schema:
    return Schema.from_descriptor(json.loads(path.read_text(encoding="utf-8")))


def read_table_stats(path: Path) -> tuple[list[str], int]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        field_names = list(reader.fieldnames or [])
        row_count = sum(1 for _ in reader)
    return field_names, row_count


def build_summary(input_path: Path, schema_path: Path) -> dict[str, object]:
    schema = load_schema(schema_path)
    resource = Resource(path=input_path.name, basepath=str(input_path.parent), schema=schema)
    report = resource.validate()
    field_names, row_count = read_table_stats(input_path)
    errors = [
        {
            "row_number": row_number,
            "field_number": field_number,
            "type": error_type,
            "note": note,
        }
        for row_number, field_number, error_type, note in report.flatten(["rowNumber", "fieldNumber", "type", "note"])
    ]
    return {
        "input_path": str(input_path),
        "schema_path": str(schema_path),
        "valid": bool(report.valid),
        "field_names": field_names,
        "row_count": row_count,
        "error_count": int(report.stats["errors"]),
        "warning_count": int(report.stats["warnings"]),
        "errors": errors,
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="CSV or TSV file to validate.")
    parser.add_argument("--schema", type=Path, default=DEFAULT_SCHEMA, help="Frictionless schema descriptor.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input table not found: {args.input}")
    if not args.schema.exists():
        raise SystemExit(f"Schema descriptor not found: {args.schema}")

    payload = build_summary(args.input, args.schema)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
