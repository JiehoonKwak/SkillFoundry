#!/usr/bin/env python3
"""Harmonize tiny metadata tables into a canonical TSV and JSON summary."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def load_mapping(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No rows were found in {path}.")
    return rows


def normalize_value(field: str, value: str, mapping: dict) -> str | int:
    value = value.strip()
    if field == "age_years":
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"Invalid integer age value: {value!r}") from exc
    normalizers = mapping.get("normalizers", {})
    field_map = normalizers.get(field, {})
    return field_map.get(value.lower(), value)


def harmonize_rows(input_paths: list[Path], mapping: dict) -> list[dict[str, str | int]]:
    canonical_fields = mapping["canonical_fields"]
    source_specs = mapping["sources"]
    harmonized: list[dict[str, str | int]] = []
    for input_path in input_paths:
        source_name = input_path.name
        if source_name not in source_specs:
            raise KeyError(f"No source mapping was provided for {source_name}.")
        source_map = source_specs[source_name]
        rows = read_rows(input_path)
        for row in rows:
            normalized: dict[str, str | int] = {}
            for canonical in canonical_fields:
                if canonical == "source_file":
                    normalized[canonical] = source_name
                    continue
                source_field = source_map.get(canonical)
                if source_field is None:
                    raise KeyError(f"Missing source field mapping for {canonical} in {source_name}.")
                if source_field not in row:
                    raise KeyError(f"{source_name} is missing mapped field {source_field}.")
                normalized[canonical] = normalize_value(canonical, row[source_field], mapping)
            harmonized.append(normalized)
    harmonized.sort(key=lambda item: str(item["sample_id"]))
    return harmonized


def write_tsv(rows: list[dict[str, str | int]], out_path: Path, fieldnames: list[str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_summary(rows: list[dict[str, str | int]], out_path: Path) -> dict:
    source_counts = Counter(str(row["source_file"]) for row in rows)
    return {
        "harmonized_tsv": str(out_path.resolve()),
        "row_count": len(rows),
        "field_count": len(rows[0]) if rows else 0,
        "fields": list(rows[0].keys()) if rows else [],
        "source_counts": dict(sorted(source_counts.items())),
        "distinct_conditions": sorted({str(row["condition"]) for row in rows}),
        "distinct_sex": sorted({str(row["sex"]) for row in rows}),
        "first_row": rows[0] if rows else None,
    }


def write_json(payload: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True, help="Input TSV file. Repeat for multiple sources.")
    parser.add_argument("--mapping", type=Path, required=True, help="JSON mapping from source columns to canonical fields.")
    parser.add_argument("--out-tsv", type=Path, required=True, help="Output harmonized TSV.")
    parser.add_argument("--summary-out", type=Path, required=True, help="Output JSON summary path.")
    args = parser.parse_args()

    mapping = load_mapping(args.mapping)
    rows = harmonize_rows(args.input, mapping)
    write_tsv(rows, args.out_tsv, mapping["canonical_fields"])
    summary = build_summary(rows, args.out_tsv)
    write_json(summary, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
