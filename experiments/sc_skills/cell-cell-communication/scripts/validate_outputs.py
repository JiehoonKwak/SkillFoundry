#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


SKILL_DIR = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_markdown_sections(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def validate_tsv(path: Path, required_columns: list[str]) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AssertionError(f"Missing TSV header in {path}")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise AssertionError(f"Missing TSV columns in {path}: {missing}")
        first_row = next(reader, None)
        if first_row is None:
            raise AssertionError(f"Empty TSV deliverable: {path}")


def validate_json(path: Path, required_keys: list[str]) -> None:
    payload = load_json(path)
    missing = [key for key in required_keys if key not in payload]
    if missing:
        raise AssertionError(f"Missing JSON keys in {path}: {missing}")
    if not payload["written_files"]:
        raise AssertionError(f"written_files must not be empty in {path}")


def validate_outputs(skill_dir: Path, outdir: Path) -> None:
    metadata = load_json(skill_dir / "metadata.yaml")
    for item in metadata["deliverables"]:
        path = outdir / item["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        if item["kind"] == "tsv":
            validate_tsv(path, item.get("required_columns", []))
        elif item["kind"] == "md":
            validate_markdown_sections(path, item.get("required_sections", []))
        else:
            raise AssertionError(f"Unsupported deliverable kind: {item['kind']}")

    for item in metadata.get("qc_artifacts", []):
        path = outdir / item["path"]
        if not path.exists():
            raise AssertionError(f"Missing QC artifact: {path}")
        if item["kind"] == "tsv":
            validate_tsv(path, item.get("required_columns", []))
        elif item["kind"] == "json":
            validate_json(path, item.get("required_keys", []))
        else:
            raise AssertionError(f"Unsupported QC artifact kind: {item['kind']}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate the cell-cell communication starter outputs.")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    validate_outputs(SKILL_DIR, args.outdir)
    print(json.dumps({"validated_outdir": str(args.outdir)}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
