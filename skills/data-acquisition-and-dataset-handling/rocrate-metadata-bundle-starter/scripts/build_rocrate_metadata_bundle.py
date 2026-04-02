#!/usr/bin/env python3
"""Build a tiny RO-Crate bundle from one local data file."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from rocrate.rocrate import ROCrate


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "toy_measurements.csv"


def build_crate(input_path: Path, crate_dir: Path, name: str, description: str, technique: str) -> dict[str, object]:
    if crate_dir.exists():
        shutil.rmtree(crate_dir)
    crate_dir.parent.mkdir(parents=True, exist_ok=True)

    crate = ROCrate()
    crate.name = name
    crate.description = description
    crate.root_dataset["measurementTechnique"] = technique
    crate.add_file(
        str(input_path),
        properties={
            "name": input_path.name,
            "description": "Toy tabular measurements bundled into an RO-Crate.",
            "encodingFormat": "text/csv",
        },
    )
    crate.write(str(crate_dir))

    metadata_path = crate_dir / "ro-crate-metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    graph = metadata.get("@graph", [])
    data_entities = [entry for entry in graph if entry.get("@type") == "File"]
    context_entities = [entry for entry in graph if entry.get("@type") not in {"Dataset", "File"}]
    root_dataset = next(entry for entry in graph if entry.get("@id") == "./")
    return {
        "input_path": str(input_path),
        "crate_dir": str(crate_dir),
        "has_metadata_file": metadata_path.exists(),
        "metadata_path": str(metadata_path),
        "root_dataset_name": root_dataset.get("name"),
        "measurement_technique": root_dataset.get("measurementTechnique"),
        "data_entity_count": len(data_entities),
        "context_entity_count": len(context_entities),
        "graph_entity_count": len(graph),
        "bundled_files": [entry.get("@id") for entry in data_entities],
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input data file to package.")
    parser.add_argument("--crate-dir", type=Path, required=True, help="Output RO-Crate directory.")
    parser.add_argument("--name", default="Toy Measurement Bundle", help="Root dataset name.")
    parser.add_argument(
        "--description",
        default="A deterministic RO-Crate bundle for local provenance testing.",
        help="Root dataset description.",
    )
    parser.add_argument(
        "--measurement-technique",
        default="synthetic parameter sweep",
        help="Lightweight measurement technique label.",
    )
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input file not found: {args.input}")

    payload = build_crate(args.input, args.crate_dir, args.name, args.description, args.measurement_technique)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
