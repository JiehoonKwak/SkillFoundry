#!/usr/bin/env python3
"""Generate a deterministic tiny dataset bundle for tests and demos."""

from __future__ import annotations

import argparse
import csv
import json
import random
from pathlib import Path


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_bundle(sample_count: int, feature_count: int, seed: int, out_dir: Path) -> dict:
    if sample_count < 2:
        raise ValueError("sample_count must be at least 2")
    if feature_count < 1:
        raise ValueError("feature_count must be at least 1")

    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_rows = []
    for idx in range(sample_count):
        sample_rows.append(
            {
                "sample_id": f"S{idx + 1:02d}",
                "condition": "control" if idx % 2 == 0 else "treated",
                "batch": f"B{(idx % 2) + 1}",
            }
        )

    feature_rows = []
    for idx in range(feature_count):
        feature_rows.append(
            {
                "feature_id": f"F{idx + 1:02d}",
                "gene_symbol": f"GENE{idx + 1}",
                "pathway": "signal" if idx % 2 == 0 else "metabolism",
            }
        )

    matrix_rows: list[dict[str, object]] = []
    matrix_checksum = 0
    for feature in feature_rows:
        row: dict[str, object] = {"feature_id": feature["feature_id"]}
        for sample in sample_rows:
            value = rng.randint(10, 99)
            matrix_checksum += value
            row[str(sample["sample_id"])] = value
        matrix_rows.append(row)

    write_tsv(out_dir / "samples.tsv", ["sample_id", "condition", "batch"], sample_rows)
    write_tsv(out_dir / "features.tsv", ["feature_id", "gene_symbol", "pathway"], feature_rows)
    write_tsv(out_dir / "matrix.tsv", ["feature_id", *[row["sample_id"] for row in sample_rows]], matrix_rows)

    manifest = {
        "seed": seed,
        "files": {
            "samples": "samples.tsv",
            "features": "features.tsv",
            "matrix": "matrix.tsv",
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "bundle_dir": str(out_dir.resolve()),
        "seed": seed,
        "sample_count": sample_count,
        "feature_count": feature_count,
        "matrix_shape": [feature_count, sample_count],
        "matrix_checksum": matrix_checksum,
        "samples_preview": sample_rows[:2],
        "features_preview": feature_rows[:2],
        "manifest_files": manifest["files"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-count", type=int, required=True, help="Number of synthetic samples.")
    parser.add_argument("--feature-count", type=int, required=True, help="Number of synthetic features.")
    parser.add_argument("--seed", type=int, required=True, help="Random seed for deterministic values.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory for generated TSV/JSON files.")
    parser.add_argument("--summary-out", type=Path, required=True, help="Output JSON summary path.")
    args = parser.parse_args()

    summary = build_bundle(args.sample_count, args.feature_count, args.seed, args.out_dir)
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
