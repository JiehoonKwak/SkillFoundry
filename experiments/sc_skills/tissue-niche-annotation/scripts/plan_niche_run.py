#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def build_plan(
    *,
    adata: str,
    cell_type_column: str,
    x_key: str,
    y_key: str,
    sample_key: str | None,
    graph_mode: str,
    k: int | None,
    radius: float | None,
) -> dict[str, object]:
    if graph_mode not in {"knn", "radius"}:
        raise ValueError("graph_mode must be 'knn' or 'radius'.")
    if graph_mode == "knn" and (k is None or k < 1):
        raise ValueError("kNN mode requires --k >= 1.")
    if graph_mode == "radius" and (radius is None or radius <= 0):
        raise ValueError("radius mode requires --radius > 0.")

    return {
        "adata_path": adata,
        "cell_type_column": cell_type_column,
        "coordinate_keys": [x_key, y_key],
        "sample_key": sample_key,
        "graph": {
            "mode": graph_mode,
            "k": k if graph_mode == "knn" else None,
            "radius": radius if graph_mode == "radius" else None,
        },
        "expected_outputs": [
            "00_runtime_check.json",
            "01_dataset_profile.json",
            "02_niche_run_plan.json",
            "03_spatial_graph.tsv",
            "04_neighborhood_profiles.tsv",
            "05_niche_candidates.tsv",
            "06_niche_review.md",
            "niche_labels.tsv",
            "niche_markers.tsv",
            "tissue_niche_report.md",
        ],
        "decision_rules": [
            "Use kNN by default when cell density varies or coordinate units are uncertain.",
            "Use radius only when the coordinate system is calibrated and a biologically meaningful distance threshold is known.",
            "If the cell-type column is weak or missing, return to the annotation workflow first.",
            "If multiple slides or samples are present, keep a sample key and avoid cross-slide edges.",
        ],
        "stop_conditions": [
            "missing spatial coordinates",
            "missing cell-type/state annotation column",
            "spot-based platform instead of single-cell spatial",
            "multiple slides but no sample/slide key",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create a concrete run-plan scaffold for tissue niche annotation.")
    parser.add_argument("--adata", required=True, help="Input h5ad path.")
    parser.add_argument("--cell-type-column", required=True)
    parser.add_argument("--x-key", required=True)
    parser.add_argument("--y-key", required=True)
    parser.add_argument("--sample-key", default=None)
    parser.add_argument("--graph-mode", choices=["knn", "radius"], default="knn")
    parser.add_argument("--k", type=int, default=12)
    parser.add_argument("--radius", type=float, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    plan = build_plan(
        adata=args.adata,
        cell_type_column=args.cell_type_column,
        x_key=args.x_key,
        y_key=args.y_key,
        sample_key=args.sample_key,
        graph_mode=args.graph_mode,
        k=args.k,
        radius=args.radius,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(str(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
