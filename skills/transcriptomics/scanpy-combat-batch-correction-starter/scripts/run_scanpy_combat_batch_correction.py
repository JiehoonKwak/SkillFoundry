#!/usr/bin/env python3
"""Run a deterministic Scanpy ComBat batch-correction pass on a toy matrix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_scanpy_stack():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - exercised through runtime tests
        raise SystemExit(
            "This script requires the dedicated Scanpy prefix at slurm/envs/scanpy. "
            "Run it with slurm/envs/scanpy/bin/python."
        ) from exc
    return ad, np, pd, sc


def load_counts(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    rows = path.read_text(encoding="utf-8").strip().splitlines()
    header = rows[0].split("\t")
    cells = header[1:]
    genes: list[str] = []
    values: list[list[float]] = []
    for raw_row in rows[1:]:
        fields = raw_row.split("\t")
        if len(fields) != len(header):
            raise ValueError(f"Ragged row in {path}: {raw_row}")
        genes.append(fields[0])
        values.append([float(value) for value in fields[1:]])
    return cells, genes, values


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No metadata rows found in {path}")
    required = {"cell", "batch", "cell_type"}
    fieldnames = set(rows[0].keys())
    missing = sorted(required - fieldnames)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    return {row["cell"]: row for row in rows}


def build_anndata(cells: list[str], genes: list[str], values: list[list[float]], metadata: dict[str, dict[str, str]]):
    ad, np, pd, _ = load_scanpy_stack()
    missing = [cell for cell in cells if cell not in metadata]
    if missing:
        raise ValueError(f"Missing metadata rows for cells: {', '.join(missing)}")
    matrix = np.asarray(values, dtype=np.float32).T
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs["batch"] = pd.Categorical([metadata[cell]["batch"] for cell in cells])
    adata.obs["cell_type"] = pd.Categorical([metadata[cell]["cell_type"] for cell in cells])
    return adata


def _group_means(adata, column: str) -> dict[str, list[float]]:
    _, np, _, _ = load_scanpy_stack()
    groups = [str(value) for value in adata.obs[column].cat.categories]
    return {
        group: [round(float(value), 6) for value in np.asarray(adata[adata.obs[column] == group].X.mean(axis=0)).ravel()]
        for group in groups
    }


def _mean_abs_diff_by_batch(adata) -> float:
    _, np, _, _ = load_scanpy_stack()
    batch_groups = [str(value) for value in adata.obs["batch"].cat.categories]
    if len(batch_groups) != 2:
        raise ValueError("This starter expects exactly two batches in the toy metadata.")
    batch_means = {
        group: np.asarray(adata[adata.obs["batch"] == group].X.mean(axis=0)).ravel()
        for group in batch_groups
    }
    return round(float(np.abs(batch_means[batch_groups[0]] - batch_means[batch_groups[1]]).mean()), 6)


def _cell_type_centroid_distance(adata) -> float | None:
    _, np, _, _ = load_scanpy_stack()
    cell_types = [str(value) for value in adata.obs["cell_type"].cat.categories]
    if len(cell_types) < 2:
        return None
    centroids = {
        cell_type: np.asarray(adata[adata.obs["cell_type"] == cell_type].X.mean(axis=0)).ravel()
        for cell_type in cell_types
    }
    return round(float(np.linalg.norm(centroids[cell_types[0]] - centroids[cell_types[1]])), 6)


def summarize_batch_correction(adata) -> dict:
    _, _, _, sc = load_scanpy_stack()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    pre_gap = _mean_abs_diff_by_batch(adata)
    pre_means = _group_means(adata, "batch")
    sc.pp.combat(adata, key="batch")
    post_gap = _mean_abs_diff_by_batch(adata)
    post_means = _group_means(adata, "batch")
    centroid_distance = _cell_type_centroid_distance(adata)
    return {
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "batches": [str(value) for value in adata.obs["batch"].cat.categories],
        "cell_types": [str(value) for value in adata.obs["cell_type"].cat.categories],
        "pre_batch_mean_abs_diff": pre_gap,
        "post_batch_mean_abs_diff": post_gap,
        "batch_gap_ratio": round(post_gap / pre_gap, 6) if pre_gap else 0.0,
        "cell_type_centroid_distance_post": centroid_distance,
        "pre_batch_gene_means": pre_means,
        "post_batch_gene_means": post_means,
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
    parser.add_argument("--counts", type=Path, required=True, help="Tab-separated count matrix.")
    parser.add_argument("--metadata", type=Path, required=True, help="Tab-separated cell metadata.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    cells, genes, values = load_counts(args.counts)
    metadata = load_metadata(args.metadata)
    adata = build_anndata(cells, genes, values, metadata)
    payload = summarize_batch_correction(adata)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
