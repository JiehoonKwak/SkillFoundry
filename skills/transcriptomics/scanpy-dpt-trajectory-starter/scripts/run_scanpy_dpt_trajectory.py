#!/usr/bin/env python3
"""Compute a deterministic toy diffusion-pseudotime trajectory with Scanpy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_scanpy_stack():
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - exercised by runtime tests
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


def load_expected_order(path: Path | None) -> list[str] | None:
    if path is None:
        return None
    values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not values:
        raise ValueError(f"Expected-order file {path} must not be empty.")
    return values


def build_anndata(cells: list[str], genes: list[str], values: list[list[float]]):
    ad, np, _, _ = load_scanpy_stack()
    matrix = np.asarray(values, dtype=np.float32).T
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes
    return adata


def summarize_trajectory(counts_path: Path, root_cell: str, expected_order_path: Path | None) -> dict:
    _, _, _, sc = load_scanpy_stack()
    cells, genes, values = load_counts(counts_path)
    expected_order = load_expected_order(expected_order_path)
    adata = build_anndata(cells, genes, values)

    if root_cell not in adata.obs_names:
        raise ValueError(f"Root cell {root_cell!r} is not present in the count matrix.")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_pcs = max(2, min(3, adata.n_vars - 1, adata.n_obs - 1))
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=min(3, adata.n_obs - 1), n_pcs=n_pcs)
    n_comps = max(2, min(5, adata.n_obs - 1))
    sc.tl.diffmap(adata, n_comps=n_comps)
    adata.uns["iroot"] = int(adata.obs_names.get_loc(root_cell))
    sc.tl.dpt(adata, n_dcs=n_comps)

    pseudotime = {str(cell): round(float(value), 6) for cell, value in adata.obs["dpt_pseudotime"].items()}
    trajectory_order = [cell for cell, _ in sorted(pseudotime.items(), key=lambda item: item[1])]
    payload = {
        "cells": len(cells),
        "genes": len(genes),
        "root_cell": root_cell,
        "trajectory_order": trajectory_order,
        "pseudotime": [{"cell": cell, "pseudotime": pseudotime[cell]} for cell in trajectory_order],
        "pseudotime_span": {
            "min": pseudotime[trajectory_order[0]],
            "max": pseudotime[trajectory_order[-1]],
        },
    }
    if expected_order is not None:
        missing = [cell for cell in expected_order if cell not in pseudotime]
        if missing:
            raise ValueError(f"Expected-order cells missing from matrix: {', '.join(missing)}")
        monotonic = all(
            pseudotime[expected_order[index]] <= pseudotime[expected_order[index + 1]]
            for index in range(len(expected_order) - 1)
        )
        payload["expected_order"] = expected_order
        payload["monotonic_in_expected_order"] = monotonic
    return payload


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True, help="Tab-separated genes-by-cells count matrix.")
    parser.add_argument("--root-cell", required=True, help="Root cell for DPT.")
    parser.add_argument("--expected-order", type=Path, default=None, help="Optional file with expected cell order.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    payload = summarize_trajectory(args.counts, args.root_cell, args.expected_order)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
