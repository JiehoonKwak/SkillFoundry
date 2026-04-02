#!/usr/bin/env python3
"""Run a deterministic Scanpy marker-ranking pass on a toy grouped count matrix."""

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
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
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


def load_groups(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No group rows found in {path}")
    fieldnames = set(rows[0].keys())
    required = {"cell", "group"}
    if not required.issubset(fieldnames):
        missing = sorted(required - fieldnames)
        raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
    return {row["cell"]: row["group"] for row in rows}


def build_anndata(cells: list[str], genes: list[str], values: list[list[float]], groups: dict[str, str]):
    ad, np, pd, _ = load_scanpy_stack()
    missing = [cell for cell in cells if cell not in groups]
    if missing:
        raise ValueError(f"Missing group labels for cells: {', '.join(missing)}")
    matrix = np.asarray(values, dtype=np.float32).T
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes
    adata.obs["group"] = pd.Categorical([groups[cell] for cell in cells])
    return adata


def summarize_ranked_genes(adata, method: str, top_n: int) -> dict:
    _, _, _, sc = load_scanpy_stack()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(adata, groupby="group", method=method)
    groups = [str(group) for group in adata.obs["group"].cat.categories]
    top_markers_by_group = {}
    for group in groups:
        ranked = sc.get.rank_genes_groups_df(adata, group=group).head(top_n)
        top_markers_by_group[group] = [
            {
                "gene": row["names"],
                "score": round(float(row["scores"]), 6),
                "logfoldchange": round(float(row.get("logfoldchanges", 0.0)), 6),
                "pvals_adj": round(float(row.get("pvals_adj", 1.0)), 6),
            }
            for _, row in ranked.iterrows()
        ]
    return {
        "method": method,
        "top_n": top_n,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "groups": groups,
        "top_markers_by_group": top_markers_by_group,
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
    parser.add_argument("--input", type=Path, required=True, help="Tab-separated count matrix.")
    parser.add_argument("--groups", type=Path, required=True, help="Tab-separated cell-to-group labels.")
    parser.add_argument("--method", default="t-test", help="Ranking method, default: t-test.")
    parser.add_argument("--top-n", type=int, default=3, help="Number of genes to report per group.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.top_n <= 0:
        raise SystemExit("top-n must be a positive integer.")

    cells, genes, values = load_counts(args.input)
    groups = load_groups(args.groups)
    adata = build_anndata(cells, genes, values, groups)
    payload = summarize_ranked_genes(adata, method=args.method, top_n=args.top_n)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
