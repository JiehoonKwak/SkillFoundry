#!/usr/bin/env python3
"""Run a minimal Scanpy QC pass on a toy count matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def load_counts(path: Path) -> tuple[list[str], list[str], np.ndarray]:
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
    matrix = np.asarray(values, dtype=np.float32).T
    return cells, genes, matrix


def build_anndata(cells: list[str], genes: list[str], matrix: np.ndarray) -> ad.AnnData:
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes
    adata.var["mt"] = pd.Index(adata.var_names).str.upper().str.startswith("MT-")
    return adata


def summarize_qc(adata: ad.AnnData) -> dict:
    obs = adata.obs.copy()
    per_cell = (
        obs.loc[:, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]]
        .reset_index(names="cell")
        .to_dict(orient="records")
    )
    return {
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "total_counts_sum": int(round(float(obs["total_counts"].sum()))),
        "median_total_counts": round(float(obs["total_counts"].median()), 4),
        "median_genes_by_counts": round(float(obs["n_genes_by_counts"].median()), 4),
        "max_pct_counts_mt": round(float(obs["pct_counts_mt"].max()), 4),
        "per_cell": [
            {
                "cell": record["cell"],
                "total_counts": round(float(record["total_counts"]), 4),
                "n_genes_by_counts": round(float(record["n_genes_by_counts"]), 4),
                "pct_counts_mt": round(float(record["pct_counts_mt"]), 4),
            }
            for record in per_cell
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Tab-separated count matrix")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional summary JSON output path")
    parser.add_argument("--h5ad-out", type=Path, default=None, help="Optional AnnData output path")
    args = parser.parse_args()

    cells, genes, matrix = load_counts(args.input)
    adata = build_anndata(cells, genes, matrix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    payload = summarize_qc(adata)
    text = json.dumps(payload, indent=2, sort_keys=True)

    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)

    if args.h5ad_out is not None:
        args.h5ad_out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.h5ad_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
