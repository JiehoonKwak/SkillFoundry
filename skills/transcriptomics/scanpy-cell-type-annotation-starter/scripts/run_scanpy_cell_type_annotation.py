#!/usr/bin/env python3
"""Assign toy cell-type labels with Scanpy marker-program scoring."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
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


def load_markers(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not payload:
        raise ValueError(f"Marker mapping in {path} must be a non-empty JSON object.")
    markers: dict[str, list[str]] = {}
    for label, genes in payload.items():
        if not isinstance(label, str) or not label.strip():
            raise ValueError(f"Marker label in {path} must be a non-empty string.")
        if not isinstance(genes, list) or not genes:
            raise ValueError(f"Marker list for {label!r} in {path} must be a non-empty list.")
        normalized = [gene.strip() for gene in genes if isinstance(gene, str) and gene.strip()]
        if not normalized:
            raise ValueError(f"Marker list for {label!r} in {path} must contain non-empty gene names.")
        markers[label.strip()] = normalized
    return markers


def load_truth(path: Path | None) -> dict[str, str] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No truth rows found in {path}")
    if {"cell", "truth_label"} - set(rows[0].keys()):
        raise ValueError(f"Truth file {path} must contain columns: cell, truth_label")
    return {row["cell"]: row["truth_label"] for row in rows}


def build_anndata(cells: list[str], genes: list[str], values: list[list[float]]):
    ad, np, _, _ = load_scanpy_stack()
    matrix = np.asarray(values, dtype=np.float32).T
    adata = ad.AnnData(X=matrix)
    adata.obs_names = cells
    adata.var_names = genes
    return adata


def score_labels(adata, markers: dict[str, list[str]]) -> tuple[dict[str, str], dict[str, dict[str, float]]]:
    _, np, _, sc = load_scanpy_stack()
    score_columns: dict[str, str] = {}
    for label, genes in markers.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if not present:
            raise ValueError(f"No marker genes from label {label!r} were found in the count matrix.")
        score_name = f"{label.lower().replace(' ', '_')}_score"
        sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=False)
        score_columns[label] = score_name

    predictions: dict[str, str] = {}
    score_table: dict[str, dict[str, float]] = {}
    for cell in adata.obs_names:
        cell_scores = {
            label: round(float(adata.obs.loc[cell, score_name]), 6)
            for label, score_name in score_columns.items()
        }
        score_table[str(cell)] = cell_scores
        predictions[str(cell)] = max(cell_scores.items(), key=lambda item: item[1])[0]
    return predictions, score_table


def summarize_annotation(counts_path: Path, markers_path: Path, truth_path: Path | None) -> dict:
    _, _, _, sc = load_scanpy_stack()
    cells, genes, values = load_counts(counts_path)
    markers = load_markers(markers_path)
    truth = load_truth(truth_path)
    adata = build_anndata(cells, genes, values)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    predictions, score_table = score_labels(adata, markers)
    predicted_counts = dict(sorted(Counter(predictions.values()).items()))
    payload = {
        "cells": len(cells),
        "genes": len(genes),
        "labels": list(markers.keys()),
        "predicted_label_counts": predicted_counts,
        "predictions": [
            {
                "cell": cell,
                "predicted_label": predictions[cell],
                "truth_label": truth.get(cell) if truth is not None else None,
                "scores": score_table[cell],
            }
            for cell in cells
        ],
        "mean_label_scores": {
            label: round(sum(score_table[cell][label] for cell in cells) / len(cells), 6)
            for label in markers
        },
    }
    if truth is not None:
        missing = [cell for cell in cells if cell not in truth]
        if missing:
            raise ValueError(f"Truth labels missing cells: {', '.join(missing)}")
        accuracy = sum(predictions[cell] == truth[cell] for cell in cells) / len(cells)
        payload["accuracy"] = round(float(accuracy), 6)
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
    parser.add_argument("--markers", type=Path, required=True, help="JSON mapping label to marker genes.")
    parser.add_argument("--truth", type=Path, default=None, help="Optional TSV with columns: cell, truth_label.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    payload = summarize_annotation(args.counts, args.markers, args.truth)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
