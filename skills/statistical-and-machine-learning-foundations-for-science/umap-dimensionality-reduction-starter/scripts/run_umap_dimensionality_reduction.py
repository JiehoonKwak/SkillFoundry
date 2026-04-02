#!/usr/bin/env python3
"""Run a deterministic UMAP embedding over a tiny labeled matrix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import umap


def load_embedding_input(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    sample_ids: list[str] = []
    labels: list[str] = []
    features: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_ids.append(row["sample_id"])
            labels.append(row["label"])
            features.append([float(row[key]) for key in ("f1", "f2", "f3", "f4")])
    return sample_ids, labels, np.asarray(features, dtype=float)


def run_embedding(input_path: Path) -> dict:
    sample_ids, labels, matrix = load_embedding_input(input_path)
    reducer = umap.UMAP(n_neighbors=3, min_dist=0.01, random_state=42, metric="euclidean")
    embedding = reducer.fit_transform(matrix)
    unique_labels = sorted(set(labels))
    centroids = {}
    for label in unique_labels:
        coords = embedding[[index for index, value in enumerate(labels) if value == label]]
        centroids[label] = [round(float(value), 6) for value in coords.mean(axis=0)]
    centroid_distance = float(np.linalg.norm(np.asarray(centroids[unique_labels[0]]) - np.asarray(centroids[unique_labels[1]])))
    points = []
    for sample_id, label, coords in zip(sample_ids, labels, embedding):
        points.append({"sample_id": sample_id, "label": label, "x": round(float(coords[0]), 6), "y": round(float(coords[1]), 6)})
    return {
        "input_path": str(input_path.resolve()),
        "sample_count": len(sample_ids),
        "labels": unique_labels,
        "points": points,
        "centroids": centroids,
        "centroid_distance": round(centroid_distance, 6),
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_embedding(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
