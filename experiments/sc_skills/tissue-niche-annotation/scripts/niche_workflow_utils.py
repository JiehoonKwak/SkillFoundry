#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


COORDINATE_OBS_PAIRS = [
    ("x", "y"),
    ("x_coord", "y_coord"),
    ("x_centroid", "y_centroid"),
    ("center_x", "center_y"),
    ("spatial_x", "spatial_y"),
    ("X_centroid", "Y_centroid"),
]
LABEL_HINTS = ["cell_type", "celltype", "annotation", "state", "subtype", "class", "cluster"]
SAMPLE_HINTS = ["sample", "slide", "library", "fov", "roi", "section", "image"]
BATCH_HINTS = ["batch", "patient", "donor", "replicate", "condition"]


def load_anndata(path: Path) -> ad.AnnData:
    return ad.read_h5ad(path)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sanitize_label(label: str) -> str:
    cleaned = "".join(character if character.isalnum() else "_" for character in str(label).lower())
    return cleaned.strip("_") or "label"


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    totals = np.where(totals <= 0, 1.0, totals)
    return matrix / totals


def softmax_rows(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def probability_margin(probabilities: np.ndarray) -> np.ndarray:
    if probabilities.shape[1] == 1:
        return np.ones(probabilities.shape[0], dtype=float)
    top_two = np.partition(probabilities, -2, axis=1)[:, -2:]
    return top_two[:, 1] - top_two[:, 0]


def normalized_entropy(probabilities: np.ndarray) -> np.ndarray:
    if probabilities.shape[1] == 1:
        return np.zeros(probabilities.shape[0], dtype=float)
    safe = np.clip(probabilities, 1e-12, 1.0)
    entropy = -np.sum(safe * np.log(safe), axis=1)
    return entropy / math.log(probabilities.shape[1])


def _string_series(series: pd.Series) -> pd.Series:
    values = series.astype("string").fillna(pd.NA)
    return values


def _column_name_score(name: str, hints: list[str]) -> int:
    lowered = name.lower()
    return sum(3 for hint in hints if hint in lowered)


def infer_candidate_columns(
    obs: pd.DataFrame,
    *,
    hints: list[str],
    max_unique: int,
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    total = max(len(obs), 1)
    for column in obs.columns:
        series = _string_series(obs[column])
        non_missing = series.dropna()
        unique_count = int(non_missing.nunique())
        if unique_count == 0 or unique_count > max_unique:
            continue
        name_score = _column_name_score(str(column), hints)
        if name_score <= 0 and unique_count > min(50, max(3, total // 2)):
            continue
        missing_fraction = float(series.isna().mean())
        dominant_fraction = float(non_missing.value_counts(normalize=True).iloc[0]) if not non_missing.empty else 1.0
        likely_annotation = (
            name_score > 0
            and 2 <= unique_count <= min(100, max(2, total // 2))
            and missing_fraction <= 0.25
        )
        candidates.append(
            {
                "column": str(column),
                "name_score": name_score,
                "unique_count": unique_count,
                "missing_fraction": round(missing_fraction, 6),
                "dominant_fraction": round(dominant_fraction, 6),
                "likely_annotation": likely_annotation,
            }
        )
    return sorted(
        candidates,
        key=lambda item: (
            not item["likely_annotation"],
            -int(item["name_score"]),
            float(item["missing_fraction"]),
            int(item["unique_count"]),
            item["column"],
        ),
    )


def choose_default_candidate(candidates: list[dict[str, Any]]) -> str | None:
    if not candidates:
        return None
    return str(candidates[0]["column"])


def detect_coordinates(
    adata: ad.AnnData,
    *,
    x_key: str | None = None,
    y_key: str | None = None,
    spatial_key: str = "spatial",
) -> dict[str, Any]:
    obs = adata.obs
    if x_key is not None and y_key is not None and x_key in obs.columns and y_key in obs.columns:
        coords = obs[[x_key, y_key]].to_numpy(dtype=float)
        return {
            "found": True,
            "source": "obs",
            "x_key": x_key,
            "y_key": y_key,
            "spatial_key": None,
            "coords": coords,
        }

    if spatial_key in adata.obsm and np.asarray(adata.obsm[spatial_key]).shape[1] >= 2:
        coords = np.asarray(adata.obsm[spatial_key])[:, :2].astype(float)
        return {
            "found": True,
            "source": "obsm",
            "x_key": f"{spatial_key}[:,0]",
            "y_key": f"{spatial_key}[:,1]",
            "spatial_key": spatial_key,
            "coords": coords,
        }

    for obs_x, obs_y in COORDINATE_OBS_PAIRS:
        if obs_x in obs.columns and obs_y in obs.columns:
            coords = obs[[obs_x, obs_y]].to_numpy(dtype=float)
            return {
                "found": True,
                "source": "obs",
                "x_key": obs_x,
                "y_key": obs_y,
                "spatial_key": None,
                "coords": coords,
            }

    return {
        "found": False,
        "source": None,
        "x_key": None,
        "y_key": None,
        "spatial_key": spatial_key if spatial_key in adata.obsm else None,
        "coords": None,
    }


def resolve_sample_ids(obs: pd.DataFrame, sample_key: str | None) -> tuple[np.ndarray, str]:
    if sample_key is not None and sample_key in obs.columns:
        values = _string_series(obs[sample_key]).fillna("missing_sample").to_numpy(dtype=str)
        return values, sample_key
    return np.asarray(["sample_001"] * len(obs), dtype=str), "__synthetic_single_sample__"


def resolve_label_values(obs: pd.DataFrame, label_key: str | None) -> tuple[np.ndarray | None, list[str]]:
    if label_key is None or label_key not in obs.columns:
        return None, []
    series = _string_series(obs[label_key])
    values = series.fillna("missing_label").to_numpy(dtype=str)
    categories = sorted(pd.Index(values).unique().tolist())
    return values, categories


def nearest_neighbor_summary(coords: np.ndarray, *, neighbor_k: int) -> dict[str, Any]:
    if len(coords) < 2:
        return {
            "median_nn_distance": None,
            "cv_nn_distance": None,
            "median_k_distance": None,
            "p90_k_distance": None,
        }
    effective_k = max(1, min(int(neighbor_k), len(coords) - 1))
    model = NearestNeighbors(n_neighbors=effective_k + 1, metric="euclidean")
    model.fit(coords)
    distances, _ = model.kneighbors(coords)
    nearest = distances[:, 1]
    kth = distances[:, effective_k]
    mean_nearest = float(nearest.mean())
    return {
        "median_nn_distance": round(float(np.median(nearest)), 6),
        "cv_nn_distance": round(float(nearest.std() / mean_nearest), 6) if mean_nearest > 0 else 0.0,
        "median_k_distance": round(float(np.median(kth)), 6),
        "p90_k_distance": round(float(np.percentile(kth, 90)), 6),
    }


def build_sample_summaries(
    *,
    sample_ids: np.ndarray,
    coords: np.ndarray,
    label_values: np.ndarray | None,
    neighbor_k: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sample in sorted(pd.Index(sample_ids).unique().tolist()):
        mask = sample_ids == sample
        sample_coords = coords[mask]
        sample_labels = label_values[mask] if label_values is not None else None
        label_counts = (
            pd.Series(sample_labels).value_counts().to_dict()
            if sample_labels is not None
            else {}
        )
        distance_summary = nearest_neighbor_summary(sample_coords, neighbor_k=neighbor_k)
        row = {
            "sample_id": sample,
            "cell_count": int(mask.sum()),
            "x_min": round(float(sample_coords[:, 0].min()), 6),
            "x_max": round(float(sample_coords[:, 0].max()), 6),
            "y_min": round(float(sample_coords[:, 1].min()), 6),
            "y_max": round(float(sample_coords[:, 1].max()), 6),
            "unique_label_count": int(len(label_counts)),
            "top_labels": [
                {"label": str(label), "count": int(count)}
                for label, count in list(label_counts.items())[:5]
            ],
        }
        row.update(distance_summary)
        rows.append(row)
    return rows


def recommend_graph_parameters(
    sample_summaries: list[dict[str, Any]],
    *,
    n_obs: int,
) -> dict[str, Any]:
    suggested_k = 6 if n_obs <= 5000 else 8
    median_k_distances = [
        float(item["median_k_distance"])
        for item in sample_summaries
        if item["median_k_distance"] is not None and float(item["median_k_distance"]) > 0
    ]
    nn_cvs = [
        float(item["cv_nn_distance"])
        for item in sample_summaries
        if item["cv_nn_distance"] is not None
    ]
    density_ratio = (
        max(median_k_distances) / min(median_k_distances)
        if len(median_k_distances) >= 2 and min(median_k_distances) > 0
        else 1.0
    )
    overall_cv = float(np.mean(nn_cvs)) if nn_cvs else 0.0
    if len(sample_summaries) > 1 and density_ratio >= 1.5:
        mode = "knn"
        reason = "per-sample spatial density differs enough that a shared radius would connect samples unevenly"
    elif overall_cv >= 0.45:
        mode = "knn"
        reason = "local point spacing is heterogeneous, so a fixed neighbor count is safer than a single radius"
    else:
        mode = "radius"
        reason = "local spacing is reasonably stable, so a shared radius can preserve geometric meaning"

    radius = None
    if median_k_distances:
        radius = round(float(np.median(median_k_distances) * 1.05), 6)
    return {
        "mode": mode,
        "reason": reason,
        "neighbor_k": suggested_k,
        "radius": radius,
        "sample_density_ratio": round(float(density_ratio), 6),
        "mean_nn_cv": round(overall_cv, 6),
    }


def recommend_composition_mode(
    *,
    label_values: np.ndarray | None,
    graph_mode: str,
) -> dict[str, str]:
    if label_values is None:
        return {
            "primary_view": "smoothed",
            "reason": "no label column was selected, so neighborhood features cannot be trusted yet",
        }
    counts = pd.Series(label_values).value_counts(normalize=True)
    top_fraction = float(counts.iloc[0]) if not counts.empty else 1.0
    unique_count = int(counts.shape[0])
    if graph_mode == "knn" or unique_count >= 4 or top_fraction <= 0.7:
        return {
            "primary_view": "smoothed",
            "reason": "mixed or multi-compartment neighborhoods benefit from graph smoothing after saving the raw composition first",
        }
    return {
        "primary_view": "raw",
        "reason": "label composition is simple enough that the unsmoothed neighborhood fractions should remain interpretable",
    }


def build_weighted_spatial_graph(
    *,
    cell_ids: list[str],
    sample_ids: list[str],
    coords: np.ndarray,
    mode: str,
    neighbor_k: int,
    radius: float | None,
    self_weight: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], dict[str, Any]]:
    if not 0.0 <= self_weight < 1.0:
        raise AssertionError("self_weight must be in [0, 1).")
    if mode not in {"knn", "radius"}:
        raise AssertionError("mode must be 'knn' or 'radius'.")

    n_cells = len(cell_ids)
    neighbor_graph = np.zeros((n_cells, n_cells), dtype=float)
    combined_graph = np.zeros((n_cells, n_cells), dtype=float)
    rows: list[dict[str, Any]] = []
    isolated_cells = 0
    fallback_cells = 0
    sample_array = np.asarray(sample_ids, dtype=str)

    for sample in sorted(pd.Index(sample_array).unique().tolist()):
        sample_indices = np.where(sample_array == sample)[0]
        sample_coords = coords[sample_indices]
        if len(sample_indices) < 2:
            raise AssertionError(f"Sample {sample} needs at least two cells for a spatial graph.")

        for local_position, global_index in enumerate(sample_indices):
            source_coord = sample_coords[local_position]
            distances = np.linalg.norm(sample_coords - source_coord, axis=1)
            non_self = [index for index in np.argsort(distances) if index != local_position]
            effective_k = max(1, min(int(neighbor_k), len(non_self)))

            if mode == "knn":
                chosen_local = non_self[:effective_k]
            else:
                if radius is None or radius <= 0:
                    raise AssertionError("radius mode requires a positive radius.")
                chosen_local = [
                    index for index in non_self
                    if float(distances[index]) <= float(radius)
                ]
                if not chosen_local:
                    chosen_local = non_self[:effective_k]
                    fallback_cells += 1
                if len(chosen_local) < effective_k:
                    isolated_cells += 1

            chosen_distances = np.asarray([float(distances[index]) for index in chosen_local], dtype=float)
            if chosen_distances.size == 0:
                raise AssertionError(f"Sample {sample} produced an empty neighborhood for cell {cell_ids[global_index]}.")
            if np.allclose(chosen_distances, 0.0):
                raw_weights = np.full(chosen_distances.shape, 1.0 / chosen_distances.size, dtype=float)
            else:
                raw_weights = 1.0 / np.clip(chosen_distances, 1e-9, None)
                raw_weights = raw_weights / raw_weights.sum()

            combined_graph[global_index, global_index] = float(self_weight)
            rows.append(
                {
                    "source_cell": cell_ids[global_index],
                    "target_cell": cell_ids[global_index],
                    "sample_id": sample,
                    "relation": "self",
                    "distance": 0.0,
                    "weight": round(float(self_weight), 6),
                    "graph_mode": mode,
                }
            )
            for target_local, normalized_weight in zip(chosen_local, raw_weights, strict=True):
                target_index = int(sample_indices[target_local])
                neighbor_graph[global_index, target_index] = float(normalized_weight)
                combined_graph[global_index, target_index] = float((1.0 - self_weight) * normalized_weight)
                rows.append(
                    {
                        "source_cell": cell_ids[global_index],
                        "target_cell": cell_ids[target_index],
                        "sample_id": sample,
                        "relation": "neighbor",
                        "distance": round(float(distances[target_local]), 6),
                        "weight": round(float((1.0 - self_weight) * normalized_weight), 6),
                        "graph_mode": mode,
                    }
                )

    if not np.allclose(neighbor_graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Neighbor graph rows must sum to 1.")
    if not np.allclose(combined_graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Combined graph rows must sum to 1.")

    graph_qc = {
        "graph_mode": mode,
        "neighbor_k": int(neighbor_k),
        "radius": None if radius is None else round(float(radius), 6),
        "self_weight": round(float(self_weight), 6),
        "isolated_or_sparse_radius_cells": int(isolated_cells),
        "radius_backstop_cells": int(fallback_cells),
    }
    return neighbor_graph, combined_graph, rows, graph_qc


def compute_neighborhood_profiles(
    *,
    label_values: list[str],
    label_categories: list[str],
    neighbor_graph: np.ndarray,
    combined_graph: np.ndarray,
    smoothing_weight: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    label_index = {label: index for index, label in enumerate(label_categories)}
    one_hot = np.zeros((len(label_values), len(label_categories)), dtype=float)
    for row_index, label in enumerate(label_values):
        one_hot[row_index, label_index[label]] = 1.0

    local_profile = neighbor_graph @ one_hot
    smoothed_profile = (
        (1.0 - float(smoothing_weight)) * local_profile
        + float(smoothing_weight) * (combined_graph @ local_profile)
    )
    smoothed_profile = row_normalize(smoothed_profile)
    dominant_neighbor_indices = local_profile.argmax(axis=1)
    dominant_neighbor_type = [label_categories[index] for index in dominant_neighbor_indices]
    top_neighbor_fraction = local_profile.max(axis=1)
    boundary_score = normalized_entropy(smoothed_profile)
    return local_profile, smoothed_profile, top_neighbor_fraction, boundary_score, dominant_neighbor_type

