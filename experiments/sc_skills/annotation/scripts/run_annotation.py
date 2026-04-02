#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from reference_search import (
    build_reference_query_string,
    build_reference_query_variants,
    build_reference_selection_summary,
    rank_reference_candidates,
)


SKILL_DIR = Path(__file__).resolve().parents[1]
ANNOTATION_COLUMNS = [
    "cell_id",
    "sample_id",
    "x_coord",
    "y_coord",
    "predicted_label",
    "label_source",
    "label_confidence",
    "broad_label",
    "marker_status",
    "notes",
]
QC_FILES = [
    "qc_dataset_exploration.json",
    "qc_reference_search.tsv",
    "qc_reference_selection.json",
    "qc_label_transfer.tsv",
    "run_summary.json",
]
LABEL_TRANSFER_COLUMNS = [
    "cell_id",
    "sample_id",
    "nearest_reference_ids",
    "nearest_reference_broad_labels",
    "nearest_reference_subtypes",
    "nearest_reference_batches",
    "broad_supports",
    "broad_margin",
    "broad_label",
    "subtype_supports",
    "subtype_margin",
    "predicted_label",
    "label_confidence",
    "broad_marker_gap",
    "subtype_marker_gap",
    "marker_status",
    "latent_1",
    "latent_2",
    "latent_3",
    "notes",
]
UTAG_COLUMNS = [
    "cell_id",
    "sample_id",
    "neighbor_ids",
    "neighbor_broad_labels",
    "neighbor_weight_sum",
    "local_broad_profile",
    "smoothed_broad_profile",
    "boundary_entropy",
    "label_only_niche",
    "local_best_niche",
    "utag_niche",
    "graph_shifted",
    "niche_confidence",
    "prototype_scores",
    "niche_reason",
]
METHOD_STEPS = [
    "dataset_exploration_and_platform_detection",
    "normalize_total_log1p_and_pca_surrogate",
    "cellxgene_reference_search_and_materialization",
    "harmony_surrogate_label_transfer",
    "hierarchical_marker_review",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_metadata(skill_dir: Path) -> dict[str, Any]:
    return load_json(skill_dir / "metadata.yaml")


def validate_markdown_sections(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def write_markdown(path: Path, title: str, sections: list[dict[str, Any]]) -> None:
    lines = [f"# {title}", ""]
    for section in sections:
        lines.append(f"## {section['name']}")
        lines.append("")
        for paragraph in section.get("paragraphs", []):
            lines.append(paragraph)
            lines.append("")
        for bullet in section.get("bullets", []):
            lines.append(f"- {bullet}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def format_counts(counts: dict[str, int]) -> str:
    return ", ".join(f"{key}={counts[key]}" for key in sorted(counts))


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise AssertionError("Every row must have positive mass.")
    return matrix / totals


def softmax_rows(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def normalized_entropy(probabilities: np.ndarray) -> np.ndarray:
    safe = np.clip(probabilities, 1e-9, 1.0)
    entropy = -np.sum(safe * np.log(safe), axis=1)
    return entropy / math.log(probabilities.shape[1])


def pairwise_euclidean(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    deltas = left[:, None, :] - right[None, :, :]
    return np.sqrt(np.square(deltas).sum(axis=2))


def normalize_log1p(matrix: np.ndarray, target_sum: float) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = matrix.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise AssertionError("Every row must have positive total counts.")
    normalized = (matrix / library_sizes[:, None]) * float(target_sum)
    return np.log1p(normalized), library_sizes


def build_shared_gene_order(query_genes: list[str], reference_genes: list[str]) -> list[str]:
    reference_set = set(reference_genes)
    return [gene for gene in query_genes if gene in reference_set]


def matrix_from_rows(rows: list[dict[str, Any]], value_key: str) -> np.ndarray:
    return np.asarray([row[value_key] for row in rows], dtype=float)


def build_pca_surrogate(reference_matrix: np.ndarray, query_matrix: np.ndarray, latent_dim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    combined = np.vstack([reference_matrix, query_matrix])
    centered = combined - combined.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    max_dim = min(latent_dim, vt.shape[0], max(combined.shape[0] - 1, 1))
    if max_dim < 1:
        raise AssertionError("At least one latent dimension is required.")
    components = vt[:max_dim].copy()
    latent = centered @ components.T
    for dim in range(max_dim):
        if latent[0, dim] < 0:
            latent[:, dim] *= -1.0
            components[dim, :] *= -1.0
    variance = singular_values**2
    total_variance = float(variance.sum())
    explained = np.zeros(max_dim, dtype=float) if total_variance <= 0 else variance[:max_dim] / total_variance
    reference_latent = latent[: reference_matrix.shape[0]].copy()
    query_latent = latent[reference_matrix.shape[0] :].copy()
    return reference_latent, query_latent, explained[:max_dim], components


def harmony_surrogate(latent: np.ndarray, batches: list[str], shrinkage: float) -> tuple[np.ndarray, dict[str, float]]:
    corrected = latent.copy()
    overall_mean = latent.mean(axis=0, keepdims=True)
    batch_shift_norms: dict[str, float] = {}
    for batch in sorted(set(batches)):
        batch_indices = [index for index, value in enumerate(batches) if value == batch]
        batch_mean = latent[batch_indices].mean(axis=0, keepdims=True)
        shift = float(shrinkage) * (batch_mean - overall_mean)
        corrected[batch_indices] = latent[batch_indices] - shift
        batch_shift_norms[batch] = round(float(np.linalg.norm(shift)), 6)
    return corrected, batch_shift_norms


def ranked_string(items: list[tuple[str, float]]) -> str:
    return ";".join(f"{label}={value:.6f}" for label, value in items)


def weighted_vote(
    *,
    distance_row: np.ndarray,
    reference_ids: list[str],
    labels: list[str],
    reference_batches: list[str],
    k: int,
    allowed_indices: list[int] | None = None,
) -> dict[str, Any]:
    if allowed_indices is None:
        allowed_indices = list(range(len(reference_ids)))
    if not allowed_indices:
        raise AssertionError("Weighted vote requires at least one candidate reference.")
    ordered_indices = sorted(
        allowed_indices,
        key=lambda index: (float(distance_row[index]), reference_ids[index]),
    )[: min(int(k), len(allowed_indices))]
    label_weights: dict[str, float] = {}
    for index in ordered_indices:
        label = labels[index]
        weight = 1.0 / max(float(distance_row[index]), 1e-9)
        label_weights[label] = label_weights.get(label, 0.0) + weight
    total_weight = float(sum(label_weights.values()))
    ranked_votes = sorted(
        ((label, weight / total_weight) for label, weight in label_weights.items()),
        key=lambda item: (-item[1], item[0]),
    )
    predicted_label = ranked_votes[0][0]
    margin = float(ranked_votes[0][1] - (ranked_votes[1][1] if len(ranked_votes) > 1 else 0.0))
    return {
        "neighbor_ids": [reference_ids[index] for index in ordered_indices],
        "neighbor_labels": [labels[index] for index in ordered_indices],
        "neighbor_batches": [reference_batches[index] for index in ordered_indices],
        "supports": ranked_votes,
        "predicted_label": predicted_label,
        "margin": margin,
    }


def score_marker_gap(
    *,
    row_vector: np.ndarray,
    predicted_label: str,
    marker_map: dict[str, list[str]],
    shared_gene_to_index: dict[str, int],
    competitor_labels: list[str],
) -> dict[str, Any]:
    predicted_markers = [gene for gene in marker_map[predicted_label] if gene in shared_gene_to_index]
    if not predicted_markers:
        raise AssertionError(f"No shared markers available for {predicted_label}.")
    predicted_score = float(np.mean([row_vector[shared_gene_to_index[gene]] for gene in predicted_markers]))
    if competitor_labels:
        competitor_items = []
        for label in competitor_labels:
            genes = [gene for gene in marker_map[label] if gene in shared_gene_to_index]
            if not genes:
                continue
            score = float(np.mean([row_vector[shared_gene_to_index[gene]] for gene in genes]))
            competitor_items.append((label, score))
        if competitor_items:
            competitor_label, competitor_score = max(competitor_items, key=lambda item: (item[1], item[0]))
        else:
            competitor_label, competitor_score = "none", 0.0
    else:
        competitor_label, competitor_score = "none", 0.0
    return {
        "predicted_markers": predicted_markers,
        "predicted_score": predicted_score,
        "competitor_label": competitor_label,
        "competitor_score": competitor_score,
        "gap": float(predicted_score - competitor_score),
    }


def detect_platform_and_explore(payload: dict[str, Any]) -> dict[str, Any]:
    metadata = payload["query_metadata"]
    query = payload["query"]
    cells = query["cells"]
    genes = query["genes"]
    counts = matrix_from_rows(cells, "counts")
    obs_columns = sorted(key for key in cells[0] if key != "counts")
    var_columns = sorted(query["var_rows"][0].keys())
    coordinate_keys = [key for key in metadata["coordinate_keys"] if key in obs_columns]
    candidate_rows = []
    for column in obs_columns:
        if column in {"cell_id", *coordinate_keys, "compartment_hint"}:
            continue
        values = [str(cell[column]) for cell in cells]
        unique_values = sorted(set(values))
        unique_count = len(unique_values)
        if unique_count <= 1 or unique_count >= len(cells):
            continue
        score = 0.0
        lowered = column.lower()
        if "sample" in lowered:
            score += 4.0
        if "batch" in lowered:
            score += 3.0
        if "patient" in lowered or "donor" in lowered:
            score += 2.0
        if unique_count <= max(3, len(cells) // 2):
            score += 1.0
        if unique_count == 2:
            score += 0.5
        candidate_rows.append(
            {
                "column": column,
                "unique_count": unique_count,
                "values": unique_values,
                "score": round(score, 6),
            }
        )
    candidate_rows = sorted(candidate_rows, key=lambda row: (-float(row["score"]), row["unique_count"], row["column"]))
    selected_batch_column = metadata.get("suggested_sample_column")
    if selected_batch_column not in {row["column"] for row in candidate_rows} and candidate_rows:
        selected_batch_column = candidate_rows[0]["column"]
    elif not candidate_rows:
        selected_batch_column = metadata.get("suggested_sample_column", "sample_id")

    supported_platforms = {item.lower() for item in metadata["supported_single_cell_platforms"]}
    spot_platforms = {item.lower() for item in metadata["spot_platforms"]}
    platform_name = str(metadata["platform_name"])
    platform_name_lower = platform_name.lower()
    reasons = []
    is_supported = False
    if platform_name_lower in supported_platforms:
        reasons.append(f"platform_name={platform_name} is a supported single-cell platform")
        is_supported = True
    if platform_name_lower in spot_platforms:
        reasons.append(f"platform_name={platform_name} is a known spot-based platform")
        is_supported = False
    if str(metadata["observation_unit"]).lower() == "cell":
        reasons.append("observation_unit is cell")
        is_supported = is_supported or True
    if len(coordinate_keys) == 2:
        reasons.append("x/y coordinates are present")
    if len(genes) <= 1000:
        reasons.append("targeted gene panel is compatible with single-cell spatial assays")
    if counts.max() > float(payload["preprocessing"]["normalize_if_max_gt"]):
        reasons.append("raw count range indicates normalization is still required")
    platform_classification = "single_cell_spatial" if platform_name_lower in supported_platforms and str(metadata["observation_unit"]).lower() == "cell" else "spot_based_or_unsupported"

    query_string = build_reference_query_string(
        str(metadata["species"]),
        str(metadata["tissue"]),
        str(metadata["condition"]),
    )
    batch_values = [str(cell[selected_batch_column]) for cell in cells]
    return {
        "run_label": payload["run_label"],
        "dataset_id": metadata["dataset_id"],
        "platform_name": platform_name,
        "platform_classification": platform_classification,
        "n_obs": len(cells),
        "n_vars": len(genes),
        "obs_columns": obs_columns,
        "var_columns": var_columns,
        "coordinate_keys": coordinate_keys,
        "data_min": float(counts.min()),
        "data_max": float(counts.max()),
        "integer_like_counts": bool(np.allclose(counts, np.round(counts))),
        "normalize_applied": bool(counts.max() > float(payload["preprocessing"]["normalize_if_max_gt"])),
        "candidate_batch_columns": candidate_rows,
        "selected_batch_column": selected_batch_column,
        "sample_counts": dict(Counter(batch_values)),
        "species": metadata["species"],
        "tissue": metadata["tissue"],
        "condition": metadata["condition"],
        "query_string": query_string,
        "query_variants": build_reference_query_variants(
            str(metadata["species"]),
            str(metadata["tissue"]),
            str(metadata["condition"]),
        ),
        "platform_detection_reasons": reasons,
    }


def rank_reference_catalog(payload: dict[str, Any], query_string: str) -> pd.DataFrame:
    query_metadata = payload["query_metadata"]
    return rank_reference_candidates(
        payload["toy_czi_catalog"],
        {
            "species": str(query_metadata["species"]),
            "tissue": str(query_metadata["tissue"]),
            "condition": str(query_metadata["condition"]),
            "query_string": query_string,
            "min_cells": int(payload.get("reference_selection", {}).get("min_cells", 1000)),
        },
    )


def build_label_transfer(
    *,
    payload: dict[str, Any],
    query_cells: list[dict[str, Any]],
    query_matrix_norm: np.ndarray,
    query_latent: np.ndarray,
    reference_cells: list[dict[str, Any]],
    reference_latent: np.ndarray,
    shared_genes: list[str],
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, int]]:
    preprocessing = payload["preprocessing"]
    marker_panels = payload["marker_panels"]
    reference_ids = [cell["reference_cell_id"] for cell in reference_cells]
    reference_batches = [cell["donor_id"] for cell in reference_cells]
    broad_labels = [cell["broad_label"] for cell in reference_cells]
    subtype_labels = [cell["subtype_label"] for cell in reference_cells]
    unique_broad_labels = sorted(set(broad_labels))
    subtype_to_broad = {entry["label"]: entry["broad_label"] for entry in marker_panels["subtype"]}
    subtype_labels_by_broad: dict[str, list[str]] = {}
    for subtype, broad in subtype_to_broad.items():
        subtype_labels_by_broad.setdefault(broad, []).append(subtype)

    shared_gene_to_index = {gene: index for index, gene in enumerate(shared_genes)}
    broad_marker_map = {entry["label"]: entry["markers"] for entry in marker_panels["broad"]}
    subtype_marker_map = {entry["label"]: entry["markers"] for entry in marker_panels["subtype"]}
    distance_matrix = pairwise_euclidean(query_latent, reference_latent)
    threshold_config = payload["marker_thresholds"]

    label_rows: list[dict[str, Any]] = []
    annotation_seed_rows: list[dict[str, Any]] = []
    broad_counts: dict[str, int] = Counter()
    for query_index, query_cell in enumerate(query_cells):
        broad_vote = weighted_vote(
            distance_row=distance_matrix[query_index],
            reference_ids=reference_ids,
            labels=broad_labels,
            reference_batches=reference_batches,
            k=int(preprocessing["reference_k"]),
        )
        predicted_broad = str(broad_vote["predicted_label"])
        broad_counts[predicted_broad] += 1
        broad_margin = float(broad_vote["margin"])

        subtype_indices = [index for index, broad in enumerate(broad_labels) if broad == predicted_broad]
        subtype_vote = weighted_vote(
            distance_row=distance_matrix[query_index],
            reference_ids=reference_ids,
            labels=subtype_labels,
            reference_batches=reference_batches,
            k=int(preprocessing["reference_k"]),
            allowed_indices=subtype_indices,
        )
        predicted_subtype = str(subtype_vote["predicted_label"])
        subtype_margin = float(subtype_vote["margin"])
        label_confidence = round((broad_margin + subtype_margin) / 2.0, 6)

        broad_gap = score_marker_gap(
            row_vector=query_matrix_norm[query_index],
            predicted_label=predicted_broad,
            marker_map=broad_marker_map,
            shared_gene_to_index=shared_gene_to_index,
            competitor_labels=[label for label in unique_broad_labels if label != predicted_broad],
        )
        subtype_gap = score_marker_gap(
            row_vector=query_matrix_norm[query_index],
            predicted_label=predicted_subtype,
            marker_map=subtype_marker_map,
            shared_gene_to_index=shared_gene_to_index,
            competitor_labels=[
                label
                for label in subtype_labels_by_broad[predicted_broad]
                if label != predicted_subtype
            ],
        )

        if (
            broad_gap["gap"] >= float(threshold_config["broad_supported_min"])
            and subtype_gap["gap"] >= float(threshold_config["subtype_supported_min"])
        ):
            marker_status = "supported"
        elif (
            broad_gap["gap"] >= float(threshold_config["broad_borderline_min"])
            and subtype_gap["gap"] >= float(threshold_config["subtype_borderline_min"])
        ):
            marker_status = "borderline"
        else:
            marker_status = "conflict"

        notes = []
        if broad_margin <= float(threshold_config["broad_margin_warning_max"]):
            notes.append("low_broad_margin")
        if subtype_margin <= float(threshold_config["subtype_margin_warning_max"]):
            notes.append("low_subtype_margin")
        if marker_status != "supported":
            notes.append(f"marker_{marker_status}")
        label_rows.append(
            {
                "cell_id": query_cell["cell_id"],
                "sample_id": query_cell["sample_id"],
                "nearest_reference_ids": ";".join(subtype_vote["neighbor_ids"]),
                "nearest_reference_broad_labels": ";".join(broad_vote["neighbor_labels"]),
                "nearest_reference_subtypes": ";".join(subtype_vote["neighbor_labels"]),
                "nearest_reference_batches": ";".join(subtype_vote["neighbor_batches"]),
                "broad_supports": ranked_string(broad_vote["supports"]),
                "broad_margin": round(broad_margin, 6),
                "broad_label": predicted_broad,
                "subtype_supports": ranked_string(subtype_vote["supports"]),
                "subtype_margin": round(subtype_margin, 6),
                "predicted_label": predicted_subtype,
                "label_confidence": label_confidence,
                "broad_marker_gap": round(float(broad_gap["gap"]), 6),
                "subtype_marker_gap": round(float(subtype_gap["gap"]), 6),
                "marker_status": marker_status,
                "latent_1": round(float(query_latent[query_index, 0]), 6),
                "latent_2": round(float(query_latent[query_index, 1]) if query_latent.shape[1] > 1 else 0.0, 6),
                "latent_3": round(float(query_latent[query_index, 2]) if query_latent.shape[1] > 2 else 0.0, 6),
                "notes": ";".join(notes) if notes else "stable",
            }
        )
        annotation_seed_rows.append(
            {
                "cell_id": query_cell["cell_id"],
                "sample_id": query_cell["sample_id"],
                "x_coord": float(query_cell["x_coord"]),
                "y_coord": float(query_cell["y_coord"]),
                "predicted_label": predicted_subtype,
                "label_source": "harmony_surrogate_hierarchical_transfer",
                "label_confidence": label_confidence,
                "broad_label": predicted_broad,
                "marker_status": marker_status,
                "notes": notes,
            }
        )
    return (
        pd.DataFrame(label_rows, columns=LABEL_TRANSFER_COLUMNS),
        annotation_seed_rows,
        {key: int(value) for key, value in broad_counts.items()},
    )


def build_utag_niches(
    *,
    payload: dict[str, Any],
    annotation_seed_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, dict[str, str], dict[str, str], list[str]]:
    broad_labels = sorted({row["broad_label"] for row in annotation_seed_rows})
    broad_to_index = {label: index for index, label in enumerate(broad_labels)}
    utag_config = payload["utag"]
    prototype_labels = [item["niche_label"] for item in utag_config["niche_prototypes"]]
    prototype_matrix = np.asarray(
        [
            [float(item["composition"][broad_label]) for broad_label in broad_labels]
            for item in utag_config["niche_prototypes"]
        ],
        dtype=float,
    )
    prototype_matrix = row_normalize(prototype_matrix)
    boundary_targets = np.asarray(
        [float(item["boundary_target"]) for item in utag_config["niche_prototypes"]],
        dtype=float,
    )

    label_only_map = dict(utag_config["label_only_niche_map"])
    broad_probability = np.zeros((len(annotation_seed_rows), len(broad_labels)), dtype=float)
    for index, row in enumerate(annotation_seed_rows):
        broad_probability[index, broad_to_index[row["broad_label"]]] = 1.0

    rows = []
    cell_to_niche: dict[str, str] = {}
    cell_to_evidence: dict[str, str] = {}
    graph_shifted_cells: list[str] = []
    sample_groups: dict[str, list[int]] = {}
    for index, row in enumerate(annotation_seed_rows):
        sample_groups.setdefault(row["sample_id"], []).append(index)

    for sample_id, indices in sorted(sample_groups.items()):
        coordinates = np.asarray(
            [[annotation_seed_rows[index]["x_coord"], annotation_seed_rows[index]["y_coord"]] for index in indices],
            dtype=float,
        )
        distances = pairwise_euclidean(coordinates, coordinates)
        np.fill_diagonal(distances, np.inf)
        neighbor_k = min(int(utag_config["neighbor_k"]), len(indices) - 1)
        if neighbor_k < 1:
            raise AssertionError("UTAG surrogate requires at least two cells per sample.")
        neighbor_graph = np.zeros((len(indices), len(indices)), dtype=float)
        neighbor_id_lists: list[list[str]] = []
        neighbor_label_lists: list[list[str]] = []
        for source_local_index, source_global_index in enumerate(indices):
            order = np.argsort(distances[source_local_index])[:neighbor_k]
            weights = np.asarray(
                [1.0 / max(float(distances[source_local_index, target_local_index]), 1e-9) for target_local_index in order],
                dtype=float,
            )
            weights = weights / weights.sum()
            neighbor_graph[source_local_index, order] = weights
            neighbor_id_lists.append([annotation_seed_rows[indices[target_local_index]]["cell_id"] for target_local_index in order])
            neighbor_label_lists.append(
                [annotation_seed_rows[indices[target_local_index]]["broad_label"] for target_local_index in order]
            )

        local_profile = neighbor_graph @ broad_probability[indices]
        combined_profile = row_normalize(
            float(utag_config["self_weight"]) * broad_probability[indices]
            + (1.0 - float(utag_config["self_weight"])) * local_profile
        )
        smoothed_profile = row_normalize(
            (1.0 - float(utag_config["smoothing_weight"])) * combined_profile
            + float(utag_config["smoothing_weight"]) * (neighbor_graph @ combined_profile)
        )
        boundary_entropy = normalized_entropy(smoothed_profile)

        local_scores = combined_profile @ prototype_matrix.T
        boundary_alignment = 1.0 - np.abs(boundary_entropy[:, None] - boundary_targets[None, :])
        final_scores = local_scores + float(utag_config["boundary_weight"]) * boundary_alignment
        probabilities = softmax_rows(final_scores * float(utag_config["temperature"]))
        local_best_indices = local_scores.argmax(axis=1)
        final_best_indices = probabilities.argmax(axis=1)
        top_two = np.partition(probabilities, -2, axis=1)[:, -2:]
        confidence_margin = top_two[:, 1] - top_two[:, 0]

        for local_index, global_index in enumerate(indices):
            cell_id = annotation_seed_rows[global_index]["cell_id"]
            default_niche = label_only_map[annotation_seed_rows[global_index]["broad_label"]]
            local_best_niche = prototype_labels[int(local_best_indices[local_index])]
            final_niche = prototype_labels[int(final_best_indices[local_index])]
            graph_shifted = final_niche != default_niche
            if graph_shifted:
                graph_shifted_cells.append(cell_id)
            profile_string = ";".join(
                f"{broad_labels[broad_index]}={smoothed_profile[local_index, broad_index]:.6f}"
                for broad_index in range(len(broad_labels))
            )
            local_profile_string = ";".join(
                f"{broad_labels[broad_index]}={combined_profile[local_index, broad_index]:.6f}"
                for broad_index in range(len(broad_labels))
            )
            prototype_score_string = ";".join(
                f"{prototype_labels[prototype_index]}={probabilities[local_index, prototype_index]:.6f}"
                for prototype_index in range(len(prototype_labels))
            )
            reason = (
                f"default={default_niche}; local={local_best_niche}; boundary={boundary_entropy[local_index]:.3f}; "
                f"smoothed={profile_string}"
            )
            rows.append(
                {
                    "cell_id": cell_id,
                    "sample_id": sample_id,
                    "neighbor_ids": ";".join(neighbor_id_lists[local_index]),
                    "neighbor_broad_labels": ";".join(neighbor_label_lists[local_index]),
                    "neighbor_weight_sum": round(float(neighbor_graph[local_index].sum()), 6),
                    "local_broad_profile": local_profile_string,
                    "smoothed_broad_profile": profile_string,
                    "boundary_entropy": round(float(boundary_entropy[local_index]), 6),
                    "label_only_niche": default_niche,
                    "local_best_niche": local_best_niche,
                    "utag_niche": final_niche,
                    "graph_shifted": bool(graph_shifted),
                    "niche_confidence": round(float(confidence_margin[local_index]), 6),
                    "prototype_scores": prototype_score_string,
                    "niche_reason": reason,
                }
            )
            cell_to_niche[cell_id] = final_niche
            cell_to_evidence[cell_id] = reason

    return pd.DataFrame(rows, columns=UTAG_COLUMNS), cell_to_niche, cell_to_evidence, sorted(graph_shifted_cells)


def validate_outputs(skill_dir: Path, outdir: Path, *, expected: dict[str, Any] | None = None) -> dict[str, Any]:
    metadata = load_metadata(skill_dir)
    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        if deliverable["kind"] == "tsv":
            frame = pd.read_csv(path, sep="\t")
            missing = [column for column in deliverable.get("required_columns", []) if column not in frame.columns]
            if missing:
                raise AssertionError(f"Missing TSV columns in {path}: {missing}")
            if frame.empty:
                raise AssertionError(f"Empty TSV deliverable: {path}")
        elif deliverable["kind"] == "md":
            validate_markdown_sections(path, deliverable.get("required_sections", []))
        else:
            raise AssertionError(f"Unsupported deliverable kind in annotation starter: {deliverable['kind']}")

    starter_qc_files = metadata.get("starter_qc_files", QC_FILES)
    for filename in starter_qc_files:
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    annotation = pd.read_csv(outdir / "annotation_table.tsv", sep="\t").set_index("cell_id")
    label_transfer = pd.read_csv(outdir / "qc_label_transfer.tsv", sep="\t").set_index("cell_id")
    search = pd.read_csv(outdir / "qc_reference_search.tsv", sep="\t")
    summary = load_json(outdir / "run_summary.json")
    exploration = load_json(outdir / "qc_dataset_exploration.json")
    selection = load_json(outdir / "qc_reference_selection.json")

    if search.empty:
        raise AssertionError("Reference search ranking must not be empty.")
    required_search_columns = {
        "accept_for_transfer",
        "selection_notes",
        "rejection_reasons",
        "species_match",
        "tissue_match",
        "condition_score",
    }
    if not required_search_columns.issubset(set(search.columns)):
        raise AssertionError("qc_reference_search.tsv is missing reference-ranking decision columns.")
    if not bool(search.iloc[0]["selected"]):
        raise AssertionError("Top-ranked reference row must be marked as selected.")
    if not bool(search.iloc[0]["accept_for_transfer"]):
        raise AssertionError("Top-ranked reference row must pass the transfer gates.")
    if str(search.iloc[0]["dataset_id"]) != str(selection["selected_dataset_id"]):
        raise AssertionError("Selected reference must match the top-ranked catalog row.")
    if not all(column in label_transfer.columns for column in LABEL_TRANSFER_COLUMNS[1:]):
        raise AssertionError("qc_label_transfer.tsv is missing required transfer columns.")
    if annotation.shape[0] != label_transfer.shape[0]:
        raise AssertionError("Deliverable and QC row counts must agree.")
    if annotation.index.tolist() != label_transfer.index.tolist():
        raise AssertionError("Annotation and transfer QC must keep the same cell ordering.")
    if not (annotation["predicted_label"] == label_transfer["predicted_label"]).all():
        raise AssertionError("annotation_table.tsv predicted labels must match qc_label_transfer.tsv.")
    if not (annotation["broad_label"] == label_transfer["broad_label"]).all():
        raise AssertionError("annotation_table.tsv broad labels must match qc_label_transfer.tsv.")
    if set(summary["written_files"]) != set([item["path"] for item in metadata["deliverables"]] + starter_qc_files):
        raise AssertionError("run_summary.json written_files must list all deliverables and QC artifacts exactly.")
    if summary["method_steps"] != METHOD_STEPS:
        raise AssertionError("run_summary.json must record the exact starter method order.")
    if summary["selected_reference_dataset_id"] != selection["selected_dataset_id"]:
        raise AssertionError("run_summary selected reference must match selection QC.")
    if summary["platform_classification"] != exploration["platform_classification"]:
        raise AssertionError("run_summary platform classification must match exploration QC.")
    if exploration["platform_classification"] != "single_cell_spatial":
        raise AssertionError("Toy dataset must be classified as single_cell_spatial.")
    if "niche" in " ".join(annotation.columns).lower():
        raise AssertionError("annotation_table.tsv must not contain tissue-niche columns.")

    if expected is not None:
        if int(exploration["n_obs"]) != int(expected["query_cell_count"]):
            raise AssertionError("Unexpected query cell count.")
        if int(search.shape[0]) != int(expected["catalog_candidate_count"]):
            raise AssertionError("Unexpected reference catalog size.")
        if exploration["query_string"] != expected["query_string"]:
            raise AssertionError("Unexpected reference search query string.")
        if exploration["selected_batch_column"] != expected["selected_batch_column"]:
            raise AssertionError("Unexpected selected batch column.")
        if selection["selected_dataset_id"] != expected["selected_reference_dataset_id"]:
            raise AssertionError("Unexpected selected reference dataset.")
        if int(selection["shared_gene_count"]) != int(expected["shared_gene_count"]):
            raise AssertionError("Unexpected shared gene count.")
        if sorted(selection["reference_only_genes"]) != sorted(expected["reference_only_genes"]):
            raise AssertionError("Unexpected reference-only genes.")
        if summary["prediction_counts"] != expected["predicted_label_counts"]:
            raise AssertionError("Unexpected predicted label counts.")
        if summary["broad_label_counts"] != expected["broad_label_counts"]:
            raise AssertionError("Unexpected broad label counts.")
        if summary["marker_status_counts"] != expected["marker_status_counts"]:
            raise AssertionError("Unexpected marker status counts.")
        borderline_like = int((annotation["marker_status"] != "supported").sum())
        if borderline_like < int(expected["minimum_borderline_cells"]):
            raise AssertionError("Expected at least one non-supported marker review cell.")

    return {
        "annotation_rows": int(annotation.shape[0]),
        "review_flagged_cells": int((annotation["marker_status"] != "supported").sum()),
        "selected_reference_dataset_id": str(selection["selected_dataset_id"]),
    }


def build_annotation(payload: dict[str, Any], *, input_path: Path, outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)

    exploration = detect_platform_and_explore(payload)
    if exploration["platform_classification"] != "single_cell_spatial":
        raise AssertionError(
            "This annotation starter only supports single-cell-resolution spatial data. "
            "Route spot-based or unsupported inputs to the spatial-deconvolution workflow instead."
        )
    search_frame = rank_reference_catalog(payload, exploration["query_string"])
    selected_rows = search_frame.loc[search_frame["selected"]]
    if selected_rows.empty:
        search_summary = build_reference_selection_summary(
            search_frame,
            {
                "species": str(payload["query_metadata"]["species"]),
                "tissue": str(payload["query_metadata"]["tissue"]),
                "condition": str(payload["query_metadata"]["condition"]),
                "min_cells": int(payload.get("reference_selection", {}).get("min_cells", 1000)),
            },
        )
        raise AssertionError(str(search_summary["stop_reason"]))
    selected_catalog_row = selected_rows.iloc[0].to_dict()
    reference_payload = payload["reference_datasets"][str(selected_catalog_row["reference_key"])]

    query = payload["query"]
    query_cells = query["cells"]
    query_genes = query["genes"]
    reference_cells = reference_payload["cells"]
    reference_genes = reference_payload["genes"]
    shared_genes = build_shared_gene_order(query_genes, reference_genes)
    if not shared_genes:
        raise AssertionError("No shared genes between query and reference.")

    query_index = {gene: index for index, gene in enumerate(query_genes)}
    reference_index = {gene: index for index, gene in enumerate(reference_genes)}
    query_gene_indices = [query_index[gene] for gene in shared_genes]
    reference_gene_indices = [reference_index[gene] for gene in shared_genes]
    reference_only_genes = [gene for gene in reference_genes if gene not in shared_genes]
    query_only_genes = [gene for gene in query_genes if gene not in shared_genes]

    query_matrix_raw = matrix_from_rows(query_cells, "counts")[:, query_gene_indices]
    reference_matrix_raw = matrix_from_rows(reference_cells, "counts")[:, reference_gene_indices]
    target_sum = float(payload["preprocessing"]["target_sum"])
    query_matrix_norm, query_library_sizes = normalize_log1p(query_matrix_raw, target_sum)
    reference_matrix_norm, reference_library_sizes = normalize_log1p(reference_matrix_raw, target_sum)

    reference_latent, query_latent, explained_variance, _ = build_pca_surrogate(
        reference_matrix_norm,
        query_matrix_norm,
        int(payload["preprocessing"]["latent_dim"]),
    )
    reference_latent, reference_batch_shift = harmony_surrogate(
        reference_latent,
        [cell["donor_id"] for cell in reference_cells],
        float(payload["preprocessing"]["harmony_shrinkage"]),
    )
    query_latent, query_batch_shift = harmony_surrogate(
        query_latent,
        [cell[exploration["selected_batch_column"]] for cell in query_cells],
        float(payload["preprocessing"]["harmony_shrinkage"]),
    )

    label_transfer, annotation_seed_rows, broad_counts = build_label_transfer(
        payload=payload,
        query_cells=query_cells,
        query_matrix_norm=query_matrix_norm,
        query_latent=query_latent,
        reference_cells=reference_cells,
        reference_latent=reference_latent,
        shared_genes=shared_genes,
    )

    annotation_rows = [
        {
            "cell_id": row["cell_id"],
            "sample_id": row["sample_id"],
            "x_coord": row["x_coord"],
            "y_coord": row["y_coord"],
            "predicted_label": row["predicted_label"],
            "label_source": row["label_source"],
            "label_confidence": row["label_confidence"],
            "broad_label": row["broad_label"],
            "marker_status": row["marker_status"],
            "notes": ";".join(row["notes"]) if row["notes"] else "stable",
        }
        for row in annotation_seed_rows
    ]
    annotation_frame = pd.DataFrame(annotation_rows, columns=ANNOTATION_COLUMNS)

    exploration["shared_gene_count_after_selection"] = len(shared_genes)
    exploration["reference_only_genes"] = reference_only_genes
    exploration["query_only_genes"] = query_only_genes
    exploration["reference_cell_count"] = len(reference_cells)
    exploration["reference_library_sizes"] = [float(value) for value in reference_library_sizes.tolist()]
    exploration["query_library_sizes"] = [float(value) for value in query_library_sizes.tolist()]
    exploration["explained_variance"] = [round(float(value), 6) for value in explained_variance.tolist()]

    host_python_version = sys.version.split()[0]
    selection = {
        "query_string": exploration["query_string"],
        "query_variants": exploration["query_variants"],
        "selected_dataset_id": str(selected_catalog_row["dataset_id"]),
        "selected_dataset_title": str(selected_catalog_row["dataset_title"]),
        "reference_key": str(selected_catalog_row["reference_key"]),
        "dataset_h5ad_uri": str(selected_catalog_row["dataset_h5ad_uri"]),
        "selected_score": float(selected_catalog_row["score"]),
        "selected_notes": str(selected_catalog_row["selection_notes"]),
        "selected_rejection_reasons": str(selected_catalog_row["rejection_reasons"]),
        "shared_gene_count": len(shared_genes),
        "shared_genes": shared_genes,
        "reference_only_genes": reference_only_genes,
        "query_only_genes": query_only_genes,
        "materialized_reference_cell_count": len(reference_cells),
        "materialized_reference_batches": dict(Counter(cell["donor_id"] for cell in reference_cells)),
        "download_contract": {
            "dataset_lookup_call": "census['census_info']['datasets'].read().concat().to_pandas()",
            "download_call": "cellxgene_census.download_source_h5ad(dataset_id, to_path='reference.h5ad')",
            "slice_anchor": "cellxgene_census.get_anndata(...)",
            "discover_fallback": "Use the CELLxGENE Discover browser to export metadata and download the selected source H5AD when the Census Python API cannot be installed locally.",
        },
        "runtime_contract": {
            "host_python_version": host_python_version,
            "cellxgene_census_requires_python": ">=3.10,<3.13",
            "host_python_supported_for_census": sys.version_info[:2] >= (3, 10) and sys.version_info[:2] < (3, 13),
            "preferred_conda_command": "conda create -y -n sc-annotation-py312 python=3.12 && conda activate sc-annotation-py312 && pip install --upgrade pip && pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
        },
        "selection_reason": [
            "highest-scoring accepted candidate after exact species+tissue matching and healthy-reference preference",
            "broad and subtype labels available for hierarchical transfer",
            "source H5AD URI retained so a real run can materialize the same dataset in a supported Python 3.10-3.12 environment",
        ],
        "reference_batch_shift_norms": reference_batch_shift,
        "query_batch_shift_norms": query_batch_shift,
    }

    marker_counts = {key: int(value) for key, value in annotation_frame["marker_status"].value_counts().to_dict().items()}
    prediction_counts = {key: int(value) for key, value in annotation_frame["predicted_label"].value_counts().to_dict().items()}
    low_confidence_cells = annotation_frame.loc[annotation_frame["notes"].str.contains("low_"), "cell_id"].tolist()
    flagged_cells = annotation_frame.loc[annotation_frame["marker_status"] != "supported", "cell_id"].tolist()

    marker_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}.",
                f"Detected platform class: {exploration['platform_classification']} from {exploration['platform_name']}.",
                f"Query cells: {exploration['n_obs']}, genes: {exploration['n_vars']}, query string: {exploration['query_string']}.",
            ],
        },
        {
            "name": "Reference source",
            "bullets": [
                f"Selected toy reference: {selection['selected_dataset_id']} ({selection['selected_dataset_title']}).",
                f"Source H5AD URI contract: {selection['dataset_h5ad_uri']}.",
                "Real runs should either create a Python 3.10-3.12 environment for `cellxgene_census` or use the Discover metadata-export fallback before downloading the chosen source H5AD.",
            ],
        },
        {
            "name": "Marker review",
            "bullets": [
                f"Marker status counts: {format_counts(marker_counts)}.",
                f"Broad label counts: {format_counts(broad_counts)}.",
                f"Shared genes used for transfer: {', '.join(shared_genes)}.",
            ],
        },
        {
            "name": "Conflicts and resolutions",
            "bullets": [
                f"Low-confidence cells retained with explicit notes: {', '.join(low_confidence_cells) if low_confidence_cells else 'none'}.",
                f"Marker-flagged cells needing review: {', '.join(flagged_cells) if flagged_cells else 'none'}.",
                "Broad labels are assigned first from batch-corrected reference votes, then subtype labels are refined within the selected broad class.",
            ],
        },
        {
            "name": "Remaining uncertainty",
            "bullets": [
                "The toy reference is synthetic and only stands in for a real CELLxGENE atlas slice.",
                "Harmony is represented by a deterministic surrogate, so these outputs validate workflow shape rather than biological truth.",
                "This package stops at cell-type annotation. Tissue niche analysis should be done in the separate `tissue-niche-annotation` package.",
            ],
        },
    ]

    annotation_frame.to_csv(outdir / "annotation_table.tsv", sep="\t", index=False)
    label_transfer.to_csv(outdir / "qc_label_transfer.tsv", sep="\t", index=False)
    search_frame.to_csv(outdir / "qc_reference_search.tsv", sep="\t", index=False)
    write_json(outdir / "qc_dataset_exploration.json", exploration)
    write_json(outdir / "qc_reference_selection.json", selection)
    write_markdown(outdir / "marker_evidence.md", "Marker Evidence", marker_sections)

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path.resolve()),
        "outdir": str(outdir.resolve()),
        "method_steps": METHOD_STEPS,
        "platform_classification": exploration["platform_classification"],
        "query_string": exploration["query_string"],
        "selected_reference_dataset_id": selection["selected_dataset_id"],
        "cell_count": int(annotation_frame.shape[0]),
        "shared_gene_count": len(shared_genes),
        "prediction_counts": prediction_counts,
        "broad_label_counts": broad_counts,
        "marker_status_counts": marker_counts,
        "low_confidence_cells": low_confidence_cells,
        "review_flagged_cells": flagged_cells,
        "written_files": sorted([item["path"] for item in load_metadata(SKILL_DIR)["deliverables"]] + QC_FILES),
    }
    write_json(outdir / "run_summary.json", summary)
    validate_outputs(SKILL_DIR, outdir, expected=payload.get("expected_invariants"))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the spatial annotation experiment starter on the bundled toy data.")
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    result = build_annotation(payload, input_path=args.input.resolve(), outdir=args.outdir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
