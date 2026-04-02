#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


STRING_DTYPE = h5py.string_dtype(encoding="utf-8")
REPORT_TITLE = "Multimodal Integration Starter Report"
MODALITY_ORDER = ["rna", "protein", "chromatin"]
MODEL_WEIGHTS = {
    "multivi_style": {"rna": 0.45, "protein": 0.20, "chromatin": 0.35},
    "totalvi_style": {"rna": 0.45, "protein": 0.35, "chromatin": 0.20},
}
NORMALIZATION_NOTES = {
    "rna": "library-size normalize to 1e4 then log1p",
    "protein": "CLR-like log1p centering per cell",
    "chromatin": "TF-IDF-like scaling then log1p",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata(skill_dir: Path) -> dict[str, Any]:
    return load_json(skill_dir / "metadata.yaml")


def default_example_path(skill_dir: Path) -> Path:
    return skill_dir / "examples" / "toy_input.json"


def require_keys(mapping: dict[str, Any], required: list[str], *, name: str) -> None:
    missing = [key for key in required if key not in mapping]
    if missing:
        raise AssertionError(f"Missing keys in {name}: {missing}")


def to_string_array(values: list[str]) -> np.ndarray:
    return np.asarray(values, dtype=STRING_DTYPE)


def validate_input(payload: dict[str, Any]) -> None:
    require_keys(
        payload,
        [
            "run_label",
            "model_choice",
            "latent_dim",
            "knn_k",
            "rna_features",
            "protein_features",
            "chromatin_features",
            "cells",
            "expected_invariants",
        ],
        name="toy input",
    )
    if payload["model_choice"] not in MODEL_WEIGHTS:
        raise AssertionError(f"Unsupported model choice: {payload['model_choice']}")
    if int(payload["latent_dim"]) < 2:
        raise AssertionError("latent_dim must be at least 2 for this starter.")

    cells = payload["cells"]
    if len(cells) < 8:
        raise AssertionError("The multimodal starter expects at least eight cells.")
    cell_ids = [str(cell["cell_id"]) for cell in cells]
    if len(cell_ids) != len(set(cell_ids)):
        raise AssertionError("Cell IDs must be unique.")

    feature_lengths = {
        "rna": len(payload["rna_features"]),
        "protein": len(payload["protein_features"]),
        "chromatin": len(payload["chromatin_features"]),
    }
    batches: set[str] = set()
    labels: set[str] = set()
    splits: set[str] = set()
    for cell in cells:
        require_keys(cell, ["cell_id", "batch", "split", "true_label", "rna", "protein", "chromatin"], name="cell")
        batches.add(str(cell["batch"]))
        labels.add(str(cell["true_label"]))
        splits.add(str(cell["split"]))
        for modality, feature_count in feature_lengths.items():
            values = cell[modality]
            if len(values) != feature_count:
                raise AssertionError(
                    f"Cell {cell['cell_id']} has {len(values)} {modality} values, expected {feature_count}."
                )
            if any(float(value) < 0 for value in values):
                raise AssertionError(f"Cell {cell['cell_id']} has negative {modality} values.")

    if len(batches) < 2:
        raise AssertionError("The starter expects at least two batches.")
    if len(labels) < 2:
        raise AssertionError("The starter expects at least two labels.")
    if not {"reference", "query"}.issubset(splits):
        raise AssertionError("The starter expects both reference and query splits.")

    expected = payload["expected_invariants"]
    require_keys(expected, ["query_labels", "min_batch_mixing", "min_query_accuracy"], name="expected_invariants")


def matrix_from_cells(cells: list[dict[str, Any]], modality: str) -> np.ndarray:
    return np.asarray([cell[modality] for cell in cells], dtype=float)


def library_normalize_log1p(matrix: np.ndarray, *, target_sum: float) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = matrix.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise AssertionError("Every row must have a positive library size.")
    normalized = (matrix / library_sizes[:, None]) * float(target_sum)
    return np.log1p(normalized), library_sizes


def clr_normalize(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = matrix.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise AssertionError("Protein rows must have positive totals.")
    logged = np.log1p(matrix)
    centered = logged - logged.mean(axis=1, keepdims=True)
    return centered, library_sizes


def tfidf_log1p(matrix: np.ndarray, *, scale: float) -> tuple[np.ndarray, np.ndarray]:
    totals = matrix.sum(axis=1)
    if np.any(totals <= 0):
        raise AssertionError("Chromatin rows must have positive totals.")
    tf = matrix / totals[:, None]
    present = (matrix > 0).sum(axis=0)
    idf = np.log1p(matrix.shape[0] / (1.0 + present))
    transformed = np.log1p(tf * idf[None, :] * float(scale))
    return transformed, totals


def zscore_columns(matrix: np.ndarray) -> np.ndarray:
    means = matrix.mean(axis=0, keepdims=True)
    stds = matrix.std(axis=0, keepdims=True)
    stds[stds < 1e-8] = 1.0
    return (matrix - means) / stds


def compute_weighted_latent(
    normalized: dict[str, np.ndarray],
    weights: dict[str, float],
    *,
    latent_dim: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    blocks = [math.sqrt(weights[modality]) * zscore_columns(normalized[modality]) for modality in MODALITY_ORDER]
    joint = np.hstack(blocks)
    u, singular_values, _ = np.linalg.svd(joint, full_matrices=False)
    latent = u[:, :latent_dim] * singular_values[:latent_dim]
    for dim in range(latent.shape[1]):
        if latent[0, dim] < 0:
            latent[:, dim] *= -1.0
    variance = singular_values**2
    explained = variance[:latent_dim] / variance.sum()
    return latent, joint, explained


def pairwise_distances(matrix: np.ndarray) -> np.ndarray:
    deltas = matrix[:, None, :] - matrix[None, :, :]
    return np.sqrt(np.square(deltas).sum(axis=2))


def build_knn_graph(latent: np.ndarray, *, k: int) -> tuple[np.ndarray, np.ndarray]:
    distances = pairwise_distances(latent)
    np.fill_diagonal(distances, np.inf)
    indices = np.argsort(distances, axis=1)[:, :k]
    neighbor_distances = np.take_along_axis(distances, indices, axis=1)
    return indices, neighbor_distances


def average_centroid_distance(matrix: np.ndarray, groups: list[str]) -> float:
    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        return 0.0
    centroids = [matrix[np.asarray(groups) == group].mean(axis=0) for group in unique_groups]
    pair_distances: list[float] = []
    for left in range(len(centroids)):
        for right in range(left + 1, len(centroids)):
            pair_distances.append(float(np.linalg.norm(centroids[left] - centroids[right])))
    return float(np.mean(pair_distances))


def inverse_distance_vote(
    latent: np.ndarray,
    *,
    cell_ids: list[str],
    true_labels: list[str],
    splits: list[str],
    k: int,
) -> tuple[list[str], np.ndarray, list[dict[str, Any]]]:
    reference_indices = [index for index, split in enumerate(splits) if split == "reference"]
    query_indices = [index for index, split in enumerate(splits) if split == "query"]
    if len(reference_indices) < k:
        raise AssertionError("Reference set must have at least k cells for label transfer.")

    predicted_labels = list(true_labels)
    confidences = np.ones(len(cell_ids), dtype=float)
    query_rows: list[dict[str, Any]] = []
    for query_index in query_indices:
        distances = []
        for reference_index in reference_indices:
            distance = float(np.linalg.norm(latent[query_index] - latent[reference_index]))
            distances.append((reference_index, distance))
        distances.sort(key=lambda item: (item[1], cell_ids[item[0]]))
        top_neighbors = distances[:k]
        label_weights: dict[str, float] = {}
        for reference_index, distance in top_neighbors:
            weight = 1.0 / max(distance, 1e-6)
            label = true_labels[reference_index]
            label_weights[label] = label_weights.get(label, 0.0) + weight
        ranked = sorted(label_weights.items(), key=lambda item: (-item[1], item[0]))
        predicted = ranked[0][0]
        confidence = ranked[0][1] / sum(label_weights.values())
        predicted_labels[query_index] = predicted
        confidences[query_index] = confidence
        query_rows.append(
            {
                "cell_id": cell_ids[query_index],
                "predicted_label": predicted,
                "true_label": true_labels[query_index],
                "confidence": round(float(confidence), 6),
                "neighbor_ids": [cell_ids[index] for index, _ in top_neighbors],
                "neighbor_labels": [true_labels[index] for index, _ in top_neighbors],
                "neighbor_distances": [round(float(distance), 6) for _, distance in top_neighbors],
            }
        )
    return predicted_labels, confidences, query_rows


def build_modality_qc_rows(
    *,
    normalized: dict[str, np.ndarray],
    raw: dict[str, np.ndarray],
    batches: list[str],
    labels: list[str],
    knn_indices: np.ndarray,
    predicted_labels: list[str],
    query_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for modality in MODALITY_ORDER:
        rows.append(
            {
                "metric": "mean_library_size",
                "modality": modality,
                "value": round(float(raw[modality].sum(axis=1).mean()), 6),
            }
        )
        rows.append(
            {
                "metric": "mean_nonzero_fraction",
                "modality": modality,
                "value": round(float((raw[modality] > 0).mean()), 6),
            }
        )
        rows.append(
            {
                "metric": "between_label_centroid_distance",
                "modality": modality,
                "value": round(average_centroid_distance(normalized[modality], labels), 6),
            }
        )
        rows.append(
            {
                "metric": "between_batch_centroid_distance",
                "modality": modality,
                "value": round(average_centroid_distance(normalized[modality], batches), 6),
            }
        )

    batch_mixing_values: list[float] = []
    label_consistency_values: list[float] = []
    for row_index in range(knn_indices.shape[0]):
        neighbors = knn_indices[row_index]
        batch_mixing_values.append(
            float(np.mean([batches[neighbor] != batches[row_index] for neighbor in neighbors]))
        )
        label_consistency_values.append(
            float(np.mean([predicted_labels[neighbor] == predicted_labels[row_index] for neighbor in neighbors]))
        )

    accuracy = float(
        np.mean([row["predicted_label"] == row["true_label"] for row in query_rows])
    )
    mean_confidence = float(np.mean([row["confidence"] for row in query_rows]))
    rows.extend(
        [
            {
                "metric": "batch_mixing",
                "modality": "joint",
                "value": round(float(np.mean(batch_mixing_values)), 6),
            },
            {
                "metric": "label_consistency",
                "modality": "joint",
                "value": round(float(np.mean(label_consistency_values)), 6),
            },
            {
                "metric": "query_label_accuracy",
                "modality": "joint",
                "value": round(accuracy, 6),
            },
            {
                "metric": "query_confidence_mean",
                "modality": "joint",
                "value": round(mean_confidence, 6),
            },
        ]
    )
    return rows


def write_string_dataset(group: h5py.Group, name: str, values: list[str]) -> None:
    group.create_dataset(name, data=to_string_array(values))


def write_modality_group(
    group: h5py.Group,
    *,
    cell_ids: list[str],
    features: list[str],
    raw_matrix: np.ndarray,
    normalized_matrix: np.ndarray,
    normalization_note: str,
) -> None:
    write_string_dataset(group, "obs_names", cell_ids)
    write_string_dataset(group, "var_names", features)
    group.create_dataset("raw_counts", data=raw_matrix)
    group.create_dataset("normalized", data=normalized_matrix)
    group.attrs["normalization"] = normalization_note


def write_h5mu_surrogate(
    path: Path,
    *,
    payload: dict[str, Any],
    cell_ids: list[str],
    batches: list[str],
    splits: list[str],
    true_labels: list[str],
    predicted_labels: list[str],
    confidences: np.ndarray,
    raw: dict[str, np.ndarray],
    normalized: dict[str, np.ndarray],
    latent: np.ndarray,
    explained: np.ndarray,
    knn_indices: np.ndarray,
    knn_distances: np.ndarray,
    query_rows: list[dict[str, Any]],
    qc_rows: list[dict[str, Any]],
) -> None:
    with h5py.File(path, "w") as handle:
        handle.attrs["format"] = "mudata_like_h5mu_surrogate"
        handle.attrs["layout_note"] = "Portable starter layout inspired by MuData; not a full mudata serialization."
        handle.attrs["model_choice"] = payload["model_choice"]

        mod_group = handle.create_group("mod")
        mod_group.attrs["mod-order"] = to_string_array(MODALITY_ORDER)
        write_modality_group(
            mod_group.create_group("rna"),
            cell_ids=cell_ids,
            features=payload["rna_features"],
            raw_matrix=raw["rna"],
            normalized_matrix=normalized["rna"],
            normalization_note=NORMALIZATION_NOTES["rna"],
        )
        write_modality_group(
            mod_group.create_group("protein"),
            cell_ids=cell_ids,
            features=payload["protein_features"],
            raw_matrix=raw["protein"],
            normalized_matrix=normalized["protein"],
            normalization_note=NORMALIZATION_NOTES["protein"],
        )
        write_modality_group(
            mod_group.create_group("chromatin"),
            cell_ids=cell_ids,
            features=payload["chromatin_features"],
            raw_matrix=raw["chromatin"],
            normalized_matrix=normalized["chromatin"],
            normalization_note=NORMALIZATION_NOTES["chromatin"],
        )

        obs_group = handle.create_group("obs")
        write_string_dataset(obs_group, "cell_id", cell_ids)
        write_string_dataset(obs_group, "batch", batches)
        write_string_dataset(obs_group, "split", splits)
        write_string_dataset(obs_group, "true_label", true_labels)
        write_string_dataset(obs_group, "predicted_label", predicted_labels)
        obs_group.create_dataset("label_transfer_confidence", data=confidences)

        obsm_group = handle.create_group("obsm")
        obsm_group.create_dataset("X_latent", data=latent)
        obsm_group.create_dataset("variance_explained", data=explained)

        obsp_group = handle.create_group("obsp")
        obsp_group.create_dataset("knn_indices", data=knn_indices)
        obsp_group.create_dataset("knn_distances", data=knn_distances)

        label_transfer_group = handle.create_group("label_transfer")
        write_string_dataset(label_transfer_group, "query_ids", [row["cell_id"] for row in query_rows])
        write_string_dataset(label_transfer_group, "predicted_label", [row["predicted_label"] for row in query_rows])
        write_string_dataset(label_transfer_group, "true_label", [row["true_label"] for row in query_rows])
        label_transfer_group.create_dataset(
            "confidence",
            data=np.asarray([row["confidence"] for row in query_rows], dtype=float),
        )
        write_string_dataset(
            label_transfer_group,
            "neighbor_labels",
            [";".join(row["neighbor_labels"]) for row in query_rows],
        )
        write_string_dataset(
            label_transfer_group,
            "neighbor_ids",
            [";".join(row["neighbor_ids"]) for row in query_rows],
        )

        qc_group = handle.create_group("qc")
        write_string_dataset(qc_group, "metric", [row["metric"] for row in qc_rows])
        write_string_dataset(qc_group, "modality", [row["modality"] for row in qc_rows])
        qc_group.create_dataset("value", data=np.asarray([row["value"] for row in qc_rows], dtype=float))

        uns_group = handle.create_group("uns")
        write_string_dataset(uns_group, "normalization_order", MODALITY_ORDER)


def write_modality_qc(path: Path, rows: list[dict[str, Any]]) -> None:
    frame = pd.DataFrame(rows, columns=["metric", "modality", "value"])
    frame.to_csv(path, sep="\t", index=False)


def write_report(
    path: Path,
    *,
    payload: dict[str, Any],
    qc_rows: list[dict[str, Any]],
    explained: np.ndarray,
    query_rows: list[dict[str, Any]],
) -> None:
    qc_lookup = {(row["metric"], row["modality"]): row["value"] for row in qc_rows}
    lines = [f"# {REPORT_TITLE}", ""]
    sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}",
                f"Model choice: {payload['model_choice']}",
                f"Cells: {len(payload['cells'])} total with {sum(cell['split'] == 'reference' for cell in payload['cells'])} reference and {sum(cell['split'] == 'query' for cell in payload['cells'])} query cells.",
            ],
        },
        {
            "name": "Modalities",
            "bullets": [
                f"RNA: {len(payload['rna_features'])} genes with {NORMALIZATION_NOTES['rna']}.",
                f"Protein: {len(payload['protein_features'])} antibodies with {NORMALIZATION_NOTES['protein']}.",
                f"Chromatin: {len(payload['chromatin_features'])} peaks with {NORMALIZATION_NOTES['chromatin']}.",
            ],
        },
        {
            "name": "Latent representation",
            "bullets": [
                f"Built a {int(payload['latent_dim'])}-dimensional weighted SVD latent that approximates {'MultiVI' if payload['model_choice'] == 'multivi_style' else 'totalVI'} emphasis without neural training.",
                f"Variance explained proxy: component_1={explained[0]:.6f}, component_2={explained[1]:.6f}.",
                f"Latent QC: batch_mixing={qc_lookup[('batch_mixing', 'joint')]:.6f}, label_consistency={qc_lookup[('label_consistency', 'joint')]:.6f}.",
            ],
        },
        {
            "name": "Label transfer",
            "bullets": [
                f"Reference-to-query kNN vote accuracy: {qc_lookup[('query_label_accuracy', 'joint')]:.6f}.",
                f"Mean query confidence: {qc_lookup[('query_confidence_mean', 'joint')]:.6f}.",
                *[
                    f"{row['cell_id']} -> {row['predicted_label']} (true={row['true_label']}, confidence={row['confidence']:.6f})"
                    for row in query_rows
                ],
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "Computed directly in this starter: modality-wise normalization, weighted latent embedding, kNN graph, label transfer, and QC metrics.",
                "Approximated: totalVI and MultiVI training are replaced by deterministic weighted fusion over tiny synthetic inputs.",
                "Surrogate boundary: integrated_latent.h5mu is a MuData-inspired HDF5 layout, not a full mudata serialization written by muon.",
            ],
        },
    ]
    for section in sections:
        lines.append(f"## {section['name']}")
        lines.append("")
        for bullet in section["bullets"]:
            lines.append(f"- {bullet}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def validate_markdown_sections(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def validate_h5mu_structure(path: Path) -> None:
    with h5py.File(path, "r") as handle:
        fmt = handle.attrs.get("format", "")
        if isinstance(fmt, bytes):
            fmt = fmt.decode("utf-8")
        if fmt != "mudata_like_h5mu_surrogate":
            raise AssertionError(f"Unexpected h5mu surrogate format in {path}: {fmt}")
        for group_name in ["mod", "obs", "obsm", "obsp", "label_transfer", "qc", "uns"]:
            if group_name not in handle:
                raise AssertionError(f"Missing group {group_name} in {path}")
        latent = handle["obsm"]["X_latent"][()]
        if latent.shape[0] == 0 or latent.shape[1] < 2:
            raise AssertionError(f"Latent space is unexpectedly small in {path}")
        for modality in MODALITY_ORDER:
            group = handle["mod"][modality]
            raw_counts = group["raw_counts"][()]
            normalized = group["normalized"][()]
            if raw_counts.shape != normalized.shape:
                raise AssertionError(f"Mismatched raw/normalized shapes for {modality} in {path}")
            if raw_counts.shape[0] != latent.shape[0]:
                raise AssertionError(f"Observation mismatch for {modality} in {path}")
        knn_indices = handle["obsp"]["knn_indices"][()]
        knn_distances = handle["obsp"]["knn_distances"][()]
        if knn_indices.shape != knn_distances.shape:
            raise AssertionError(f"kNN graph shape mismatch in {path}")


def validate_outputs(skill_dir: Path, outdir: Path) -> None:
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
                raise AssertionError(f"TSV deliverable is empty: {path}")
        elif deliverable["kind"] == "md":
            validate_markdown_sections(path, deliverable.get("required_sections", []))
        elif deliverable["kind"] == "h5mu":
            validate_h5mu_structure(path)
        else:
            raise AssertionError(f"Unsupported deliverable kind: {deliverable['kind']}")

    summary_path = outdir / "run_summary.json"
    if not summary_path.exists():
        raise AssertionError(f"Missing run summary: {summary_path}")


def run_starter(skill_dir: Path, *, input_path: Path, outdir: Path) -> dict[str, Any]:
    payload = load_json(input_path)
    validate_input(payload)
    outdir.mkdir(parents=True, exist_ok=True)

    cells = payload["cells"]
    cell_ids = [str(cell["cell_id"]) for cell in cells]
    batches = [str(cell["batch"]) for cell in cells]
    splits = [str(cell["split"]) for cell in cells]
    true_labels = [str(cell["true_label"]) for cell in cells]

    raw = {modality: matrix_from_cells(cells, modality) for modality in MODALITY_ORDER}
    normalized_rna, _ = library_normalize_log1p(raw["rna"], target_sum=10000.0)
    normalized_protein, _ = clr_normalize(raw["protein"])
    normalized_chromatin, _ = tfidf_log1p(raw["chromatin"], scale=100.0)
    normalized = {
        "rna": normalized_rna,
        "protein": normalized_protein,
        "chromatin": normalized_chromatin,
    }

    weights = MODEL_WEIGHTS[payload["model_choice"]]
    latent, _, explained = compute_weighted_latent(
        normalized,
        weights,
        latent_dim=int(payload["latent_dim"]),
    )
    knn_indices, knn_distances = build_knn_graph(latent, k=int(payload["knn_k"]))
    predicted_labels, confidences, query_rows = inverse_distance_vote(
        latent,
        cell_ids=cell_ids,
        true_labels=true_labels,
        splits=splits,
        k=int(payload["knn_k"]),
    )
    qc_rows = build_modality_qc_rows(
        normalized=normalized,
        raw=raw,
        batches=batches,
        labels=true_labels,
        knn_indices=knn_indices,
        predicted_labels=predicted_labels,
        query_rows=query_rows,
    )

    write_h5mu_surrogate(
        outdir / "integrated_latent.h5mu",
        payload=payload,
        cell_ids=cell_ids,
        batches=batches,
        splits=splits,
        true_labels=true_labels,
        predicted_labels=predicted_labels,
        confidences=confidences,
        raw=raw,
        normalized=normalized,
        latent=latent,
        explained=explained,
        knn_indices=knn_indices,
        knn_distances=knn_distances,
        query_rows=query_rows,
        qc_rows=qc_rows,
    )
    write_modality_qc(outdir / "modality_qc.tsv", qc_rows)
    write_report(
        outdir / "integration_report.md",
        payload=payload,
        qc_rows=qc_rows,
        explained=explained,
        query_rows=query_rows,
    )

    summary = {
        "run_label": payload["run_label"],
        "model_choice": payload["model_choice"],
        "counts": {
            "cells": len(cells),
            "reference_cells": splits.count("reference"),
            "query_cells": splits.count("query"),
            "rna_features": len(payload["rna_features"]),
            "protein_features": len(payload["protein_features"]),
            "chromatin_features": len(payload["chromatin_features"]),
        },
        "qc_metrics": {f"{row['metric']}::{row['modality']}": row["value"] for row in qc_rows},
        "query_predictions": [
            {
                "cell_id": row["cell_id"],
                "predicted_label": row["predicted_label"],
                "true_label": row["true_label"],
                "confidence": row["confidence"],
            }
            for row in query_rows
        ],
    }
    summary["written_files"] = sorted(
        [str(path.relative_to(outdir)) for path in outdir.rglob("*") if path.is_file()] + ["run_summary.json"]
    )
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    validate_outputs(skill_dir, outdir)
    return summary


def run_cli(skill_dir: Path, argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"Run the multimodal integration starter for {skill_dir.name}.")
    parser.add_argument("--input", type=Path, default=default_example_path(skill_dir))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = run_starter(skill_dir, input_path=args.input, outdir=args.outdir.resolve())
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0
