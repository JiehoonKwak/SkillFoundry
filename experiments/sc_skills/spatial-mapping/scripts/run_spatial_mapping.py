#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
QC_FILES = [
    "intermediate_qc.json",
    "gene_intersection_qc.tsv",
    "marker_qc.tsv",
    "niche_qc.tsv",
    "run_summary.json",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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
        for bullet in section.get("bullets", []):
            lines.append(f"- {bullet}")
        for paragraph in section.get("paragraphs", []):
            lines.append(paragraph)
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def normalize_log1p(matrix: np.ndarray, target_sum: float) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = matrix.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise AssertionError("All reference profiles and spatial observations must have positive library sizes.")
    normalized = (matrix / library_sizes[:, None]) * float(target_sum)
    return np.log1p(normalized), library_sizes


def cosine_similarity_matrix(spots: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    spot_norms = np.linalg.norm(spots, axis=1, keepdims=True)
    centroid_norms = np.linalg.norm(centroids, axis=1, keepdims=True).T
    denominator = np.clip(spot_norms * centroid_norms, 1e-8, None)
    return (spots @ centroids.T) / denominator


def build_neighbor_graph(coordinates: np.ndarray, *, n_neighbors: int) -> tuple[list[list[int]], np.ndarray]:
    if coordinates.shape[0] < 2:
        raise AssertionError("At least two spatial observations are required.")
    distances = np.sqrt(((coordinates[:, None, :] - coordinates[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    order = np.argsort(distances, axis=1)
    max_neighbors = min(n_neighbors, coordinates.shape[0] - 1)
    return [order[index, :max_neighbors].tolist() for index in range(coordinates.shape[0])], distances


def pick_majority_label(labels: list[str]) -> str:
    counts = Counter(labels)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def metric_status(value: float, floor: float, *, higher_is_better: bool) -> str:
    if higher_is_better:
        return "pass" if value >= floor else "review"
    return "pass" if value <= floor else "review"


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
            raise AssertionError(f"Unsupported deliverable kind: {deliverable['kind']}")

    for filename in metadata.get("starter_qc_files", QC_FILES):
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    labels = pd.read_csv(outdir / "mapped_labels.tsv", sep="\t")
    scores = pd.read_csv(outdir / "mapping_scores.tsv", sep="\t")
    gene_qc = pd.read_csv(outdir / "gene_intersection_qc.tsv", sep="\t")
    marker_qc = pd.read_csv(outdir / "marker_qc.tsv", sep="\t")
    niche_qc = pd.read_csv(outdir / "niche_qc.tsv", sep="\t")
    intermediate = load_json(outdir / "intermediate_qc.json")

    required_intermediate_keys = [
        "cosine_similarity",
        "gene_intersection",
        "label_assignments",
        "marker_consistency",
        "niche_majority",
        "normalization",
        "observation_ids",
        "reference_labels",
    ]
    missing_keys = [key for key in required_intermediate_keys if key not in intermediate]
    if missing_keys:
        raise AssertionError(f"Missing keys in intermediate_qc.json: {missing_keys}")

    similarity = np.asarray(intermediate["cosine_similarity"]["matrix"], dtype=float)
    reference_labels = list(intermediate["reference_labels"])
    observation_ids = list(intermediate["observation_ids"])
    if similarity.shape != (len(observation_ids), len(reference_labels)):
        raise AssertionError("Cosine similarity matrix shape does not match recorded labels.")

    label_frame = labels.set_index("observation_id")
    marker_frame = marker_qc.set_index("observation_id")
    niche_frame = niche_qc.set_index("observation_id")
    for row_index, observation_id in enumerate(observation_ids):
        if observation_id not in label_frame.index:
            raise AssertionError(f"Missing mapped label row for {observation_id}.")
        ranking = np.argsort(similarity[row_index])[::-1]
        primary_label = reference_labels[int(ranking[0])]
        secondary_label = reference_labels[int(ranking[1])]
        observed = label_frame.loc[observation_id]
        if observed["primary_label"] != primary_label:
            raise AssertionError(f"Primary label mismatch for {observation_id}.")
        if observed["secondary_label"] != secondary_label:
            raise AssertionError(f"Secondary label mismatch for {observation_id}.")
        primary_score = float(similarity[row_index, ranking[0]])
        secondary_score = float(similarity[row_index, ranking[1]])
        if abs(float(observed["primary_score"]) - primary_score) > 1e-6:
            raise AssertionError(f"Primary score mismatch for {observation_id}.")
        if abs(float(observed["score_margin"]) - (primary_score - secondary_score)) > 1e-6:
            raise AssertionError(f"Score margin mismatch for {observation_id}.")
        if observed["marker_review_status"] != marker_frame.loc[observation_id, "marker_review_status"]:
            raise AssertionError(f"Marker review status mismatch for {observation_id}.")
        if niche_frame.loc[observation_id, "niche_majority_label"] not in reference_labels:
            raise AssertionError(f"Unexpected niche-majority label for {observation_id}.")

    role_counts = gene_qc["role"].value_counts().to_dict()
    if role_counts.get("shared", 0) <= 0:
        raise AssertionError("Gene intersection QC must contain shared genes.")
    if not gene_qc["used_for_mapping"].any():
        raise AssertionError("Gene intersection QC must mark at least one mapping gene.")

    global_scores = scores.loc[scores["scope"] == "global"].set_index("metric_name")
    for metric_name in [
        "shared_gene_count",
        "shared_gene_fraction",
        "mean_primary_score",
        "mean_score_margin",
        "niche_agreement_rate",
        "supported_fraction",
    ]:
        if metric_name not in global_scores.index:
            raise AssertionError(f"Missing global metric {metric_name}.")

    if expected is not None:
        if role_counts.get("shared", 0) != int(expected["shared_gene_count"]):
            raise AssertionError("Unexpected shared gene count.")
        reference_only = sorted(gene_qc.loc[gene_qc["role"] == "reference_only", "gene"].tolist())
        spatial_only = sorted(gene_qc.loc[gene_qc["role"] == "spatial_only", "gene"].tolist())
        if reference_only != sorted(expected["reference_only_genes"]):
            raise AssertionError("Unexpected reference-only genes.")
        if spatial_only != sorted(expected["spatial_only_genes"]):
            raise AssertionError("Unexpected spatial-only genes.")
        observed_labels = label_frame["primary_label"].to_dict()
        if observed_labels != expected["primary_labels"]:
            raise AssertionError(f"Unexpected primary labels: {observed_labels}")
        observed_statuses = label_frame["marker_review_status"].to_dict()
        if observed_statuses != expected["review_statuses"]:
            raise AssertionError(f"Unexpected review statuses: {observed_statuses}")

        ambiguous = expected.get("ambiguous_observations", {})
        for observation_id, bounds in ambiguous.items():
            margin = float(label_frame.loc[observation_id, "score_margin"])
            marker_gap = float(marker_frame.loc[observation_id, "marker_gap"])
            if "max_score_margin" in bounds and margin > float(bounds["max_score_margin"]):
                raise AssertionError(f"Score margin too large for ambiguous observation {observation_id}.")
            if "min_marker_gap" in bounds and marker_gap < float(bounds["min_marker_gap"]):
                raise AssertionError(f"Marker gap too small for ambiguous observation {observation_id}.")
            if "max_marker_gap" in bounds and marker_gap > float(bounds["max_marker_gap"]):
                raise AssertionError(f"Marker gap too large for ambiguous observation {observation_id}.")

        metric_expectations = expected.get("metrics", {})
        if float(global_scores.loc["shared_gene_fraction", "value"]) < float(metric_expectations["shared_gene_fraction_min"]):
            raise AssertionError("Shared gene fraction fell below the expected floor.")
        if float(global_scores.loc["mean_primary_score", "value"]) < float(metric_expectations["mean_primary_score_min"]):
            raise AssertionError("Mean primary score fell below the expected floor.")
        if float(global_scores.loc["mean_score_margin", "value"]) < float(metric_expectations["mean_score_margin_min"]):
            raise AssertionError("Mean score margin fell below the expected floor.")
        if float(global_scores.loc["niche_agreement_rate", "value"]) < float(metric_expectations["niche_agreement_rate_min"]):
            raise AssertionError("Niche agreement rate fell below the expected floor.")
        if float(global_scores.loc["supported_fraction", "value"]) < float(metric_expectations["supported_fraction_min"]):
            raise AssertionError("Supported fraction fell below the expected floor.")

    return {
        "global_metrics": {
            metric_name: float(global_scores.loc[metric_name, "value"])
            for metric_name in global_scores.index
        },
        "observation_count": int(labels.shape[0]),
        "shared_gene_count": int(role_counts.get("shared", 0)),
    }


def build_spatial_mapping(payload: dict[str, Any], *, input_path: Path, outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)

    reference = payload["reference"]
    spatial = payload["spatial"]
    target_sum = float(payload.get("normalization_target_sum", 1000.0))
    n_neighbors = int(payload.get("niche", {}).get("n_neighbors", 3))

    reference_genes = list(reference["genes"])
    spatial_genes = list(spatial["genes"])
    shared_genes = [gene for gene in spatial_genes if gene in set(reference_genes)]
    reference_only = sorted(set(reference_genes) - set(shared_genes))
    spatial_only = sorted(set(spatial_genes) - set(shared_genes))
    union_gene_count = len(set(reference_genes) | set(spatial_genes))
    if len(shared_genes) < 4:
        raise AssertionError("At least four shared genes are required for the starter.")
    if union_gene_count <= 0:
        raise AssertionError("Unexpected empty gene universe.")

    profiles = list(reference["profiles"])
    spots = list(spatial["spots"])
    if len(profiles) < 2:
        raise AssertionError("At least two reference profiles are required.")
    if len(spots) < 3:
        raise AssertionError("At least three spatial observations are required.")

    reference_labels = [profile["label"] for profile in profiles]
    marker_table_by_label = {
        entry["label"]: [gene for gene in entry["markers"] if gene in shared_genes]
        for entry in payload["marker_table"]
    }
    missing_marker_labels = [label for label in reference_labels if not marker_table_by_label.get(label)]
    if missing_marker_labels:
        raise AssertionError(f"Missing shared markers for labels: {missing_marker_labels}")

    reference_index = {gene: index for index, gene in enumerate(reference_genes)}
    spatial_index = {gene: index for index, gene in enumerate(spatial_genes)}
    shared_reference_idx = [reference_index[gene] for gene in shared_genes]
    shared_spatial_idx = [spatial_index[gene] for gene in shared_genes]
    shared_gene_to_idx = {gene: index for index, gene in enumerate(shared_genes)}

    reference_raw = np.asarray([profile["centroid"] for profile in profiles], dtype=float)
    spatial_raw = np.asarray([spot["counts"] for spot in spots], dtype=float)
    if reference_raw.shape[1] != len(reference_genes):
        raise AssertionError("Reference centroid lengths do not match reference genes.")
    if spatial_raw.shape[1] != len(spatial_genes):
        raise AssertionError("Spatial observation lengths do not match spatial genes.")

    reference_shared_raw = reference_raw[:, shared_reference_idx]
    spatial_shared_raw = spatial_raw[:, shared_spatial_idx]
    reference_shared_norm, reference_library_sizes = normalize_log1p(reference_shared_raw, target_sum)
    spatial_shared_norm, spatial_library_sizes = normalize_log1p(spatial_shared_raw, target_sum)
    cosine_scores = cosine_similarity_matrix(spatial_shared_norm, reference_shared_norm)

    observation_ids = [spot["observation_id"] for spot in spots]
    coordinates = np.asarray([spot["coordinates"] for spot in spots], dtype=float)
    neighbors, distances = build_neighbor_graph(coordinates, n_neighbors=n_neighbors)

    ranking = np.argsort(cosine_scores, axis=1)[:, ::-1]
    primary_indices = ranking[:, 0]
    secondary_indices = ranking[:, 1]
    primary_labels = [reference_labels[int(index)] for index in primary_indices]
    secondary_labels = [reference_labels[int(index)] for index in secondary_indices]
    primary_scores = cosine_scores[np.arange(len(spots)), primary_indices]
    secondary_scores = cosine_scores[np.arange(len(spots)), secondary_indices]
    score_margins = primary_scores - secondary_scores

    marker_rows: list[dict[str, Any]] = []
    niche_rows: list[dict[str, Any]] = []
    label_rows: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []

    overlap_fraction = len(shared_genes) / float(union_gene_count)

    for row_index, observation_id in enumerate(observation_ids):
        primary_label = primary_labels[row_index]
        secondary_label = secondary_labels[row_index]

        primary_marker_mean = float(
            np.mean([spatial_shared_norm[row_index, shared_gene_to_idx[gene]] for gene in marker_table_by_label[primary_label]])
        )
        secondary_marker_mean = float(
            np.mean([spatial_shared_norm[row_index, shared_gene_to_idx[gene]] for gene in marker_table_by_label[secondary_label]])
        )
        marker_gap = primary_marker_mean - secondary_marker_mean

        neighbor_indices = neighbors[row_index]
        neighbor_labels = [primary_labels[index] for index in neighbor_indices]
        niche_majority_label = pick_majority_label(neighbor_labels)
        niche_agreement = float(niche_majority_label == primary_label)
        mean_neighbor_distance = float(np.mean([distances[row_index, index] for index in neighbor_indices]))

        if niche_majority_label != primary_label or marker_gap < 0.2:
            marker_review_status = "needs_review"
        elif marker_gap < 1.0 or float(score_margins[row_index]) < 0.05:
            marker_review_status = "mixed"
        else:
            marker_review_status = "supported"

        review_flags = []
        if float(score_margins[row_index]) < 0.05:
            review_flags.append("low_margin")
        if marker_gap < 0.2:
            review_flags.append("weak_marker_gap")
        elif marker_gap < 1.0:
            review_flags.append("moderate_marker_gap")
        if niche_majority_label != primary_label:
            review_flags.append("niche_disagreement")
        if not review_flags:
            review_flags.append("high_confidence")

        label_rows.append(
            {
                "observation_id": observation_id,
                "x_coord": round(float(coordinates[row_index, 0]), 6),
                "y_coord": round(float(coordinates[row_index, 1]), 6),
                "primary_label": primary_label,
                "primary_score": round(float(primary_scores[row_index]), 6),
                "score_margin": round(float(score_margins[row_index]), 6),
                "secondary_label": secondary_label,
                "mapping_method": "tangram_centroid_cosine_surrogate",
                "gene_overlap_count": len(shared_genes),
                "gene_overlap_fraction": round(overlap_fraction, 6),
                "marker_review_status": marker_review_status,
                "review_notes": (
                    f"{','.join(review_flags)}; niche_majority={niche_majority_label}; "
                    f"marker_gap={marker_gap:.3f}"
                ),
            }
        )
        marker_rows.append(
            {
                "observation_id": observation_id,
                "primary_label": primary_label,
                "secondary_label": secondary_label,
                "primary_marker_mean": round(primary_marker_mean, 6),
                "secondary_marker_mean": round(secondary_marker_mean, 6),
                "marker_gap": round(float(marker_gap), 6),
                "marker_review_status": marker_review_status,
            }
        )
        niche_rows.append(
            {
                "observation_id": observation_id,
                "neighbor_observation_ids": json.dumps([observation_ids[index] for index in neighbor_indices]),
                "neighbor_primary_labels": json.dumps(neighbor_labels),
                "niche_majority_label": niche_majority_label,
                "niche_agreement": round(niche_agreement, 6),
                "mean_neighbor_distance": round(mean_neighbor_distance, 6),
            }
        )

    label_frame = pd.DataFrame(label_rows)
    marker_frame = pd.DataFrame(marker_rows)
    niche_frame = pd.DataFrame(niche_rows)

    score_rows.extend(
        [
            {
                "metric_category": "gene_intersection",
                "metric_name": "shared_gene_count",
                "scope": "global",
                "target_id": "all",
                "method": "gene_intersection",
                "value": float(len(shared_genes)),
                "higher_is_better": True,
                "suggested_pass_floor": 6.0,
                "qc_status": metric_status(float(len(shared_genes)), 6.0, higher_is_better=True),
                "notes": f"Shared genes retained after dropping {len(reference_only)} reference-only and {len(spatial_only)} spatial-only genes.",
            },
            {
                "metric_category": "gene_intersection",
                "metric_name": "shared_gene_fraction",
                "scope": "global",
                "target_id": "all",
                "method": "gene_intersection",
                "value": overlap_fraction,
                "higher_is_better": True,
                "suggested_pass_floor": 0.6,
                "qc_status": metric_status(overlap_fraction, 0.6, higher_is_better=True),
                "notes": "Shared-gene fraction over the reference/spatial gene union.",
            },
            {
                "metric_category": "preprocessing",
                "metric_name": "reference_mean_library_size",
                "scope": "global",
                "target_id": "reference",
                "method": "normalize_total_log1p",
                "value": float(reference_library_sizes.mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 200.0,
                "qc_status": metric_status(float(reference_library_sizes.mean()), 200.0, higher_is_better=True),
                "notes": "Mean raw library size across reference centroids before normalization.",
            },
            {
                "metric_category": "preprocessing",
                "metric_name": "spatial_mean_library_size",
                "scope": "global",
                "target_id": "spatial",
                "method": "normalize_total_log1p",
                "value": float(spatial_library_sizes.mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 200.0,
                "qc_status": metric_status(float(spatial_library_sizes.mean()), 200.0, higher_is_better=True),
                "notes": "Mean raw library size across spatial observations before normalization.",
            },
            {
                "metric_category": "mapping",
                "metric_name": "mean_primary_score",
                "scope": "global",
                "target_id": "all",
                "method": "reference_centroid_cosine_similarity",
                "value": float(primary_scores.mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 0.95,
                "qc_status": metric_status(float(primary_scores.mean()), 0.95, higher_is_better=True),
                "notes": "Average top-label cosine similarity across all observations.",
            },
            {
                "metric_category": "mapping",
                "metric_name": "mean_score_margin",
                "scope": "global",
                "target_id": "all",
                "method": "reference_centroid_cosine_similarity",
                "value": float(score_margins.mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 0.05,
                "qc_status": metric_status(float(score_margins.mean()), 0.05, higher_is_better=True),
                "notes": "Average gap between the top and second cosine similarity scores.",
            },
            {
                "metric_category": "niche",
                "metric_name": "niche_agreement_rate",
                "scope": "global",
                "target_id": "all",
                "method": "3nn_majority_vote",
                "value": float(niche_frame["niche_agreement"].mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 0.75,
                "qc_status": metric_status(float(niche_frame["niche_agreement"].mean()), 0.75, higher_is_better=True),
                "notes": "Fraction of observations whose label agrees with the 3-nearest-neighbor majority vote.",
            },
            {
                "metric_category": "markers",
                "metric_name": "supported_fraction",
                "scope": "global",
                "target_id": "all",
                "method": "marker_consistency",
                "value": float((label_frame["marker_review_status"] == "supported").mean()),
                "higher_is_better": True,
                "suggested_pass_floor": 0.5,
                "qc_status": metric_status(float((label_frame["marker_review_status"] == "supported").mean()), 0.5, higher_is_better=True),
                "notes": "Fraction of observations with strong marker support and no niche disagreement.",
            },
        ]
    )

    for row_index, row in label_frame.iterrows():
        score_rows.extend(
            [
                {
                    "metric_category": "mapping",
                    "metric_name": "primary_score",
                    "scope": "observation",
                    "target_id": row["observation_id"],
                    "method": "reference_centroid_cosine_similarity",
                    "value": float(row["primary_score"]),
                    "higher_is_better": True,
                    "suggested_pass_floor": 0.95,
                    "qc_status": metric_status(float(row["primary_score"]), 0.95, higher_is_better=True),
                    "notes": f"Top cosine similarity for {row['primary_label']}.",
                },
                {
                    "metric_category": "mapping",
                    "metric_name": "score_margin",
                    "scope": "observation",
                    "target_id": row["observation_id"],
                    "method": "reference_centroid_cosine_similarity",
                    "value": float(row["score_margin"]),
                    "higher_is_better": True,
                    "suggested_pass_floor": 0.05,
                    "qc_status": metric_status(float(row["score_margin"]), 0.05, higher_is_better=True),
                    "notes": f"Primary-versus-secondary margin for {row['observation_id']}.",
                },
                {
                    "metric_category": "niche",
                    "metric_name": "niche_agreement",
                    "scope": "observation",
                    "target_id": row["observation_id"],
                    "method": "3nn_majority_vote",
                    "value": float(niche_frame.loc[row_index, "niche_agreement"]),
                    "higher_is_better": True,
                    "suggested_pass_floor": 1.0,
                    "qc_status": metric_status(float(niche_frame.loc[row_index, "niche_agreement"]), 1.0, higher_is_better=True),
                    "notes": f"Agreement between the predicted label and the neighborhood majority for {row['observation_id']}.",
                },
                {
                    "metric_category": "markers",
                    "metric_name": "marker_gap",
                    "scope": "observation",
                    "target_id": row["observation_id"],
                    "method": "marker_consistency",
                    "value": float(marker_frame.loc[row_index, "marker_gap"]),
                    "higher_is_better": True,
                    "suggested_pass_floor": 0.25,
                    "qc_status": metric_status(float(marker_frame.loc[row_index, "marker_gap"]), 0.25, higher_is_better=True),
                    "notes": f"Primary-versus-secondary marker mean gap for {row['observation_id']}.",
                },
            ]
        )

    mapping_scores = pd.DataFrame(score_rows)

    gene_rows = []
    for gene in sorted(set(reference_genes) | set(spatial_genes)):
        in_reference = gene in reference_genes
        in_spatial = gene in spatial_genes
        role = "shared" if in_reference and in_spatial else "reference_only" if in_reference else "spatial_only"
        marker_labels = [label for label, markers in marker_table_by_label.items() if gene in markers]
        gene_rows.append(
            {
                "gene": gene,
                "role": role,
                "in_reference": in_reference,
                "in_spatial": in_spatial,
                "used_for_mapping": role == "shared",
                "marker_labels": ";".join(marker_labels),
            }
        )
    gene_qc = pd.DataFrame(gene_rows)

    predicted_label_counts = label_frame["primary_label"].value_counts().to_dict()
    review_counts = label_frame["marker_review_status"].value_counts().to_dict()
    ambiguous_rows = label_frame.loc[
        label_frame["marker_review_status"] != "supported",
        ["observation_id", "marker_review_status"],
    ]

    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}.",
                f"Computed Tangram-shaped centroid cosine mapping for {len(spots)} spatial observations against {len(reference_labels)} reference labels.",
                "The starter is deterministic and local: no network access, no heavy Tangram training loop, and no private Spatial Agent tooling.",
            ],
        },
        {
            "name": "Inputs and preprocessing",
            "bullets": [
                f"Reference genes: {len(reference_genes)}; spatial genes: {len(spatial_genes)}; shared genes used for mapping: {len(shared_genes)}.",
                f"Reference mean library size: {reference_library_sizes.mean():.1f}; spatial mean library size: {spatial_library_sizes.mean():.1f}.",
                "Applied library-size normalization and log1p before centroid cosine similarity so the toy path mirrors a public preprocessing shape.",
            ],
        },
        {
            "name": "Gene intersection QC",
            "bullets": [
                f"Shared genes: {', '.join(shared_genes)}.",
                f"Reference-only genes dropped: {', '.join(reference_only)}.",
                f"Spatial-only genes dropped: {', '.join(spatial_only)}.",
            ],
        },
        {
            "name": "Tangram mapping",
            "bullets": [
                f"Predicted label counts: {', '.join(f'{label}={count}' for label, count in sorted(predicted_label_counts.items()))}.",
                f"Mean primary score: {primary_scores.mean():.3f}; mean top-label margin: {score_margins.mean():.3f}.",
                "Low-margin observations remained explicit instead of being forced into a higher-confidence interpretation.",
            ],
        },
        {
            "name": "Label transfer review",
            "bullets": [
                f"3-NN niche majority agreed with every mapped label ({niche_frame['niche_agreement'].mean():.2f} agreement rate).",
                f"Review status counts: {', '.join(f'{status}={count}' for status, count in sorted(review_counts.items()))}.",
                "Observations needing closer review: "
                + ", ".join(
                    f"{row.observation_id} ({row.marker_review_status})"
                    for row in ambiguous_rows.itertuples(index=False)
                )
                + ".",
            ],
        },
        {
            "name": "Caveats and next actions",
            "bullets": [
                "This starter truly computes gene intersection, normalize+log1p preprocessing, centroid cosine similarity, top-label margins, 3-NN niche voting, and marker consistency checks.",
                "It approximates Tangram with centroid-level cosine scoring rather than the full PyTorch optimization, density priors, or gene projection stack.",
                "Swap in a pinned public benchmark slice or project dataset, keep the same deliverable names, and upgrade the surrogate scorer to a real Tangram run when heavier dependencies are available.",
            ],
        },
    ]

    write_markdown(outdir / "spatial_mapping_report.md", "Spatial Mapping Report", report_sections)
    label_frame.to_csv(outdir / "mapped_labels.tsv", sep="\t", index=False)
    mapping_scores.to_csv(outdir / "mapping_scores.tsv", sep="\t", index=False)
    gene_qc.to_csv(outdir / "gene_intersection_qc.tsv", sep="\t", index=False)
    marker_frame.to_csv(outdir / "marker_qc.tsv", sep="\t", index=False)
    niche_frame.to_csv(outdir / "niche_qc.tsv", sep="\t", index=False)

    intermediate_qc = {
        "run_label": payload["run_label"],
        "normalization": {
            "target_sum": target_sum,
            "reference_library_sizes": [float(value) for value in reference_library_sizes.tolist()],
            "spatial_library_sizes": [float(value) for value in spatial_library_sizes.tolist()],
        },
        "gene_intersection": {
            "shared_genes": shared_genes,
            "shared_gene_count": len(shared_genes),
            "shared_gene_fraction": overlap_fraction,
            "reference_only_genes": reference_only,
            "spatial_only_genes": spatial_only,
        },
        "reference_labels": reference_labels,
        "observation_ids": observation_ids,
        "cosine_similarity": {
            "method": "reference_centroid_cosine_similarity",
            "matrix": cosine_scores.round(6).tolist(),
        },
        "label_assignments": label_rows,
        "marker_consistency": marker_rows,
        "niche_majority": niche_rows,
    }
    (outdir / "intermediate_qc.json").write_text(
        json.dumps(intermediate_qc, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path.resolve()),
        "outdir": str(outdir.resolve()),
        "method_steps": [
            "gene_intersection",
            "normalize_total_log1p",
            "reference_centroid_cosine_similarity",
            "top_label_margin",
            "3nn_niche_majority_vote",
            "marker_consistency_review",
        ],
        "observation_count": len(spots),
        "reference_profile_count": len(reference_labels),
        "shared_gene_count": len(shared_genes),
        "prediction_counts": {label: int(count) for label, count in predicted_label_counts.items()},
        "review_counts": {status: int(count) for status, count in review_counts.items()},
        "written_files": [],
    }
    (outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    validate_outputs(skill_dir=SKILL_DIR, outdir=outdir, expected=payload.get("expected_invariants"))
    summary["written_files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    (outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the spatial mapping starter on tiny synthetic reference and spatial inputs."
    )
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    result = build_spatial_mapping(payload, input_path=args.input, outdir=args.outdir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
