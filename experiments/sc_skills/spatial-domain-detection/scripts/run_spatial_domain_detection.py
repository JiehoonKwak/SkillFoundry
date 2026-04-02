#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
QC_FILES = [
    "qc_spatial_graph.tsv",
    "qc_domain_scores.tsv",
    "qc_method_comparison.tsv",
    "run_summary.json",
]
METHODS = ("spagcn_like", "graphst_like")


@dataclass
class MethodResult:
    method: str
    assigned_indices: np.ndarray
    assigned_domains: list[str]
    probabilities: np.ndarray
    confidence_margin: np.ndarray
    graph: np.ndarray
    graph_rows: list[dict[str, Any]]
    qc_rows: list[dict[str, Any]]
    metrics: dict[str, float]
    marker_table: pd.DataFrame


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


def sanitize_label(label: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in label.lower()).strip("_")


def softmax_rows(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def probability_margin(probabilities: np.ndarray) -> np.ndarray:
    top_two = np.partition(probabilities, -2, axis=1)[:, -2:]
    return top_two[:, 1] - top_two[:, 0]


def count_matrix_from_input(payload: dict[str, Any]) -> tuple[list[str], list[str], np.ndarray, np.ndarray, np.ndarray]:
    genes = payload["genes"]
    spots = payload["spots"]
    cell_ids = [spot["cell_id"] for spot in spots]
    coordinates = np.asarray([[spot["x"], spot["y"]] for spot in spots], dtype=float)
    histology_rgb = np.asarray([spot["histology_rgb"] for spot in spots], dtype=float)
    counts = np.asarray([[spot["counts"][gene] for gene in genes] for spot in spots], dtype=float)

    if len(set(cell_ids)) != len(cell_ids):
        raise AssertionError("Cell identifiers must be unique.")
    if np.any(counts <= 0):
        raise AssertionError("Toy counts must be strictly positive for deterministic log-normalization.")
    if counts.shape[0] < 2:
        raise AssertionError("At least two spots are required.")

    return genes, cell_ids, counts, coordinates, histology_rgb


def normalize_counts(counts: np.ndarray, scale_factor: float) -> tuple[np.ndarray, np.ndarray]:
    library_size = counts.sum(axis=1, keepdims=True)
    if np.any(library_size <= 0):
        raise AssertionError("Every spot must have a positive library size.")
    normalized = np.log1p((counts / library_size) * float(scale_factor))
    return normalized, library_size[:, 0]


def build_signature_matrix(
    domain_names: list[str],
    genes: list[str],
    domain_signatures: dict[str, dict[str, float]],
) -> np.ndarray:
    matrix = np.asarray(
        [
            [float(domain_signatures[domain].get(gene, 0.0)) for gene in genes]
            for domain in domain_names
        ],
        dtype=float,
    )
    if np.any(np.all(np.isclose(matrix, 0.0), axis=1)):
        raise AssertionError("Each domain signature needs at least one non-zero gene weight.")
    return matrix


def build_spatial_graph(
    *,
    cell_ids: list[str],
    coordinates: np.ndarray,
    histology_rgb: np.ndarray | None,
    neighbor_count: int,
    self_weight: float,
    histology_alpha: float,
    method: str,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    if not 0.0 <= self_weight < 1.0:
        raise AssertionError("self_weight must be in [0, 1).")

    graph = np.zeros((len(cell_ids), len(cell_ids)), dtype=float)
    rows: list[dict[str, Any]] = []
    capped_neighbors = max(1, min(int(neighbor_count), len(cell_ids) - 1))

    for source_index, source_cell in enumerate(cell_ids):
        distances = np.linalg.norm(coordinates[source_index] - coordinates, axis=1)
        order = [index for index in np.argsort(distances) if index != source_index][:capped_neighbors]
        raw_weights: list[float] = []

        for target_index in order:
            spatial_term = 1.0 / max(float(distances[target_index]), 1e-9)
            histology_distance = (
                float(np.linalg.norm(histology_rgb[source_index] - histology_rgb[target_index]))
                if histology_rgb is not None
                else 0.0
            )
            histology_affinity = math.exp(-float(histology_alpha) * histology_distance) if histology_rgb is not None else 1.0
            raw_weights.append(spatial_term * histology_affinity)

        normalized_weights = np.asarray(raw_weights, dtype=float)
        normalized_weights = normalized_weights / normalized_weights.sum()

        graph[source_index, source_index] = float(self_weight)
        rows.append(
            {
                "method": method,
                "source_cell": source_cell,
                "target_cell": source_cell,
                "relation": "self",
                "distance": 0.0,
                "histology_distance": 0.0,
                "weight": round(float(self_weight), 6),
            }
        )

        for target_index, normalized_weight in zip(order, normalized_weights, strict=True):
            distance = float(distances[target_index])
            histology_distance = (
                float(np.linalg.norm(histology_rgb[source_index] - histology_rgb[target_index]))
                if histology_rgb is not None
                else 0.0
            )
            edge_weight = float((1.0 - self_weight) * normalized_weight)
            graph[source_index, target_index] = edge_weight
            rows.append(
                {
                    "method": method,
                    "source_cell": source_cell,
                    "target_cell": cell_ids[target_index],
                    "relation": "neighbor",
                    "distance": round(distance, 6),
                    "histology_distance": round(histology_distance, 6),
                    "weight": round(edge_weight, 6),
                }
            )

    if not np.allclose(graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError(f"{method} graph rows must sum to 1.")

    return graph, rows


def farthest_point_initialization(embedding: np.ndarray, cluster_count: int, coordinates: np.ndarray) -> np.ndarray:
    selected = [int(np.argmin(coordinates[:, 0]))]
    while len(selected) < cluster_count:
        distances = np.min(
            ((embedding[:, None, :] - embedding[np.asarray(selected)][None, :, :]) ** 2).sum(axis=2),
            axis=1,
        )
        distances[selected] = -1.0
        selected.append(int(np.argmax(distances)))
    return np.asarray(selected, dtype=int)


def kmeans(embedding: np.ndarray, cluster_count: int, coordinates: np.ndarray, max_iterations: int = 20) -> tuple[np.ndarray, np.ndarray]:
    centroids = embedding[farthest_point_initialization(embedding, cluster_count, coordinates)]
    labels = np.zeros(len(embedding), dtype=int)

    for _ in range(max_iterations):
        squared_distance = ((embedding[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
        labels = np.argmin(squared_distance, axis=1)
        updated = []
        for cluster_id in range(cluster_count):
            members = embedding[labels == cluster_id]
            if len(members) == 0:
                raise AssertionError("Deterministic k-means produced an empty cluster on the toy input.")
            updated.append(members.mean(axis=0))
        updated_centroids = np.vstack(updated)
        if np.allclose(updated_centroids, centroids):
            break
        centroids = updated_centroids

    return labels, centroids


def map_clusters_to_domains(cluster_labels: np.ndarray, base_scores: np.ndarray, domain_count: int) -> tuple[np.ndarray, dict[int, int]]:
    cluster_domain_scores = np.zeros((domain_count, domain_count), dtype=float)
    for cluster_id in range(domain_count):
        members = base_scores[cluster_labels == cluster_id]
        if len(members) == 0:
            raise AssertionError("Every GraphST-like cluster must contain at least one spot.")
        cluster_domain_scores[cluster_id] = members.mean(axis=0)

    _, best_mapping = max(
        (
            (
                float(sum(cluster_domain_scores[cluster_id, permutation[cluster_id]] for cluster_id in range(domain_count))),
                permutation,
            )
            for permutation in itertools.permutations(range(domain_count))
        ),
        key=lambda item: item[0],
    )
    return cluster_domain_scores, {cluster_id: domain_index for cluster_id, domain_index in enumerate(best_mapping)}


def within_label_histology_score(histology_rgb: np.ndarray, labels: np.ndarray, domain_count: int) -> float:
    within_domain_distances = []
    for domain_index in range(domain_count):
        members = histology_rgb[labels == domain_index]
        if len(members) == 0:
            raise AssertionError("Every domain needs at least one member for histology consistency.")
        centroid = members.mean(axis=0)
        within_domain_distances.append(float(np.mean(np.linalg.norm(members - centroid, axis=1))))
    return float(math.exp(-float(np.mean(within_domain_distances))))


def graph_edge_agreement(graph: np.ndarray, labels: np.ndarray) -> tuple[float, np.ndarray]:
    supports = []
    numerator = 0.0
    denominator = 0.0
    for source_index in range(len(labels)):
        same_label_support = 0.0
        for target_index in range(len(labels)):
            weight = float(graph[source_index, target_index])
            if source_index == target_index or weight == 0.0:
                continue
            denominator += weight
            if labels[source_index] == labels[target_index]:
                numerator += weight
                same_label_support += weight
        supports.append(same_label_support)
    return numerator / max(denominator, 1e-9), np.asarray(supports, dtype=float)


def compute_marker_table(
    *,
    expression: np.ndarray,
    genes: list[str],
    domain_names: list[str],
    assigned_indices: np.ndarray,
    top_n: int,
    method: str,
) -> tuple[pd.DataFrame, float]:
    rows: list[dict[str, Any]] = []
    strongest_effects: list[float] = []

    for domain_index, domain_name in enumerate(domain_names):
        in_domain = assigned_indices == domain_index
        out_domain = assigned_indices != domain_index
        if in_domain.sum() == 0 or out_domain.sum() == 0:
            raise AssertionError(f"Domain {domain_name} must have both in-domain and out-domain spots.")

        mean_in = expression[in_domain].mean(axis=0)
        mean_out = expression[out_domain].mean(axis=0)
        effect_size = mean_in - mean_out
        ordered_genes = np.argsort(effect_size)[::-1]

        added = 0
        for gene_index in ordered_genes:
            if effect_size[gene_index] <= 0:
                continue
            rows.append(
                {
                    "domain_label": domain_name,
                    "marker_gene": genes[gene_index],
                    "effect_size": round(float(effect_size[gene_index]), 6),
                    "mean_in_domain": round(float(mean_in[gene_index]), 6),
                    "mean_out_domain": round(float(mean_out[gene_index]), 6),
                    "method": method,
                    "rank": added + 1,
                }
            )
            added += 1
            if added >= top_n:
                break

        if added == 0:
            raise AssertionError(f"No positive marker effect was found for domain {domain_name}.")
        strongest_effects.append(float(effect_size[ordered_genes[0]]))

    marker_table = pd.DataFrame(rows)
    marker_table = marker_table.sort_values(["domain_label", "rank", "marker_gene"]).reset_index(drop=True)
    return marker_table, float(np.mean(strongest_effects))


def build_method_metrics(
    *,
    graph: np.ndarray,
    histology_rgb: np.ndarray,
    assigned_indices: np.ndarray,
    probabilities: np.ndarray,
    marker_strength_mean: float,
    shifted_from_base: np.ndarray,
    domain_count: int,
) -> tuple[dict[str, float], np.ndarray]:
    edge_agreement, local_support = graph_edge_agreement(graph, assigned_indices)
    histology_consistency = within_label_histology_score(histology_rgb, assigned_indices, domain_count)
    mean_confidence_margin = float(np.mean(probability_margin(probabilities)))
    marker_strength_normalized = float(1.0 - math.exp(-marker_strength_mean))
    composite_score = (
        0.4 * edge_agreement
        + 0.35 * histology_consistency
        + 0.25 * marker_strength_normalized
    )
    return (
        {
            "weighted_edge_agreement": float(edge_agreement),
            "histology_consistency": float(histology_consistency),
            "mean_confidence_margin": float(mean_confidence_margin),
            "mean_top_marker_effect": float(marker_strength_mean),
            "graph_shift_fraction": float(np.mean(shifted_from_base.astype(float))),
            "composite_score": float(composite_score),
        },
        local_support,
    )


def run_spagcn_like(
    *,
    cell_ids: list[str],
    genes: list[str],
    coordinates: np.ndarray,
    histology_rgb: np.ndarray,
    expression: np.ndarray,
    base_scores: np.ndarray,
    base_best_indices: np.ndarray,
    domain_names: list[str],
    neighbor_count: int,
    params: dict[str, float],
    marker_top_n: int,
) -> MethodResult:
    graph, graph_rows = build_spatial_graph(
        cell_ids=cell_ids,
        coordinates=coordinates,
        histology_rgb=histology_rgb,
        neighbor_count=neighbor_count,
        self_weight=float(params["self_weight"]),
        histology_alpha=float(params["histology_alpha"]),
        method="spagcn_like",
    )

    smoothed_scores = base_scores.copy()
    for _ in range(int(params["smoothing_rounds"])):
        smoothed_scores = (
            (1.0 - float(params["smoothing_weight"])) * smoothed_scores
            + float(params["smoothing_weight"]) * (graph @ smoothed_scores)
        )

    probabilities = softmax_rows(smoothed_scores)
    assigned_indices = probabilities.argmax(axis=1)
    assigned_domains = [domain_names[index] for index in assigned_indices]
    confidence_margin = probability_margin(probabilities)
    shifted_from_base = assigned_indices != base_best_indices
    marker_table, marker_strength_mean = compute_marker_table(
        expression=expression,
        genes=genes,
        domain_names=domain_names,
        assigned_indices=assigned_indices,
        top_n=marker_top_n,
        method="spagcn_like",
    )
    metrics, local_support = build_method_metrics(
        graph=graph,
        histology_rgb=histology_rgb,
        assigned_indices=assigned_indices,
        probabilities=probabilities,
        marker_strength_mean=marker_strength_mean,
        shifted_from_base=shifted_from_base,
        domain_count=len(domain_names),
    )

    qc_rows = []
    for row_index, cell_id in enumerate(cell_ids):
        row: dict[str, Any] = {
            "cell_id": cell_id,
            "method": "spagcn_like",
            "base_best_domain": domain_names[base_best_indices[row_index]],
            "assigned_domain": assigned_domains[row_index],
            "domain_score": round(float(probabilities[row_index].max()), 6),
            "confidence_margin": round(float(confidence_margin[row_index]), 6),
            "graph_shifted": bool(shifted_from_base[row_index]),
            "neighbor_same_domain_weight": round(float(local_support[row_index]), 6),
            "cluster_id": -1,
            "embedding_1": np.nan,
            "embedding_2": np.nan,
        }
        for domain_index, domain_name in enumerate(domain_names):
            slug = sanitize_label(domain_name)
            row[f"base_{slug}_score"] = round(float(base_scores[row_index, domain_index]), 6)
            row[f"final_{slug}_score"] = round(float(smoothed_scores[row_index, domain_index]), 6)
            row[f"{slug}_probability"] = round(float(probabilities[row_index, domain_index]), 6)
        qc_rows.append(row)

    return MethodResult(
        method="spagcn_like",
        assigned_indices=assigned_indices,
        assigned_domains=assigned_domains,
        probabilities=probabilities,
        confidence_margin=confidence_margin,
        graph=graph,
        graph_rows=graph_rows,
        qc_rows=qc_rows,
        metrics=metrics,
        marker_table=marker_table,
    )


def run_graphst_like(
    *,
    cell_ids: list[str],
    genes: list[str],
    coordinates: np.ndarray,
    histology_rgb: np.ndarray,
    expression: np.ndarray,
    base_scores: np.ndarray,
    base_best_indices: np.ndarray,
    domain_names: list[str],
    neighbor_count: int,
    params: dict[str, float],
    marker_top_n: int,
) -> MethodResult:
    cluster_count = int(params["cluster_count"])
    if cluster_count != len(domain_names):
        raise AssertionError("The GraphST-like starter expects cluster_count to match the candidate domain count.")

    graph, graph_rows = build_spatial_graph(
        cell_ids=cell_ids,
        coordinates=coordinates,
        histology_rgb=None,
        neighbor_count=neighbor_count,
        self_weight=float(params["self_weight"]),
        histology_alpha=0.0,
        method="graphst_like",
    )

    feature_std = expression.std(axis=0)
    feature_std[feature_std == 0] = 1.0
    scaled_expression = (expression - expression.mean(axis=0, keepdims=True)) / feature_std

    coordinate_range = np.ptp(coordinates, axis=0)
    coordinate_range[coordinate_range == 0] = 1.0
    scaled_coordinates = (coordinates - coordinates.min(axis=0, keepdims=True)) / coordinate_range

    latent = np.column_stack([scaled_expression, scaled_coordinates])
    latent = (
        (1.0 - float(params["latent_mix"])) * latent
        + float(params["latent_mix"]) * (graph @ latent)
    )
    latent = latent - latent.mean(axis=0, keepdims=True)

    embedding_dims = min(int(params["embedding_dims"]), latent.shape[0], latent.shape[1])
    left_singular_vectors, singular_values, _ = np.linalg.svd(latent, full_matrices=False)
    embedding = left_singular_vectors[:, :embedding_dims] * singular_values[:embedding_dims]
    if embedding.shape[1] == 1:
        embedding = np.column_stack([embedding[:, 0], np.zeros(len(embedding), dtype=float)])

    cluster_labels, centroids = kmeans(embedding, cluster_count, coordinates)
    _, cluster_to_domain = map_clusters_to_domains(cluster_labels, base_scores, len(domain_names))

    squared_distance = ((embedding[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    cluster_probabilities = softmax_rows(-squared_distance)
    domain_probabilities = np.zeros((len(cell_ids), len(domain_names)), dtype=float)
    for cluster_id, domain_index in cluster_to_domain.items():
        domain_probabilities[:, domain_index] = cluster_probabilities[:, cluster_id]

    assigned_indices = domain_probabilities.argmax(axis=1)
    assigned_domains = [domain_names[index] for index in assigned_indices]
    confidence_margin = probability_margin(domain_probabilities)
    shifted_from_base = assigned_indices != base_best_indices

    marker_table, marker_strength_mean = compute_marker_table(
        expression=expression,
        genes=genes,
        domain_names=domain_names,
        assigned_indices=assigned_indices,
        top_n=marker_top_n,
        method="graphst_like",
    )
    metrics, local_support = build_method_metrics(
        graph=graph,
        histology_rgb=histology_rgb,
        assigned_indices=assigned_indices,
        probabilities=domain_probabilities,
        marker_strength_mean=marker_strength_mean,
        shifted_from_base=shifted_from_base,
        domain_count=len(domain_names),
    )

    qc_rows = []
    for row_index, cell_id in enumerate(cell_ids):
        row = {
            "cell_id": cell_id,
            "method": "graphst_like",
            "base_best_domain": domain_names[base_best_indices[row_index]],
            "assigned_domain": assigned_domains[row_index],
            "domain_score": round(float(domain_probabilities[row_index].max()), 6),
            "confidence_margin": round(float(confidence_margin[row_index]), 6),
            "graph_shifted": bool(shifted_from_base[row_index]),
            "neighbor_same_domain_weight": round(float(local_support[row_index]), 6),
            "cluster_id": int(cluster_labels[row_index]),
            "embedding_1": round(float(embedding[row_index, 0]), 6),
            "embedding_2": round(float(embedding[row_index, 1]), 6),
        }
        for domain_index, domain_name in enumerate(domain_names):
            slug = sanitize_label(domain_name)
            row[f"base_{slug}_score"] = round(float(base_scores[row_index, domain_index]), 6)
            row[f"final_{slug}_score"] = round(float(domain_probabilities[row_index, domain_index]), 6)
            row[f"{slug}_probability"] = round(float(domain_probabilities[row_index, domain_index]), 6)
        qc_rows.append(row)

    return MethodResult(
        method="graphst_like",
        assigned_indices=assigned_indices,
        assigned_domains=assigned_domains,
        probabilities=domain_probabilities,
        confidence_margin=confidence_margin,
        graph=graph,
        graph_rows=graph_rows,
        qc_rows=qc_rows,
        metrics=metrics,
        marker_table=marker_table,
    )


def choose_selected_method(
    requested_method: str,
    method_results: dict[str, MethodResult],
) -> str:
    if requested_method != "best":
        return requested_method

    ranked = sorted(
        method_results.values(),
        key=lambda result: (
            result.metrics["composite_score"],
            result.metrics["weighted_edge_agreement"],
            result.metrics["histology_consistency"],
        ),
        reverse=True,
    )
    return ranked[0].method


def write_outputs(
    *,
    outdir: Path,
    payload: dict[str, Any],
    cell_ids: list[str],
    genes: list[str],
    library_size: np.ndarray,
    base_best_indices: np.ndarray,
    domain_names: list[str],
    method_results: dict[str, MethodResult],
    selected_method: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    selected_result = method_results[selected_method]
    selected_labels = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "domain_label": selected_result.assigned_domains,
            "domain_score": np.round(selected_result.probabilities.max(axis=1), 6),
            "selected_method": selected_method,
            "confidence_margin": np.round(selected_result.confidence_margin, 6),
            "base_best_domain": [domain_names[index] for index in base_best_indices],
            "spagcn_like_label": method_results["spagcn_like"].assigned_domains,
            "graphst_like_label": method_results["graphst_like"].assigned_domains,
            "library_size": np.round(library_size, 6),
        }
    )
    selected_labels.to_csv(outdir / "domain_labels.tsv", sep="\t", index=False)

    selected_markers = method_results[selected_method].marker_table.drop(columns=["method", "rank"])
    selected_markers.to_csv(outdir / "domain_markers.tsv", sep="\t", index=False)

    graph_rows = pd.DataFrame(
        [row for result in method_results.values() for row in result.graph_rows]
    )
    graph_rows.to_csv(outdir / "qc_spatial_graph.tsv", sep="\t", index=False)

    domain_score_rows = pd.DataFrame(
        [row for result in method_results.values() for row in result.qc_rows]
    )
    domain_score_rows.to_csv(outdir / "qc_domain_scores.tsv", sep="\t", index=False)

    comparison_rows = []
    for method, result in method_results.items():
        for metric, value in result.metrics.items():
            comparison_rows.append({"method": method, "metric": metric, "value": round(float(value), 6)})
        comparison_rows.append(
            {
                "method": method,
                "metric": "selected_output",
                "value": 1.0 if method == selected_method else 0.0,
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(outdir / "qc_method_comparison.tsv", sep="\t", index=False)

    boundary_spots = selected_labels.loc[
        selected_labels["domain_label"] != selected_labels["base_best_domain"],
        "cell_id",
    ].tolist()
    domain_counts = selected_labels["domain_label"].value_counts().sort_index().to_dict()

    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}",
                f"Toy input size: {len(cell_ids)} spots, {len(genes)} genes, {len(domain_names)} candidate domains.",
                (
                    f"Selected partition: {selected_method} with composite score "
                    f"{method_results[selected_method].metrics['composite_score']:.3f} "
                    f"versus {method_results['graphst_like' if selected_method == 'spagcn_like' else 'spagcn_like'].method} "
                    f"{method_results['graphst_like' if selected_method == 'spagcn_like' else 'spagcn_like'].metrics['composite_score']:.3f}."
                ),
            ],
        },
        {
            "name": "Domain labels",
            "bullets": [
                "Domain counts: "
                + ", ".join(f"{domain}={count}" for domain, count in sorted(domain_counts.items())),
                (
                    "Boundary evidence: "
                    + ", ".join(
                        f"{cell_id} moved from {selected_labels.set_index('cell_id').loc[cell_id, 'base_best_domain']} "
                        f"to {selected_labels.set_index('cell_id').loc[cell_id, 'domain_label']}"
                        for cell_id in boundary_spots
                    )
                    if boundary_spots
                    else "No spots changed relative to the raw signature ranking."
                ),
                (
                    "Method disagreement remained localized to the bridge region: "
                    + ", ".join(
                        f"{cell_id}={method_results['spagcn_like'].assigned_domains[cell_ids.index(cell_id)]}/"
                        f"{method_results['graphst_like'].assigned_domains[cell_ids.index(cell_id)]}"
                        for cell_id in boundary_spots
                    )
                    if boundary_spots
                    else "SpaGCN-like and GraphST-like assignments matched for all spots."
                ),
            ],
        },
        {
            "name": "Markers",
            "bullets": [
                (
                    f"{domain}: "
                    + ", ".join(
                        f"{row.marker_gene} ({row.effect_size:.2f})"
                        for row in selected_markers.loc[selected_markers["domain_label"] == domain].itertuples()
                    )
                )
                for domain in domain_names
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "SpaGCN-like scoring uses a toy RGB histology proxy and graph-smoothed signature scores instead of a trained graph convolutional network.",
                "GraphST-like scoring uses deterministic graph diffusion, SVD, and k-means instead of self-supervised contrastive learning.",
                "Marker interpretation uses domain-vs-rest effect sizes on normalized toy counts rather than full differential-expression testing.",
            ],
        },
    ]
    write_markdown(outdir / "domain_detection_report.md", "Domain Detection Report", report_sections)

    summary = {
        "boundary_spots": boundary_spots,
        "domain_counts": domain_counts,
        "gene_count": len(genes),
        "method_steps": [
            "spatial_graph_construction",
            "signature_scoring",
            "spagcn_like_refinement",
            "graphst_like_partitioning",
        ],
        "run_label": payload["run_label"],
        "selected_method": selected_method,
        "spot_count": len(cell_ids),
    }
    for method, result in method_results.items():
        summary[f"{method}_composite_score"] = round(float(result.metrics["composite_score"]), 6)
    write_json(outdir / "run_summary.json", summary)

    validate_outputs(SKILL_DIR, outdir)


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
                raise AssertionError(f"Empty TSV deliverable: {path}")
        elif deliverable["kind"] == "md":
            validate_markdown_sections(path, deliverable.get("required_sections", []))
        else:
            raise AssertionError(f"Unsupported deliverable kind: {deliverable['kind']}")

    for filename in metadata.get("starter_qc_files", QC_FILES):
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    labels = pd.read_csv(outdir / "domain_labels.tsv", sep="\t")
    if not labels["domain_score"].between(0.0, 1.0).all():
        raise AssertionError("domain_score must stay in [0, 1].")
    if labels["cell_id"].duplicated().any():
        raise AssertionError("domain_labels.tsv must have unique cell_id values.")

    markers = pd.read_csv(outdir / "domain_markers.tsv", sep="\t")
    if (markers["effect_size"] <= 0).any():
        raise AssertionError("All toy marker effect sizes must be positive.")

    graph = pd.read_csv(outdir / "qc_spatial_graph.tsv", sep="\t")
    if set(graph["method"]) != set(METHODS):
        raise AssertionError("qc_spatial_graph.tsv must include both starter methods.")
    row_sums = graph.groupby(["method", "source_cell"])["weight"].sum().to_numpy(dtype=float)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise AssertionError("Each source row in the spatial graph QC must sum to 1.")

    score_frame = pd.read_csv(outdir / "qc_domain_scores.tsv", sep="\t")
    if set(score_frame["method"]) != set(METHODS):
        raise AssertionError("qc_domain_scores.tsv must include both starter methods.")
    counts_per_method = score_frame.groupby("method")["cell_id"].nunique().to_dict()
    if len(set(counts_per_method.values())) != 1:
        raise AssertionError("Each method must score the same number of cells.")
    if not score_frame["domain_score"].between(0.0, 1.0).all():
        raise AssertionError("QC domain scores must stay in [0, 1].")

    comparison = pd.read_csv(outdir / "qc_method_comparison.tsv", sep="\t")
    comparison_flags = comparison.loc[comparison["metric"] == "selected_output", "value"].to_numpy(dtype=float)
    if comparison_flags.sum() != 1.0:
        raise AssertionError("Exactly one method must be marked as selected_output.")

    summary = load_json(outdir / "run_summary.json")
    if summary["selected_method"] not in METHODS:
        raise AssertionError("run_summary.json must record one of the starter methods.")
    if sum(summary["domain_counts"].values()) != len(labels):
        raise AssertionError("run_summary.json domain_counts must sum to the number of labeled spots.")

    selected_rows = score_frame.loc[score_frame["method"] == summary["selected_method"]].set_index("cell_id")
    label_rows = labels.set_index("cell_id")
    if not selected_rows.index.equals(label_rows.index):
        raise AssertionError("Selected QC rows and label rows must share the same cell order.")
    if list(selected_rows["assigned_domain"]) != list(label_rows["domain_label"]):
        raise AssertionError("domain_labels.tsv must match the selected method assignments.")


def run_skill(input_path: Path, outdir: Path, selected_method: str) -> None:
    payload = load_json(input_path)
    genes, cell_ids, counts, coordinates, histology_rgb = count_matrix_from_input(payload)
    expression, library_size = normalize_counts(counts, float(payload["parameters"]["expression_scale"]))
    domain_names = payload["candidate_domains"]
    signature_matrix = build_signature_matrix(domain_names, genes, payload["domain_signatures"])
    base_scores = expression @ signature_matrix.T
    base_best_indices = base_scores.argmax(axis=1)
    neighbor_count = int(payload["parameters"]["neighbor_count"])
    marker_top_n = int(payload["parameters"]["marker_top_n"])
    expected = payload.get("expected_invariants", {})

    if expected.get("spot_count") is not None and int(expected["spot_count"]) != len(cell_ids):
        raise AssertionError("Toy input spot_count invariant does not match the raw spot table.")
    if expected.get("gene_count") is not None and int(expected["gene_count"]) != len(genes):
        raise AssertionError("Toy input gene_count invariant does not match the raw gene list.")
    if expected.get("domain_count") is not None and int(expected["domain_count"]) != len(domain_names):
        raise AssertionError("Toy input domain_count invariant does not match the candidate domains.")

    spagcn_result = run_spagcn_like(
        cell_ids=cell_ids,
        genes=genes,
        coordinates=coordinates,
        histology_rgb=histology_rgb,
        expression=expression,
        base_scores=base_scores,
        base_best_indices=base_best_indices,
        domain_names=domain_names,
        neighbor_count=neighbor_count,
        params=payload["parameters"]["spagcn"],
        marker_top_n=marker_top_n,
    )
    graphst_result = run_graphst_like(
        cell_ids=cell_ids,
        genes=genes,
        coordinates=coordinates,
        histology_rgb=histology_rgb,
        expression=expression,
        base_scores=base_scores,
        base_best_indices=base_best_indices,
        domain_names=domain_names,
        neighbor_count=neighbor_count,
        params=payload["parameters"]["graphst"],
        marker_top_n=marker_top_n,
    )
    method_results = {
        spagcn_result.method: spagcn_result,
        graphst_result.method: graphst_result,
    }
    resolved_method = choose_selected_method(selected_method, method_results)

    if expected.get("expected_selected_method") is not None and resolved_method != expected["expected_selected_method"]:
        raise AssertionError("Toy input expected a different selected method than the computed result.")

    write_outputs(
        outdir=outdir,
        payload=payload,
        cell_ids=cell_ids,
        genes=genes,
        library_size=library_size,
        base_best_indices=base_best_indices,
        domain_names=domain_names,
        method_results=method_results,
        selected_method=resolved_method,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the spatial domain detection starter.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--selected-method", choices=["best", *METHODS], default="best")
    args = parser.parse_args()

    run_skill(args.input, args.outdir, args.selected_method)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
