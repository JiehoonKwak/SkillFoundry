#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
QC_FILES = [
    "qc_spatial_graph.tsv",
    "qc_neighborhood_profiles.tsv",
    "qc_niche_scores.tsv",
    "run_summary.json",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_metadata(skill_dir: Path) -> dict[str, Any]:
    return load_json(skill_dir / "metadata.yaml")


def sanitize_label(label: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in label.lower()).strip("_")


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


def validate_markdown_sections(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise AssertionError("Every row must have positive mass before normalization.")
    return matrix / totals


def softmax_rows(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def probability_margin(probabilities: np.ndarray) -> np.ndarray:
    top_two = np.partition(probabilities, -2, axis=1)[:, -2:]
    return top_two[:, 1] - top_two[:, 0]


def normalized_entropy(probabilities: np.ndarray) -> np.ndarray:
    safe = np.clip(probabilities, 1e-9, 1.0)
    entropy = -np.sum(safe * np.log(safe), axis=1)
    return entropy / math.log(probabilities.shape[1])


def load_toy_inputs(
    payload: dict[str, Any],
) -> tuple[list[str], list[str], list[str], list[str], np.ndarray, np.ndarray]:
    genes = payload["genes"]
    cell_types = payload["cell_types"]
    niche_names = payload["candidate_niches"]
    cells = payload["cells"]

    cell_ids = [cell["cell_id"] for cell in cells]
    if len(set(cell_ids)) != len(cell_ids):
        raise AssertionError("Cell identifiers must be unique.")

    coordinates = np.asarray([[cell["x"], cell["y"]] for cell in cells], dtype=float)
    expression = np.asarray(
        [[float(cell["expression"][gene]) for gene in genes] for cell in cells],
        dtype=float,
    )
    declared_cell_types = [cell["cell_type"] for cell in cells]
    sample_ids = [cell["sample_id"] for cell in cells]

    if np.any(expression <= 0):
        raise AssertionError("Toy expression values must be strictly positive for deterministic log-normalization.")
    if len(cell_ids) < 3:
        raise AssertionError("At least three cells are required to build a niche starter.")

    for cell in cells:
        if sorted(cell["expression"].keys()) != sorted(genes):
            raise AssertionError(f"Cell {cell['cell_id']} does not cover the full gene list.")
        if cell["cell_type"] not in cell_types:
            raise AssertionError(f"Unknown cell_type {cell['cell_type']} in toy input.")

    for niche_name in niche_names:
        prototype = payload["niche_prototypes"][niche_name]
        if sorted(prototype["composition"].keys()) != sorted(cell_types):
            raise AssertionError(f"Niche prototype {niche_name} must define every declared cell type.")

    return genes, cell_types, niche_names, sample_ids, coordinates, expression, declared_cell_types, cell_ids


def normalize_expression(expression: np.ndarray, scale_factor: float) -> tuple[np.ndarray, np.ndarray]:
    library_size = expression.sum(axis=1, keepdims=True)
    if np.any(library_size <= 0):
        raise AssertionError("Every cell must have positive expression mass.")
    normalized = np.log1p((expression / library_size) * float(scale_factor))
    return normalized, library_size[:, 0]


def build_spatial_graph(
    *,
    cell_ids: list[str],
    coordinates: np.ndarray,
    neighbor_count: int,
    self_weight: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    if not 0.0 <= self_weight < 1.0:
        raise AssertionError("self_weight must be in [0, 1).")

    neighbor_graph = np.zeros((len(cell_ids), len(cell_ids)), dtype=float)
    combined_graph = np.zeros((len(cell_ids), len(cell_ids)), dtype=float)
    rows: list[dict[str, Any]] = []
    capped_neighbors = max(1, min(int(neighbor_count), len(cell_ids) - 1))

    for source_index, source_cell in enumerate(cell_ids):
        distances = np.linalg.norm(coordinates[source_index] - coordinates, axis=1)
        order = [index for index in np.argsort(distances) if index != source_index][:capped_neighbors]
        raw_weights = np.asarray([1.0 / max(float(distances[index]), 1e-9) for index in order], dtype=float)
        raw_weights = raw_weights / raw_weights.sum()

        combined_graph[source_index, source_index] = float(self_weight)
        rows.append(
            {
                "source_cell": source_cell,
                "target_cell": source_cell,
                "relation": "self",
                "distance": 0.0,
                "weight": round(float(self_weight), 6),
            }
        )

        for target_index, normalized_weight in zip(order, raw_weights, strict=True):
            edge_weight = float((1.0 - self_weight) * normalized_weight)
            neighbor_graph[source_index, target_index] = float(normalized_weight)
            combined_graph[source_index, target_index] = edge_weight
            rows.append(
                {
                    "source_cell": source_cell,
                    "target_cell": cell_ids[target_index],
                    "relation": "neighbor",
                    "distance": round(float(distances[target_index]), 6),
                    "weight": round(edge_weight, 6),
                }
            )

    if not np.allclose(neighbor_graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Neighbor-only graph rows must sum to 1.")
    if not np.allclose(combined_graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Combined graph rows must sum to 1.")

    return neighbor_graph, combined_graph, rows


def compute_neighborhood_profiles(
    *,
    declared_cell_types: list[str],
    all_cell_types: list[str],
    neighbor_graph: np.ndarray,
    combined_graph: np.ndarray,
    smoothing_weight: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    cell_type_index = {cell_type: index for index, cell_type in enumerate(all_cell_types)}
    one_hot = np.zeros((len(declared_cell_types), len(all_cell_types)), dtype=float)
    for row_index, cell_type in enumerate(declared_cell_types):
        one_hot[row_index, cell_type_index[cell_type]] = 1.0

    local_profile = neighbor_graph @ one_hot
    smoothed_profile = (
        (1.0 - float(smoothing_weight)) * local_profile
        + float(smoothing_weight) * (combined_graph @ local_profile)
    )
    smoothed_profile = row_normalize(smoothed_profile)

    top_neighbor_fraction = local_profile.max(axis=1)
    boundary_score = normalized_entropy(smoothed_profile)
    dominant_neighbor_indices = local_profile.argmax(axis=1)
    dominant_neighbor_type = [all_cell_types[index] for index in dominant_neighbor_indices]
    return local_profile, smoothed_profile, top_neighbor_fraction, boundary_score, dominant_neighbor_type


def build_niche_score_tables(
    *,
    local_profile: np.ndarray,
    smoothed_profile: np.ndarray,
    boundary_score: np.ndarray,
    niche_names: list[str],
    cell_types: list[str],
    niche_prototypes: dict[str, dict[str, Any]],
    boundary_weight: float,
    temperature: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    composition_matrix = np.asarray(
        [
            [float(niche_prototypes[niche]["composition"][cell_type]) for cell_type in cell_types]
            for niche in niche_names
        ],
        dtype=float,
    )
    composition_matrix = row_normalize(composition_matrix)
    expected_boundary = np.asarray(
        [float(niche_prototypes[niche]["expected_boundary"]) for niche in niche_names],
        dtype=float,
    )

    base_scores = local_profile @ composition_matrix.T
    boundary_alignment = 1.0 - np.abs(boundary_score[:, None] - expected_boundary[None, :])
    final_scores = smoothed_profile @ composition_matrix.T + float(boundary_weight) * boundary_alignment
    probabilities = softmax_rows(final_scores * float(temperature))
    base_best = base_scores.argmax(axis=1)
    final_best = probabilities.argmax(axis=1)
    confidence_margin = probability_margin(probabilities)
    graph_shifted = base_best != final_best
    return base_scores, final_scores, probabilities, base_best, final_best, confidence_margin, graph_shifted


def note_for_cell(graph_shifted: bool, boundary_score: float) -> str:
    if graph_shifted:
        return "graph-smoothed niche differs from the raw neighbor mix"
    if boundary_score >= 0.72:
        return "high-boundary mixed neighborhood"
    if boundary_score <= 0.32:
        return "cohesive local neighborhood"
    return "moderately mixed local neighborhood"


def compute_marker_table(
    *,
    normalized_expression: np.ndarray,
    genes: list[str],
    niche_names: list[str],
    assigned_indices: np.ndarray,
    smoothed_profile: np.ndarray,
    cell_types: list[str],
    probabilities: np.ndarray,
    boundary_score: np.ndarray,
    niche_prototypes: dict[str, dict[str, Any]],
    top_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for niche_index, niche_name in enumerate(niche_names):
        in_niche = assigned_indices == niche_index
        out_niche = assigned_indices != niche_index
        if in_niche.sum() == 0 or out_niche.sum() == 0:
            raise AssertionError(f"Niche {niche_name} must have both in-niche and out-of-niche cells.")

        mean_in = normalized_expression[in_niche].mean(axis=0)
        mean_out = normalized_expression[out_niche].mean(axis=0)
        effect_sizes = mean_in - mean_out
        ordered = np.argsort(effect_sizes)[::-1]
        positive_genes = [genes[index] for index in ordered if effect_sizes[index] > 0]
        if not positive_genes:
            raise AssertionError(f"Niche {niche_name} has no positive marker effect in the toy starter.")

        mean_profile = smoothed_profile[in_niche].mean(axis=0)
        dominant_order = np.argsort(mean_profile)[::-1]
        dominant_types = [cell_types[index] for index in dominant_order if mean_profile[index] > 0]
        expected_boundary = float(niche_prototypes[niche_name]["expected_boundary"])
        top_effect = float(effect_sizes[ordered[0]])
        confidence_term = float(probabilities[in_niche, niche_index].mean())
        boundary_term = 1.0 - abs(float(boundary_score[in_niche].mean()) - expected_boundary)
        support_score = (
            0.5 * confidence_term
            + 0.3 * (1.0 - math.exp(-top_effect))
            + 0.2 * boundary_term
        )

        rows.append(
            {
                "niche_label": niche_name,
                "dominant_cell_types": ";".join(dominant_types[:2]),
                "marker_context": ";".join(positive_genes[: int(top_n)]),
                "support_score": round(float(min(max(support_score, 0.0), 1.0)), 6),
                "top_marker_effect": round(top_effect, 6),
                "mean_boundary_score": round(float(boundary_score[in_niche].mean()), 6),
                "mean_niche_confidence": round(confidence_term, 6),
            }
        )

    marker_table = pd.DataFrame(rows)
    return marker_table.sort_values("niche_label").reset_index(drop=True)


def validate_expected_invariants(
    payload: dict[str, Any],
    *,
    genes: list[str],
    cell_types: list[str],
    niche_names: list[str],
    cell_ids: list[str],
    graph_shifted: np.ndarray,
) -> None:
    expected = payload.get("expected_invariants", {})
    checks = {
        "cell_count": len(cell_ids),
        "gene_count": len(genes),
        "cell_type_count": len(cell_types),
        "niche_count": len(niche_names),
    }
    for key, observed in checks.items():
        if key in expected and int(expected[key]) != int(observed):
            raise AssertionError(f"Toy input invariant {key} does not match the raw starter inputs.")
    if "min_graph_shift_count" in expected and int(graph_shifted.sum()) < int(expected["min_graph_shift_count"]):
        raise AssertionError("The toy starter must include at least one graph-shifted niche assignment.")


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

    labels = pd.read_csv(outdir / "niche_labels.tsv", sep="\t")
    if labels["cell_id"].duplicated().any():
        raise AssertionError("niche_labels.tsv must have unique cell_id values.")
    if not labels["niche_confidence"].between(0.0, 1.0).all():
        raise AssertionError("niche_confidence must stay in [0, 1].")
    if not labels["boundary_score"].between(0.0, 1.0).all():
        raise AssertionError("boundary_score must stay in [0, 1].")

    markers = pd.read_csv(outdir / "niche_markers.tsv", sep="\t")
    if (markers["support_score"] <= 0).any():
        raise AssertionError("support_score must stay positive for every toy niche.")

    graph = pd.read_csv(outdir / "qc_spatial_graph.tsv", sep="\t")
    row_sums = graph.groupby("source_cell")["weight"].sum().to_numpy(dtype=float)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise AssertionError("Each source row in qc_spatial_graph.tsv must sum to 1.")
    if set(graph["relation"]) != {"neighbor", "self"}:
        raise AssertionError("qc_spatial_graph.tsv must contain self and neighbor rows.")

    neighborhood = pd.read_csv(outdir / "qc_neighborhood_profiles.tsv", sep="\t")
    local_columns = [column for column in neighborhood.columns if column.startswith("local_") and column.endswith("_fraction")]
    smoothed_columns = [column for column in neighborhood.columns if column.startswith("smoothed_") and column.endswith("_fraction")]
    if not np.allclose(neighborhood[local_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6):
        raise AssertionError("Local neighborhood fractions must sum to 1 for each cell.")
    if not np.allclose(neighborhood[smoothed_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6):
        raise AssertionError("Smoothed neighborhood fractions must sum to 1 for each cell.")

    scores = pd.read_csv(outdir / "qc_niche_scores.tsv", sep="\t")
    probability_columns = [column for column in scores.columns if column.endswith("_probability")]
    if not np.allclose(scores[probability_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6):
        raise AssertionError("Niche probabilities must sum to 1 in qc_niche_scores.tsv.")
    label_rows = labels.set_index("cell_id")
    score_rows = scores.set_index("cell_id")
    if list(label_rows.index) != list(score_rows.index):
        raise AssertionError("niche_labels.tsv and qc_niche_scores.tsv must share the same cell ordering.")
    if list(label_rows["niche_label"]) != list(score_rows["assigned_niche"]):
        raise AssertionError("niche_labels.tsv must match the assigned niche calls in qc_niche_scores.tsv.")

    summary = load_json(outdir / "run_summary.json")
    written_files = set(summary.get("written_files", []))
    expected_files = set([item["path"] for item in metadata["deliverables"]] + metadata.get("starter_qc_files", QC_FILES))
    if expected_files != written_files:
        raise AssertionError("run_summary.json written_files must list all deliverables and QC files.")


def run_skill(input_path: Path, outdir: Path) -> None:
    payload = load_json(input_path)
    (
        genes,
        cell_types,
        niche_names,
        sample_ids,
        coordinates,
        expression,
        declared_cell_types,
        cell_ids,
    ) = load_toy_inputs(payload)
    parameters = payload["parameters"]
    normalized_expression, library_size = normalize_expression(expression, float(parameters["expression_scale"]))

    neighbor_graph, combined_graph, graph_rows = build_spatial_graph(
        cell_ids=cell_ids,
        coordinates=coordinates,
        neighbor_count=int(parameters["neighbor_count"]),
        self_weight=float(parameters["self_weight"]),
    )
    local_profile, smoothed_profile, top_neighbor_fraction, boundary_score, dominant_neighbor_type = compute_neighborhood_profiles(
        declared_cell_types=declared_cell_types,
        all_cell_types=cell_types,
        neighbor_graph=neighbor_graph,
        combined_graph=combined_graph,
        smoothing_weight=float(parameters["smoothing_weight"]),
    )
    (
        base_scores,
        final_scores,
        probabilities,
        base_best,
        final_best,
        confidence_margin,
        graph_shifted,
    ) = build_niche_score_tables(
        local_profile=local_profile,
        smoothed_profile=smoothed_profile,
        boundary_score=boundary_score,
        niche_names=niche_names,
        cell_types=cell_types,
        niche_prototypes=payload["niche_prototypes"],
        boundary_weight=float(parameters["boundary_weight"]),
        temperature=float(parameters["temperature"]),
    )
    validate_expected_invariants(
        payload,
        genes=genes,
        cell_types=cell_types,
        niche_names=niche_names,
        cell_ids=cell_ids,
        graph_shifted=graph_shifted,
    )

    outdir.mkdir(parents=True, exist_ok=True)

    neighborhood_rows = []
    for row_index, cell_id in enumerate(cell_ids):
        row: dict[str, Any] = {
            "cell_id": cell_id,
            "cell_type": declared_cell_types[row_index],
            "dominant_neighbor_type": dominant_neighbor_type[row_index],
            "top_neighbor_fraction": round(float(top_neighbor_fraction[row_index]), 6),
            "boundary_score": round(float(boundary_score[row_index]), 6),
            "graph_shifted": bool(graph_shifted[row_index]),
        }
        for cell_type_index, cell_type in enumerate(cell_types):
            slug = sanitize_label(cell_type)
            row[f"local_{slug}_fraction"] = round(float(local_profile[row_index, cell_type_index]), 6)
            row[f"smoothed_{slug}_fraction"] = round(float(smoothed_profile[row_index, cell_type_index]), 6)
        neighborhood_rows.append(row)
    neighborhood_frame = pd.DataFrame(neighborhood_rows)
    neighborhood_frame.to_csv(outdir / "qc_neighborhood_profiles.tsv", sep="\t", index=False)

    score_rows = []
    for row_index, cell_id in enumerate(cell_ids):
        row = {
            "cell_id": cell_id,
            "cell_type": declared_cell_types[row_index],
            "base_best_niche": niche_names[base_best[row_index]],
            "assigned_niche": niche_names[final_best[row_index]],
            "graph_shifted": bool(graph_shifted[row_index]),
            "boundary_score": round(float(boundary_score[row_index]), 6),
            "confidence_margin": round(float(confidence_margin[row_index]), 6),
        }
        for niche_index, niche_name in enumerate(niche_names):
            slug = sanitize_label(niche_name)
            row[f"base_{slug}_score"] = round(float(base_scores[row_index, niche_index]), 6)
            row[f"final_{slug}_score"] = round(float(final_scores[row_index, niche_index]), 6)
            row[f"{slug}_probability"] = round(float(probabilities[row_index, niche_index]), 6)
        score_rows.append(row)
    score_frame = pd.DataFrame(score_rows)
    score_frame.to_csv(outdir / "qc_niche_scores.tsv", sep="\t", index=False)

    niche_labels = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "sample_id": sample_ids,
            "cell_type": declared_cell_types,
            "niche_label": [niche_names[index] for index in final_best],
            "niche_confidence": np.round(probabilities.max(axis=1), 6),
            "dominant_neighbor_type": dominant_neighbor_type,
            "boundary_score": np.round(boundary_score, 6),
            "notes": [note_for_cell(bool(shifted), float(boundary)) for shifted, boundary in zip(graph_shifted, boundary_score, strict=True)],
            "base_best_niche": [niche_names[index] for index in base_best],
            "graph_shifted": graph_shifted.astype(bool),
            "confidence_margin": np.round(confidence_margin, 6),
            "library_size": np.round(library_size, 6),
        }
    )
    niche_labels.to_csv(outdir / "niche_labels.tsv", sep="\t", index=False)

    marker_table = compute_marker_table(
        normalized_expression=normalized_expression,
        genes=genes,
        niche_names=niche_names,
        assigned_indices=final_best,
        smoothed_profile=smoothed_profile,
        cell_types=cell_types,
        probabilities=probabilities,
        boundary_score=boundary_score,
        niche_prototypes=payload["niche_prototypes"],
        top_n=int(parameters["marker_top_n"]),
    )
    marker_table.to_csv(outdir / "niche_markers.tsv", sep="\t", index=False)

    pd.DataFrame(graph_rows).to_csv(outdir / "qc_spatial_graph.tsv", sep="\t", index=False)

    niche_counts = niche_labels["niche_label"].value_counts().sort_index().to_dict()
    shifted_cells = niche_labels.loc[niche_labels["graph_shifted"], "cell_id"].tolist()
    boundary_cells = niche_labels.loc[niche_labels["boundary_score"] >= 0.70, "cell_id"].tolist()
    epithelial_multi_niche = sorted(
        niche_labels.loc[niche_labels["cell_type"] == "Epithelial", "niche_label"].unique().tolist()
    )

    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}",
                f"Toy input size: {len(cell_ids)} cells, {len(genes)} genes, {len(niche_names)} candidate niches.",
                "Starter path: spatial neighbor graph, neighborhood composition, graph-smoothed niche scoring, and niche marker summarization.",
            ],
        },
        {
            "name": "Spatial graph",
            "bullets": [
                f"kNN graph used {int(parameters['neighbor_count'])} neighbors per cell with self_weight={float(parameters['self_weight']):.2f}.",
                (
                    "Graph smoothing changed at least one niche call: " + ", ".join(shifted_cells)
                    if shifted_cells
                    else "Graph smoothing did not alter the raw composition ranking on this run."
                ),
                (
                    "High-boundary interface candidates: " + ", ".join(boundary_cells)
                    if boundary_cells
                    else "No cells crossed the high-boundary threshold."
                ),
            ],
        },
        {
            "name": "Niche definitions",
            "bullets": [
                "Niche counts: " + ", ".join(f"{label}={count}" for label, count in sorted(niche_counts.items())),
                (
                    "Epithelial cells split across neighborhood-defined niches: "
                    + ", ".join(epithelial_multi_niche)
                ),
                "Markers: "
                + "; ".join(
                    f"{row.niche_label} -> {row.marker_context} [{row.dominant_cell_types}]"
                    for row in marker_table.itertuples()
                ),
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "The starter uses deterministic prototype scoring over local cell-type composition instead of a trained CellCharter, Squidpy permutation test, or UTAG model.",
                "Boundary calls come from graph-smoothed neighborhood entropy rather than a formal statistical interface model.",
                "Marker summaries are effect-size rankings on tiny synthetic inputs and are only a portable surrogate for real niche interpretation.",
            ],
        },
    ]
    write_markdown(outdir / "tissue_niche_report.md", "Tissue Niche Report", report_sections)

    summary = {
        "boundary_cells": boundary_cells,
        "cell_count": len(cell_ids),
        "gene_count": len(genes),
        "method_steps": [
            "spatial_graph_construction",
            "neighborhood_composition",
            "graph_smoothed_niche_scoring",
            "marker_summarization",
        ],
        "niche_counts": niche_counts,
        "run_label": payload["run_label"],
        "shifted_cells": shifted_cells,
        "written_files": [
            "niche_labels.tsv",
            "niche_markers.tsv",
            "tissue_niche_report.md",
            "qc_spatial_graph.tsv",
            "qc_neighborhood_profiles.tsv",
            "qc_niche_scores.tsv",
            "run_summary.json",
        ],
    }
    write_json(outdir / "run_summary.json", summary)

    validate_outputs(SKILL_DIR, outdir)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the tissue niche annotation toy starter.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    run_skill(args.input, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
