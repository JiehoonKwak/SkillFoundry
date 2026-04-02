#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
METHOD_STEPS = [
    "spatial_neighbors",
    "nhood_enrichment",
    "co_occurrence",
    "morans_i",
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


def sanitize_label(label: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in label.lower()).strip("_")


def softmax_rows(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - matrix.max(axis=1, keepdims=True)
    exp_scores = np.exp(shifted)
    return exp_scores / exp_scores.sum(axis=1, keepdims=True)


def probability_margin(probabilities: np.ndarray) -> np.ndarray:
    top_two = np.sort(probabilities, axis=1)[:, -2:]
    return top_two[:, 1] - top_two[:, 0]


def load_payload(
    input_path: Path,
) -> tuple[dict[str, Any], list[str], list[str], list[str], list[str], np.ndarray, np.ndarray, np.ndarray]:
    payload = load_json(input_path)
    genes = payload["genes"]
    image_features = payload["image_features"]
    domains = payload["candidate_domains"]
    spots = payload["spots"]

    cell_ids = [spot["cell_id"] for spot in spots]
    coordinates = np.asarray([[spot["x"], spot["y"]] for spot in spots], dtype=float)
    counts = np.asarray([[spot["counts"][gene] for gene in genes] for spot in spots], dtype=float)
    image_values = np.asarray(
        [[spot["image_features"][feature] for feature in image_features] for spot in spots],
        dtype=float,
    )

    if len(set(cell_ids)) != len(cell_ids):
        raise AssertionError("Cell identifiers must be unique.")
    if np.any(counts <= 0):
        raise AssertionError("Toy counts must be strictly positive for deterministic normalization.")
    if np.any(image_values < 0):
        raise AssertionError("Image feature summaries must be non-negative.")
    if counts.shape[0] < 4:
        raise AssertionError("At least four spots are required for the starter graph.")

    return payload, genes, image_features, domains, cell_ids, counts, coordinates, image_values


def normalize_counts(counts: np.ndarray, scale_factor: float) -> tuple[np.ndarray, np.ndarray]:
    library_size = counts.sum(axis=1, keepdims=True)
    if np.any(library_size <= 0):
        raise AssertionError("Every spot must have a positive library size.")
    normalized = np.log1p((counts / library_size) * float(scale_factor))
    return normalized, library_size[:, 0]


def normalize_image_features(image_values: np.ndarray) -> np.ndarray:
    row_sums = image_values.sum(axis=1, keepdims=True)
    if np.any(row_sums <= 0):
        raise AssertionError("Each spot needs at least one positive image feature summary.")
    return image_values / row_sums


def build_profile_matrix(names: list[str], keys: list[str], payload: dict[str, dict[str, float]]) -> np.ndarray:
    matrix = np.asarray(
        [[float(payload[name][key]) for key in keys] for name in names],
        dtype=float,
    )
    if np.any(np.isclose(matrix.sum(axis=1), 0.0)):
        raise AssertionError("Each profile needs at least one non-zero weight.")
    return matrix / matrix.sum(axis=1, keepdims=True)


def build_spatial_graph(
    *,
    cell_ids: list[str],
    coordinates: np.ndarray,
    image_values: np.ndarray,
    neighbor_count: int,
    self_weight: float,
    image_alpha: float,
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
            image_distance = float(np.linalg.norm(image_values[source_index] - image_values[target_index]))
            raw_weights.append(spatial_term * math.exp(-float(image_alpha) * image_distance))

        normalized_weights = np.asarray(raw_weights, dtype=float)
        normalized_weights = normalized_weights / normalized_weights.sum()

        graph[source_index, source_index] = float(self_weight)
        rows.append(
            {
                "source_cell": source_cell,
                "target_cell": source_cell,
                "relation": "self",
                "distance": 0.0,
                "image_distance": 0.0,
                "weight": round(float(self_weight), 6),
            }
        )

        for target_index, normalized_weight in zip(order, normalized_weights, strict=True):
            distance = float(distances[target_index])
            image_distance = float(np.linalg.norm(image_values[source_index] - image_values[target_index]))
            edge_weight = float((1.0 - self_weight) * normalized_weight)
            graph[source_index, target_index] = edge_weight
            rows.append(
                {
                    "source_cell": source_cell,
                    "target_cell": cell_ids[target_index],
                    "relation": "neighbor",
                    "distance": round(distance, 6),
                    "image_distance": round(image_distance, 6),
                    "weight": round(edge_weight, 6),
                }
            )

    if not np.allclose(graph.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Graph rows must sum to 1.")

    return graph, rows


def score_domains(
    *,
    normalized_counts: np.ndarray,
    normalized_image: np.ndarray,
    domain_signatures: np.ndarray,
    domain_image_profiles: np.ndarray,
    expression_weight: float,
) -> tuple[np.ndarray, np.ndarray]:
    expression_scores = normalized_counts @ domain_signatures.T
    image_scores = normalized_image @ domain_image_profiles.T
    base_scores = float(expression_weight) * expression_scores + (1.0 - float(expression_weight)) * image_scores
    return base_scores, image_scores


def smooth_scores(base_scores: np.ndarray, graph: np.ndarray, rounds: int) -> np.ndarray:
    smoothed = np.asarray(base_scores, dtype=float)
    for _ in range(max(1, int(rounds))):
        smoothed = graph @ smoothed
    return smoothed


def build_label_qc(
    *,
    cell_ids: list[str],
    domains: list[str],
    base_scores: np.ndarray,
    smoothed_scores: np.ndarray,
) -> tuple[pd.DataFrame, list[str], list[str], np.ndarray, np.ndarray]:
    probabilities = softmax_rows(smoothed_scores)
    confidence_margin = probability_margin(probabilities)
    base_indices = np.argmax(base_scores, axis=1)
    assigned_indices = np.argmax(smoothed_scores, axis=1)
    base_labels = [domains[index] for index in base_indices]
    assigned_labels = [domains[index] for index in assigned_indices]
    assigned_scores = probabilities[np.arange(len(cell_ids)), assigned_indices]

    rows: list[dict[str, Any]] = []
    for index, cell_id in enumerate(cell_ids):
        row: dict[str, Any] = {
            "cell_id": cell_id,
            "base_best_domain": base_labels[index],
            "assigned_label": assigned_labels[index],
            "domain_score": round(float(assigned_scores[index]), 6),
            "confidence_margin": round(float(confidence_margin[index]), 6),
            "graph_shifted": bool(base_labels[index] != assigned_labels[index]),
        }
        for domain_index, domain in enumerate(domains):
            suffix = sanitize_label(domain)
            row[f"base_score_{suffix}"] = round(float(base_scores[index, domain_index]), 6)
            row[f"smoothed_score_{suffix}"] = round(float(smoothed_scores[index, domain_index]), 6)
        rows.append(row)

    return pd.DataFrame(rows), base_labels, assigned_labels, assigned_scores, confidence_margin


def compute_neighborhood_enrichment(
    *,
    graph: np.ndarray,
    labels: list[str],
    domains: list[str],
    pseudocount: float,
) -> pd.DataFrame:
    nonself_graph = np.asarray(graph, dtype=float).copy()
    np.fill_diagonal(nonself_graph, 0.0)

    counts = {domain: labels.count(domain) for domain in domains}
    prevalence = {domain: counts[domain] / len(labels) for domain in domains}
    rows: list[dict[str, Any]] = []

    for label_a in domains:
        source_indices = [index for index, label in enumerate(labels) if label == label_a]
        total_outgoing = float(nonself_graph[source_indices].sum()) if source_indices else 0.0

        for label_b in domains:
            target_indices = [index for index, label in enumerate(labels) if label == label_b]
            observed_weight = (
                float(nonself_graph[np.ix_(source_indices, target_indices)].sum())
                if source_indices and target_indices
                else 0.0
            )
            expected_weight = total_outgoing * prevalence[label_b]
            enrichment_score = math.log2((observed_weight + pseudocount) / (expected_weight + pseudocount))
            rows.append(
                {
                    "label_a": label_a,
                    "label_b": label_b,
                    "source_count": counts[label_a],
                    "target_count": counts[label_b],
                    "observed_weight": round(observed_weight, 6),
                    "expected_weight": round(expected_weight, 6),
                    "enrichment_score": round(enrichment_score, 6),
                }
            )

    return pd.DataFrame(rows)


def compute_co_occurrence(
    *,
    coordinates: np.ndarray,
    labels: list[str],
    domains: list[str],
    distance_threshold: float,
) -> pd.DataFrame:
    domain_order = {domain: index for index, domain in enumerate(domains)}
    near_total = 0
    total_pairs = 0
    pair_counts = {
        (label_a, label_b): {"near_pairs": 0, "total_pairs": 0}
        for label_a in domains
        for label_b in domains
    }

    for left_index, right_index in itertools.combinations(range(len(labels)), 2):
        distance = float(np.linalg.norm(coordinates[left_index] - coordinates[right_index]))
        near = distance <= float(distance_threshold)
        total_pairs += 1
        near_total += int(near)

        label_left = labels[left_index]
        label_right = labels[right_index]
        if domain_order[label_left] <= domain_order[label_right]:
            pair_key = (label_left, label_right)
        else:
            pair_key = (label_right, label_left)
        pair_counts[pair_key]["total_pairs"] += 1
        pair_counts[pair_key]["near_pairs"] += int(near)

    baseline_near_fraction = near_total / total_pairs if total_pairs else 0.0
    rows: list[dict[str, Any]] = []
    for label_a in domains:
        for label_b in domains:
            pair_key = (label_a, label_b) if domain_order[label_a] <= domain_order[label_b] else (label_b, label_a)
            total_pair_count = pair_counts[pair_key]["total_pairs"]
            near_pair_count = pair_counts[pair_key]["near_pairs"]
            near_fraction = near_pair_count / total_pair_count if total_pair_count else 0.0
            co_occurrence_score = near_fraction / baseline_near_fraction if baseline_near_fraction else 0.0
            rows.append(
                {
                    "label_a": label_a,
                    "label_b": label_b,
                    "near_pairs": near_pair_count,
                    "total_pairs": total_pair_count,
                    "near_fraction": round(near_fraction, 6),
                    "baseline_near_fraction": round(baseline_near_fraction, 6),
                    "co_occurrence_score": round(co_occurrence_score, 6),
                }
            )

    return pd.DataFrame(rows)


def morans_i(graph: np.ndarray, values: np.ndarray) -> float:
    centered = values - float(np.mean(values))
    denominator = float(np.square(centered).sum())
    if denominator <= 0.0:
        raise AssertionError("Moran's I requires non-constant values.")
    weight_sum = float(graph.sum())
    if weight_sum <= 0.0:
        raise AssertionError("Moran's I requires at least one non-self graph edge.")
    numerator = float((graph * np.outer(centered, centered)).sum())
    return (len(values) / weight_sum) * (numerator / denominator)


def compute_moran_table(
    *,
    graph: np.ndarray,
    normalized_counts: np.ndarray,
    normalized_image: np.ndarray,
    genes: list[str],
    image_features: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for index, gene in enumerate(genes):
        rows.append(
            {
                "statistic": "moran_i",
                "scope": f"gene:{gene}",
                "value": round(float(morans_i(graph, normalized_counts[:, index])), 6),
                "feature_type": "gene",
            }
        )

    for index, feature in enumerate(image_features):
        rows.append(
            {
                "statistic": "moran_i",
                "scope": f"image:{feature}",
                "value": round(float(morans_i(graph, normalized_image[:, index])), 6),
                "feature_type": "image_feature",
            }
        )

    return pd.DataFrame(rows)


def build_report_sections(
    *,
    payload: dict[str, Any],
    derived_counts: dict[str, int],
    shifted_cells: list[str],
    top_enrichment: pd.Series,
    top_co_occurrence: pd.Series,
    top_moran: pd.Series,
) -> list[dict[str, Any]]:
    parameters = payload["parameters"]
    return [
        {
            "name": "Run context",
            "bullets": [
                (
                    f"{len(payload['spots'])} synthetic spots, {len(payload['genes'])} genes, and "
                    f"{len(payload['image_features'])} per-spot image summary features were analyzed."
                ),
                "Spot domains were derived from graph-smoothed signature scores; the toy input does not provide final labels.",
                "This starter computes spatial neighbors, neighborhood enrichment, co-occurrence, and Moran's I without network access.",
            ],
        },
        {
            "name": "Spatial graph",
            "bullets": [
                (
                    "Built an image-aware KNN graph with "
                    f"neighbor_count={parameters['neighbor_count']}, self_weight={parameters['self_weight']}, "
                    f"image_alpha={parameters['image_alpha']}."
                ),
                f"Derived domain counts: {', '.join(f'{label}={count}' for label, count in derived_counts.items())}.",
                (
                    "Boundary spots shifted after graph smoothing: "
                    f"{', '.join(shifted_cells) if shifted_cells else 'none'}."
                ),
            ],
        },
        {
            "name": "Statistics",
            "bullets": [
                (
                    "Top neighborhood enrichment surrogate: "
                    f"{top_enrichment['label_a']}|{top_enrichment['label_b']}="
                    f"{float(top_enrichment['enrichment_score']):.3f}."
                ),
                (
                    "Top co-occurrence score at the near-distance threshold: "
                    f"{top_co_occurrence['label_a']}|{top_co_occurrence['label_b']}="
                    f"{float(top_co_occurrence['co_occurrence_score']):.3f}."
                ),
                f"Highest Moran's I: {top_moran['scope']}={float(top_moran['value']):.3f}.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "Neighborhood enrichment is a deterministic log2 observed-versus-expected surrogate, not Squidpy's permutation z-score.",
                "Image-aware behavior uses precomputed per-spot summary features instead of a real Squidpy ImageContainer or segmentation workflow.",
                "Moran's I values come from a tiny synthetic graph and do not include significance testing or multiple-testing correction.",
            ],
        },
    ]


def run_analysis(skill_dir: Path, input_path: Path, outdir: Path) -> dict[str, Any]:
    metadata = load_metadata(skill_dir)
    (
        payload,
        genes,
        image_features,
        domains,
        cell_ids,
        counts,
        coordinates,
        image_values,
    ) = load_payload(input_path)
    parameters = payload["parameters"]

    normalized_counts, library_sizes = normalize_counts(counts, parameters["expression_scale"])
    normalized_image = normalize_image_features(image_values)
    domain_signatures = build_profile_matrix(domains, genes, payload["domain_signatures"])
    domain_image_profiles = build_profile_matrix(domains, image_features, payload["domain_image_profiles"])

    graph, graph_rows = build_spatial_graph(
        cell_ids=cell_ids,
        coordinates=coordinates,
        image_values=normalized_image,
        neighbor_count=parameters["neighbor_count"],
        self_weight=parameters["self_weight"],
        image_alpha=parameters["image_alpha"],
    )
    base_scores, image_scores = score_domains(
        normalized_counts=normalized_counts,
        normalized_image=normalized_image,
        domain_signatures=domain_signatures,
        domain_image_profiles=domain_image_profiles,
        expression_weight=parameters["expression_weight"],
    )
    smoothed_scores = smooth_scores(base_scores, graph, parameters["smoothing_rounds"])
    label_qc, base_labels, assigned_labels, assigned_scores, confidence_margin = build_label_qc(
        cell_ids=cell_ids,
        domains=domains,
        base_scores=base_scores,
        smoothed_scores=smoothed_scores,
    )

    enrichment = compute_neighborhood_enrichment(
        graph=graph,
        labels=assigned_labels,
        domains=domains,
        pseudocount=parameters["enrichment_pseudocount"],
    )
    co_occurrence = compute_co_occurrence(
        coordinates=coordinates,
        labels=assigned_labels,
        domains=domains,
        distance_threshold=parameters["co_occurrence_distance"],
    )
    nonself_graph = graph.copy()
    np.fill_diagonal(nonself_graph, 0.0)
    moran_table = compute_moran_table(
        graph=nonself_graph,
        normalized_counts=normalized_counts,
        normalized_image=normalized_image,
        genes=genes,
        image_features=image_features,
    )

    pairwise_qc = enrichment.merge(co_occurrence, on=["label_a", "label_b"], how="inner")
    neighbor_deliverable = enrichment[
        [
            "label_a",
            "label_b",
            "enrichment_score",
            "observed_weight",
            "expected_weight",
            "source_count",
            "target_count",
        ]
    ].copy()

    co_occurrence_stats = co_occurrence.rename(columns={"co_occurrence_score": "value"}).copy()
    co_occurrence_stats.insert(0, "statistic", "co_occurrence")
    co_occurrence_stats.insert(
        1,
        "scope",
        co_occurrence_stats["label_a"] + "|" + co_occurrence_stats["label_b"],
    )
    spatial_stats = pd.concat(
        [
            co_occurrence_stats[
                [
                    "statistic",
                    "scope",
                    "value",
                    "near_pairs",
                    "total_pairs",
                    "near_fraction",
                    "baseline_near_fraction",
                ]
            ],
            moran_table[["statistic", "scope", "value", "feature_type"]],
        ],
        ignore_index=True,
        sort=False,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(graph_rows).to_csv(outdir / "qc_spatial_graph.tsv", sep="\t", index=False)

    label_qc = label_qc.copy()
    label_qc["library_size"] = np.round(library_sizes, 6)
    label_qc["assigned_domain_score"] = np.round(assigned_scores, 6)
    label_qc["confidence_margin"] = np.round(confidence_margin, 6)
    for feature_index, feature in enumerate(image_features):
        label_qc[f"image_feature_{sanitize_label(feature)}"] = np.round(normalized_image[:, feature_index], 6)
    for domain_index, domain in enumerate(domains):
        label_qc[f"image_profile_score_{sanitize_label(domain)}"] = np.round(image_scores[:, domain_index], 6)
    label_qc.to_csv(outdir / "qc_label_scores.tsv", sep="\t", index=False)
    pairwise_qc.to_csv(outdir / "qc_pairwise_stats.tsv", sep="\t", index=False)
    neighbor_deliverable.to_csv(outdir / "neighbor_enrichment.tsv", sep="\t", index=False)
    spatial_stats.to_csv(outdir / "spatial_stats.tsv", sep="\t", index=False)

    shifted_cells = sorted(label_qc.loc[label_qc["graph_shifted"], "cell_id"].tolist())
    derived_counts = {domain: int(assigned_labels.count(domain)) for domain in domains}
    top_enrichment = pairwise_qc.sort_values(["enrichment_score", "label_a", "label_b"], ascending=[False, True, True]).iloc[0]
    top_co_occurrence = pairwise_qc.sort_values(
        ["co_occurrence_score", "label_a", "label_b"],
        ascending=[False, True, True],
    ).iloc[0]
    top_moran = moran_table.sort_values(["value", "scope"], ascending=[False, True]).iloc[0]

    write_markdown(
        outdir / "squidpy_report.md",
        title="Squidpy Report",
        sections=build_report_sections(
            payload=payload,
            derived_counts=derived_counts,
            shifted_cells=shifted_cells,
            top_enrichment=top_enrichment,
            top_co_occurrence=top_co_occurrence,
            top_moran=top_moran,
        ),
    )

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path),
        "outdir": str(outdir),
        "method_steps": METHOD_STEPS,
        "spot_count": len(cell_ids),
        "gene_count": len(genes),
        "image_feature_count": len(image_features),
        "domain_count": len(domains),
        "graph_mode": "image_aware_knn",
        "label_assignment_mode": "graph_smoothed_signature_scores",
        "derived_domain_counts": derived_counts,
        "boundary_spots": shifted_cells,
        "top_neighbor_enrichment": {
            "label_a": str(top_enrichment["label_a"]),
            "label_b": str(top_enrichment["label_b"]),
            "enrichment_score": round(float(top_enrichment["enrichment_score"]), 6),
        },
        "top_co_occurrence": {
            "label_a": str(top_co_occurrence["label_a"]),
            "label_b": str(top_co_occurrence["label_b"]),
            "co_occurrence_score": round(float(top_co_occurrence["co_occurrence_score"]), 6),
        },
        "top_moran": {
            "scope": str(top_moran["scope"]),
            "value": round(float(top_moran["value"]), 6),
        },
        "written_files": [],
    }
    summary["written_files"] = sorted(
        set([item["path"] for item in metadata["deliverables"]] + metadata["starter_qc_files"] + ["run_summary.json"])
    )
    write_json(outdir / "run_summary.json", summary)

    validate_outputs(skill_dir, outdir, input_path)
    return summary


def validate_outputs(skill_dir: Path, outdir: Path, input_path: Path | None = None) -> None:
    metadata = load_metadata(skill_dir)
    payload = load_json(input_path or skill_dir / "examples" / "toy_input.json")

    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable {deliverable['path']}")
        if deliverable["kind"] == "tsv":
            table = pd.read_csv(path, sep="\t")
            if table.empty:
                raise AssertionError(f"{deliverable['path']} must not be empty.")
            missing = set(deliverable["required_columns"]) - set(table.columns)
            if missing:
                raise AssertionError(f"{deliverable['path']} is missing columns {sorted(missing)}")
        elif deliverable["kind"] == "md":
            validate_markdown_sections(path, deliverable["required_sections"])

    for qc_file in metadata["starter_qc_files"]:
        if not (outdir / qc_file).exists():
            raise AssertionError(f"Missing QC file {qc_file}")

    graph = pd.read_csv(outdir / "qc_spatial_graph.tsv", sep="\t")
    labels = pd.read_csv(outdir / "qc_label_scores.tsv", sep="\t")
    pairwise = pd.read_csv(outdir / "qc_pairwise_stats.tsv", sep="\t")
    stats = pd.read_csv(outdir / "spatial_stats.tsv", sep="\t")
    summary = load_json(outdir / "run_summary.json")

    expected_spot_count = payload["expected_invariants"]["spot_count"]
    domains = payload["candidate_domains"]
    expected_neighbors = payload["parameters"]["neighbor_count"] + 1

    if labels["cell_id"].nunique() != expected_spot_count:
        raise AssertionError("QC label table must have one unique row per spot.")
    if labels.shape[0] != expected_spot_count:
        raise AssertionError("QC label table row count must match the toy input.")
    if sorted(labels["assigned_label"].unique()) != sorted(domains):
        raise AssertionError("Derived labels must cover the declared domains.")
    if not labels["domain_score"].between(0.0, 1.0).all():
        raise AssertionError("Assigned domain scores must be probabilities in [0, 1].")

    shifted_cells = sorted(labels.loc[labels["graph_shifted"], "cell_id"].tolist())
    if shifted_cells != sorted(payload["expected_invariants"]["boundary_spots"]):
        raise AssertionError("Shifted spots do not match the toy invariants.")

    row_sums = graph.groupby("source_cell")["weight"].sum().to_numpy(dtype=float)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        raise AssertionError("Spatial graph rows must sum to 1.")
    if not graph["relation"].isin({"self", "neighbor"}).all():
        raise AssertionError("Spatial graph relations must be self or neighbor.")
    if not (graph.groupby("source_cell").size() == expected_neighbors).all():
        raise AssertionError("Each source cell must have one self edge plus the configured neighbors.")

    if pairwise.shape[0] != len(domains) ** 2:
        raise AssertionError("Pairwise QC must cover every ordered domain pair.")
    if not {"enrichment_score", "co_occurrence_score"}.issubset(pairwise.columns):
        raise AssertionError("Pairwise QC must include enrichment and co-occurrence scores.")

    statistic_names = set(stats["statistic"])
    if statistic_names != {"co_occurrence", "moran_i"}:
        raise AssertionError("Spatial stats must include only co_occurrence and moran_i rows.")
    if int((stats["statistic"] == "co_occurrence").sum()) != len(domains) ** 2:
        raise AssertionError("Spatial stats must keep one co-occurrence row per ordered domain pair.")
    if int((stats["statistic"] == "moran_i").sum()) != len(payload["genes"]) + len(payload["image_features"]):
        raise AssertionError("Spatial stats Moran rows must cover all genes and image features.")

    derived_counts = labels["assigned_label"].value_counts().to_dict()
    if summary["method_steps"] != METHOD_STEPS:
        raise AssertionError("run_summary.json must record the exact starter method order.")
    if summary["derived_domain_counts"] != derived_counts:
        raise AssertionError("run_summary.json derived domain counts must match QC labels.")
    if summary["boundary_spots"] != shifted_cells:
        raise AssertionError("run_summary.json boundary spots must match shifted QC rows.")
    if summary["spot_count"] != expected_spot_count:
        raise AssertionError("run_summary.json spot_count mismatch.")
    expected_written_files = sorted(
        set([item["path"] for item in metadata["deliverables"]] + metadata["starter_qc_files"] + ["run_summary.json"])
    )
    if sorted(summary["written_files"]) != expected_written_files:
        raise AssertionError("run_summary.json written_files must list all deliverables, QC files, and run_summary.json.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    run_analysis(SKILL_DIR, args.input, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
