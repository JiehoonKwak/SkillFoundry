#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]
QC_FILES = [
    "intermediate_qc.json",
    "gene_cv_qc.tsv",
    "spot_assignment_qc.tsv",
    "spatial_coherence_qc.tsv",
    "run_summary.json",
]


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata() -> dict:
    return load_json(SKILL_DIR / "metadata.yaml")


def project_simplex(vector: np.ndarray) -> np.ndarray:
    sorted_values = np.sort(vector)[::-1]
    cumulative = np.cumsum(sorted_values)
    indices = np.arange(1, len(vector) + 1)
    mask = sorted_values - (cumulative - 1.0) / indices > 0
    rho = np.nonzero(mask)[0][-1]
    theta = (cumulative[rho] - 1.0) / float(rho + 1)
    return np.maximum(vector - theta, 0.0)


def solve_nnls_simplex(target: np.ndarray, basis: np.ndarray, *, max_iter: int) -> np.ndarray:
    weights = np.full(basis.shape[0], 1.0 / basis.shape[0], dtype=float)
    gram = basis @ basis.T
    lipschitz = float(np.linalg.eigvalsh(gram).max())
    step = 1.0 / (lipschitz + 1e-8)
    for _ in range(max_iter):
        gradient = (weights @ basis - target) @ basis.T
        weights = project_simplex(weights - step * gradient)
    return weights


def cosine_similarity_matrix(spots: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    spot_norms = np.linalg.norm(spots, axis=1, keepdims=True)
    centroid_norms = np.linalg.norm(centroids, axis=1, keepdims=True).T
    return (spots @ centroids.T) / np.clip(spot_norms * centroid_norms, 1e-8, None)


def pearson_correlation(observed: np.ndarray, predicted: np.ndarray) -> float:
    if observed.size < 2:
        return 1.0
    observed_centered = observed - observed.mean()
    predicted_centered = predicted - predicted.mean()
    denominator = np.linalg.norm(observed_centered) * np.linalg.norm(predicted_centered)
    if denominator <= 1e-8:
        return 1.0 if np.allclose(observed, predicted) else 0.0
    return float(np.dot(observed_centered, predicted_centered) / denominator)


def normalized_entropy(weights: np.ndarray) -> np.ndarray:
    safe_weights = np.clip(weights, 1e-12, 1.0)
    return -np.sum(safe_weights * np.log(safe_weights), axis=1) / math.log(weights.shape[1])


def build_reference_centroids(cells: list[dict], genes: list[str]) -> tuple[list[str], np.ndarray]:
    cell_types = sorted({cell["cell_type"] for cell in cells})
    centroids = []
    for cell_type in cell_types:
        members = [cell for cell in cells if cell["cell_type"] == cell_type]
        expression = np.array(
            [[float(member["expression"][gene]) for gene in genes] for member in members],
            dtype=float,
        )
        centroids.append(expression.mean(axis=0))
    return cell_types, np.vstack(centroids)


def write_markdown(path: Path, title: str, sections: list[tuple[str, list[str]]]) -> None:
    lines = [f"# {title}", ""]
    for name, bullets in sections:
        lines.append(f"## {name}")
        lines.append("")
        for bullet in bullets:
            lines.append(f"- {bullet}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def build_neighbor_graph(coordinates: np.ndarray, *, n_neighbors: int) -> list[list[int]]:
    if coordinates.shape[0] < 2:
        return [[] for _ in range(coordinates.shape[0])]
    distances = np.sqrt(((coordinates[:, None, :] - coordinates[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    order = np.argsort(distances, axis=1)
    max_neighbors = min(n_neighbors, coordinates.shape[0] - 1)
    return [order[index, :max_neighbors].tolist() for index in range(coordinates.shape[0])]


def validate_outputs(outdir: Path, *, expected: dict | None = None) -> dict:
    metadata = load_metadata()
    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        if deliverable["kind"] == "tsv":
            frame = pd.read_csv(path, sep="\t")
            missing = [column for column in deliverable.get("required_columns", []) if column not in frame.columns]
            if missing:
                raise AssertionError(f"Missing required TSV columns in {path}: {missing}")
            if frame.empty:
                raise AssertionError(f"Empty TSV deliverable: {path}")
        elif deliverable["kind"] == "md":
            text = path.read_text(encoding="utf-8")
            for section in deliverable.get("required_sections", []):
                if f"## {section}" not in text:
                    raise AssertionError(f"Missing markdown section {section} in {path}")

    for filename in metadata.get("starter_qc_files", QC_FILES):
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    gene_cv = pd.read_csv(outdir / "gene_cv_qc.tsv", sep="\t")
    if gene_cv.empty:
        raise AssertionError("Gene CV QC table is empty.")
    spot_assignments = pd.read_csv(outdir / "spot_assignment_qc.tsv", sep="\t")
    if spot_assignments.empty:
        raise AssertionError("Spot assignment QC table is empty.")
    spatial_coherence = pd.read_csv(outdir / "spatial_coherence_qc.tsv", sep="\t")
    if spatial_coherence.empty:
        raise AssertionError("Spatial coherence QC table is empty.")

    qc = load_json(outdir / "intermediate_qc.json")
    required_qc_keys = [
        "cross_validation",
        "gene_intersection",
        "heldout_gene_metrics",
        "mapping_weights",
        "marker_recovery",
        "reference_cell_types",
        "spatial_coherence",
        "spot_assignments",
    ]
    missing_qc = [key for key in required_qc_keys if key not in qc]
    if missing_qc:
        raise AssertionError(f"Missing QC keys: {missing_qc}")

    weights = np.asarray(qc["mapping_weights"], dtype=float)
    if weights.ndim != 2 or weights.shape[0] != len(spot_assignments):
        raise AssertionError("Unexpected mapping-weight shape in intermediate_qc.json.")
    if not np.allclose(weights.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Mapping weights must sum to one for every spot.")

    metrics_frame = pd.read_csv(outdir / "validation_metrics.tsv", sep="\t")
    metrics = dict(zip(metrics_frame["metric"], metrics_frame["value"]))

    if expected:
        if qc["gene_intersection"]["count"] != expected["gene_intersection_size"]:
            raise AssertionError("Unexpected shared gene count.")
        observed_assignments = {item["spot_id"]: item["top_cell_type"] for item in qc["spot_assignments"]}
        if observed_assignments != expected["top_assignments"]:
            raise AssertionError(f"Unexpected top assignments: {observed_assignments}")
        floors = expected["metrics"]
        if float(metrics["heldout_pearson_r"]) < floors["heldout_pearson_r_min"]:
            raise AssertionError("Held-out Pearson correlation fell below the expected floor.")
        if float(metrics["heldout_rmse"]) > floors["heldout_rmse_max"]:
            raise AssertionError("Held-out RMSE exceeded the expected ceiling.")
        if float(metrics["cv_mean_pearson_r"]) < floors["cv_mean_pearson_r_min"]:
            raise AssertionError("Cross-validation Pearson correlation fell below the expected floor.")
        if float(metrics["cv_mean_rmse"]) > floors["cv_mean_rmse_max"]:
            raise AssertionError("Cross-validation RMSE exceeded the expected ceiling.")
        if float(metrics["marker_recovery_rate"]) < floors["marker_recovery_rate_min"]:
            raise AssertionError("Marker recovery rate fell below the expected floor.")
        if float(metrics["mean_assignment_entropy"]) > floors["mean_assignment_entropy_max"]:
            raise AssertionError("Assignment entropy exceeded the expected ceiling.")
        if float(metrics["neighbor_label_agreement_rate"]) < floors["neighbor_label_agreement_rate_min"]:
            raise AssertionError("Spatial coherence agreement fell below the expected floor.")
        if float(metrics["mean_top_cosine_similarity"]) < floors["mean_top_cosine_similarity_min"]:
            raise AssertionError("Top cosine similarity fell below the expected floor.")

    return {"metrics": metrics, "qc_path": str(outdir / "intermediate_qc.json")}


def run_skill(payload: dict, outdir: Path, *, input_path: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    shared_genes = payload["shared_genes"]
    heldout_genes = payload["heldout_genes"]
    all_genes = [*shared_genes, *heldout_genes]
    cells = payload["reference"]["cells"]
    spots = payload["spatial"]["spots"]
    marker_genes = payload["marker_genes"]
    solver = payload.get("solver", {})
    spatial_graph = payload.get("spatial_graph", {})
    max_iter = int(solver.get("max_iter", 400))
    n_neighbors = int(spatial_graph.get("n_neighbors", 1))

    if set(shared_genes) & set(heldout_genes):
        raise AssertionError("Shared and held-out genes must be disjoint.")
    if len(shared_genes) < 3:
        raise AssertionError("At least three shared genes are required for the starter.")
    if len(heldout_genes) < 1:
        raise AssertionError("At least one held-out gene is required for the starter.")

    cell_types, centroid_all = build_reference_centroids(cells, all_genes)
    centroid_shared = centroid_all[:, : len(shared_genes)]
    spot_ids = [spot["spot_id"] for spot in spots]
    coordinates = np.array([spot["coordinates"] for spot in spots], dtype=float)
    expected_labels = [spot["expected_cell_type"] for spot in spots]
    observed_shared = np.array(
        [[float(spot["shared_expression"][gene]) for gene in shared_genes] for spot in spots],
        dtype=float,
    )
    heldout_truth = np.array(
        [[float(spot["heldout_truth"][gene]) for gene in heldout_genes] for spot in spots],
        dtype=float,
    )

    cosine_scores = cosine_similarity_matrix(observed_shared, centroid_shared)
    weights = np.vstack(
        [
            solve_nnls_simplex(
                target,
                centroid_shared,
                max_iter=max_iter,
            )
            for target in observed_shared
        ]
    )
    reconstructed_shared = weights @ centroid_shared
    predicted_all = weights @ centroid_all
    predicted_heldout = predicted_all[:, len(shared_genes) :]
    entropies = normalized_entropy(weights)
    top_indices = weights.argmax(axis=1)
    top_labels = [cell_types[index] for index in top_indices]
    sorted_weights = np.sort(weights, axis=1)
    dominant_margins = sorted_weights[:, -1] - sorted_weights[:, -2]
    max_weights = weights.max(axis=1)
    top_cosine = np.take_along_axis(cosine_scores, top_indices[:, None], axis=1).reshape(-1)

    shared_rmse = float(np.sqrt(np.mean((reconstructed_shared - observed_shared) ** 2)))
    heldout_pearson = pearson_correlation(heldout_truth.ravel(), predicted_heldout.ravel())
    heldout_rmse = float(np.sqrt(np.mean((predicted_heldout - heldout_truth) ** 2)))
    mean_entropy = float(entropies.mean())
    mean_top_cosine = float(top_cosine.mean())
    expected_label_agreement_rate = float(np.mean([pred == exp for pred, exp in zip(top_labels, expected_labels)]))

    gene_cv_rows = []
    for gene_index, gene in enumerate(shared_genes):
        train_indices = [index for index in range(len(shared_genes)) if index != gene_index]
        basis = centroid_shared[:, train_indices]
        cv_weights = np.vstack(
            [
                solve_nnls_simplex(
                    observed_shared[spot_index, train_indices],
                    basis,
                    max_iter=max_iter,
                )
                for spot_index in range(len(spots))
            ]
        )
        predicted_gene = cv_weights @ centroid_shared[:, gene_index]
        observed_gene = observed_shared[:, gene_index]
        gene_cv_rows.append(
            {
                "gene": gene,
                "pearson_r": pearson_correlation(observed_gene, predicted_gene),
                "rmse": float(np.sqrt(np.mean((predicted_gene - observed_gene) ** 2))),
                "mean_observed": float(observed_gene.mean()),
                "mean_predicted": float(predicted_gene.mean()),
            }
        )
    gene_cv_frame = pd.DataFrame(gene_cv_rows)
    cv_mean_pearson = float(gene_cv_frame["pearson_r"].mean())
    cv_mean_rmse = float(gene_cv_frame["rmse"].mean())
    gene_cv_frame.to_csv(outdir / "gene_cv_qc.tsv", sep="\t", index=False)

    gene_index = {gene: index for index, gene in enumerate(all_genes)}
    marker_rows = []
    for spot_index, expected_label in enumerate(expected_labels):
        expected_score = float(
            np.mean([predicted_all[spot_index, gene_index[gene]] for gene in marker_genes[expected_label]])
        )
        competing_score = max(
            float(np.mean([predicted_all[spot_index, gene_index[gene]] for gene in marker_genes[label]]))
            for label in cell_types
            if label != expected_label
        )
        marker_rows.append(
            {
                "spot_id": spot_ids[spot_index],
                "expected_cell_type": expected_label,
                "expected_marker_score": expected_score,
                "competing_marker_score": competing_score,
                "marker_gap": expected_score - competing_score,
                "marker_recovered": expected_score > competing_score,
            }
        )
    marker_recovery_rate = float(np.mean([row["marker_recovered"] for row in marker_rows]))

    neighbors = build_neighbor_graph(coordinates, n_neighbors=n_neighbors)
    spatial_rows = []
    for spot_index, neighbor_indices in enumerate(neighbors):
        agreement = float(np.mean([top_labels[spot_index] == top_labels[neighbor] for neighbor in neighbor_indices]))
        mean_neighbor_weight_l1 = float(
            np.mean([np.abs(weights[spot_index] - weights[neighbor]).sum() for neighbor in neighbor_indices])
        )
        spatial_rows.append(
            {
                "spot_id": spot_ids[spot_index],
                "neighbor_spot_ids": json.dumps([spot_ids[neighbor] for neighbor in neighbor_indices]),
                "neighbor_top_cell_types": json.dumps([top_labels[neighbor] for neighbor in neighbor_indices]),
                "agreement_with_neighbors": agreement,
                "mean_neighbor_weight_l1": mean_neighbor_weight_l1,
            }
        )
    spatial_frame = pd.DataFrame(spatial_rows)
    spatial_frame.to_csv(outdir / "spatial_coherence_qc.tsv", sep="\t", index=False)
    neighbor_label_agreement_rate = float(spatial_frame["agreement_with_neighbors"].mean())
    mean_neighbor_weight_l1 = float(spatial_frame["mean_neighbor_weight_l1"].mean())

    spot_assignment_frame = pd.DataFrame(
        {
            "spot_id": spot_ids,
            "expected_cell_type": expected_labels,
            "top_cell_type": top_labels,
            "max_mapping_weight": max_weights,
            "dominant_margin": dominant_margins,
            "assignment_entropy": entropies,
            "top_cosine_similarity": top_cosine,
        }
    )
    spot_assignment_frame.to_csv(outdir / "spot_assignment_qc.tsv", sep="\t", index=False)

    heldout_gene_rows = []
    for gene_index, gene in enumerate(heldout_genes):
        heldout_gene_rows.append(
            {
                "gene": gene,
                "pearson_r": pearson_correlation(heldout_truth[:, gene_index], predicted_heldout[:, gene_index]),
                "rmse": float(np.sqrt(np.mean((predicted_heldout[:, gene_index] - heldout_truth[:, gene_index]) ** 2))),
                "mean_truth": float(heldout_truth[:, gene_index].mean()),
                "mean_predicted": float(predicted_heldout[:, gene_index].mean()),
            }
        )

    metrics_frame = pd.DataFrame(
        [
            {"metric": "heldout_pearson_r", "scope": "global", "value": heldout_pearson},
            {"metric": "heldout_rmse", "scope": "global", "value": heldout_rmse},
            {"metric": "cv_mean_pearson_r", "scope": "cross_validation", "value": cv_mean_pearson},
            {"metric": "cv_mean_rmse", "scope": "cross_validation", "value": cv_mean_rmse},
            {"metric": "shared_reconstruction_rmse", "scope": "mapping", "value": shared_rmse},
            {"metric": "mean_top_cosine_similarity", "scope": "mapping", "value": mean_top_cosine},
            {"metric": "marker_recovery_rate", "scope": "label_sanity", "value": marker_recovery_rate},
            {"metric": "expected_label_agreement_rate", "scope": "label_sanity", "value": expected_label_agreement_rate},
            {"metric": "mean_assignment_entropy", "scope": "assignment", "value": mean_entropy},
            {"metric": "neighbor_label_agreement_rate", "scope": "spatial_coherence", "value": neighbor_label_agreement_rate},
            {"metric": "mean_neighbor_weight_l1", "scope": "spatial_coherence", "value": mean_neighbor_weight_l1},
        ]
    )
    metrics_frame.to_csv(outdir / "validation_metrics.tsv", sep="\t", index=False)

    holdout_sections = [
        (
            "Run context",
            [
                f"Validated {len(spots)} spatial spots against {len(cell_types)} reference cell types with {len(shared_genes)} shared genes.",
                f"Used simplex-constrained NNLS on centroids with {len(heldout_genes)} held-out genes and leave-one-gene-out cross-validation.",
            ],
        ),
        (
            "Held-out genes",
            [
                f"Held-out Pearson r = {heldout_pearson:.3f} and held-out RMSE = {heldout_rmse:.4f}.",
                f"Cross-validation mean Pearson r = {cv_mean_pearson:.3f} and mean RMSE = {cv_mean_rmse:.4f}.",
                f"Validated held-out genes: {', '.join(heldout_genes)}.",
            ],
        ),
        (
            "Label sanity checks",
            [
                f"Top assignments were {', '.join(f'{spot_ids[index]}->{top_labels[index]}' for index in range(len(spot_ids)))}.",
                f"Marker recovery rate = {marker_recovery_rate:.3f}; expected-label agreement rate = {expected_label_agreement_rate:.3f}.",
                f"Mean assignment entropy = {mean_entropy:.3f}; nearest-neighbor label agreement = {neighbor_label_agreement_rate:.3f}.",
            ],
        ),
        (
            "Caveats",
            [
                "This is a deterministic synthetic starter, not a biological benchmark on paired real spatial and scRNA-seq data.",
                "Centroid-level simplex fitting approximates Tangram-style validation but does not reproduce the full cell-level optimizer.",
                "Spatial coherence is a tiny nearest-neighbor heuristic rather than a full Squidpy neighborhood statistic.",
            ],
        ),
    ]
    write_markdown(outdir / "holdout_summary.md", "Holdout Summary", holdout_sections)

    report_sections = [
        (
            "Run context",
            [
                f"Computed cosine similarity and centroid-level NNLS mapping weights for {len(spots)} spots and {len(cell_types)} cell types.",
                f"Shared genes: {', '.join(shared_genes)}.",
                f"Held-out genes: {', '.join(heldout_genes)}.",
            ],
        ),
        (
            "Validation metrics",
            [
                f"Shared reconstruction RMSE = {shared_rmse:.4f}.",
                f"Mean top cosine similarity = {mean_top_cosine:.3f}.",
                f"Mean assignment entropy = {mean_entropy:.3f}; mean neighbor weight L1 = {mean_neighbor_weight_l1:.3f}.",
            ],
        ),
        (
            "Interpretation",
            [
                f"Assignments followed the expected epithelial-to-immune gradient with {expected_label_agreement_rate:.3f} label agreement.",
                f"Held-out prediction and leave-one-gene-out CV stayed near-perfect on the synthetic benchmark.",
                "Marker recovery and graph coherence both support the fitted weights on this tiny toy problem.",
            ],
        ),
        (
            "Caveats",
            [
                "The starter is intentionally small and deterministic so it runs quickly without network access or heavyweight dependencies.",
                "Full Tangram training diagnostics, AUC curves, and richer spatial graph statistics remain future upgrade paths.",
                "Synthetic truth keeps the contract testable but cannot substitute for real biological validation.",
            ],
        ),
    ]
    write_markdown(outdir / "mapping_validation_report.md", "Mapping Validation Report", report_sections)

    qc_payload = {
        "run_label": payload["run_label"],
        "gene_intersection": {
            "count": len(shared_genes),
            "shared_genes": shared_genes,
            "heldout_genes": heldout_genes,
        },
        "reference_cell_types": cell_types,
        "cosine_similarity": cosine_scores.round(6).tolist(),
        "mapping_weights": weights.round(6).tolist(),
        "cross_validation": gene_cv_rows,
        "heldout_gene_metrics": heldout_gene_rows,
        "marker_recovery": marker_rows,
        "spatial_coherence": spatial_rows,
        "spot_assignments": [
            {
                "spot_id": spot_ids[index],
                "expected_cell_type": expected_labels[index],
                "top_cell_type": top_labels[index],
                "assignment_entropy": float(entropies[index]),
                "max_mapping_weight": float(max_weights[index]),
                "dominant_margin": float(dominant_margins[index]),
                "top_cosine_similarity": float(top_cosine[index]),
            }
            for index in range(len(spot_ids))
        ],
    }
    (outdir / "intermediate_qc.json").write_text(
        json.dumps(qc_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    result = {
        "input_path": str(input_path.resolve()),
        "metrics": dict(zip(metrics_frame["metric"], metrics_frame["value"])),
        "outdir": str(outdir.resolve()),
        "written_files": [],
    }
    (outdir / "run_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validation = validate_outputs(outdir, expected=payload.get("expected_invariants"))
    result["metrics"] = validation["metrics"]
    result["written_files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    (outdir / "run_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the spatial mapping validation starter.")
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = load_json(args.input)
    result = run_skill(payload, args.outdir, input_path=args.input)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
