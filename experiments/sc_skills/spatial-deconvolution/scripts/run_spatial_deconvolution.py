#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy.optimize import nnls


SKILL_DIR = Path(__file__).resolve().parents[1]
MODEL_KEYS = {
    "cell2location": "cell2location_like",
    "destvi": "destvi_like",
    "baseline": "uniform_baseline",
}
QC_FILES = [
    "qc_gene_intersection.json",
    "qc_neighbor_graph.tsv",
    "qc_spot_diagnostics.tsv",
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


def build_shared_gene_order(spatial_genes: list[str], reference_genes: list[str]) -> list[str]:
    reference_set = set(reference_genes)
    return [gene for gene in spatial_genes if gene in reference_set]


def normalize_rows(matrix: np.ndarray, target_sum: float) -> tuple[np.ndarray, np.ndarray]:
    row_sums = matrix.sum(axis=1)
    if np.any(row_sums <= 0):
        raise AssertionError("All profiles and spots must have positive total counts.")
    normalized = (matrix / row_sums[:, None]) * float(target_sum)
    return normalized, row_sums


def safe_pearson(observed: np.ndarray, predicted: np.ndarray) -> float:
    if np.allclose(observed, observed[0]) or np.allclose(predicted, predicted[0]):
        return 1.0 if np.allclose(observed, predicted) else 0.0
    return float(np.corrcoef(observed, predicted)[0, 1])


def shannon_entropy(probabilities: np.ndarray) -> float:
    clipped = np.clip(probabilities, 1e-12, None)
    return float(-(clipped * np.log2(clipped)).sum())


def compute_nnls_weights(observed_matrix: np.ndarray, reference_matrix: np.ndarray) -> np.ndarray:
    weight_rows: list[np.ndarray] = []
    for row_index, observed in enumerate(observed_matrix):
        weights, residual = nnls(reference_matrix.T, observed)
        if residual < 0:
            raise AssertionError(f"Unexpected negative residual for row index {row_index}.")
        if weights.sum() <= 0:
            raise AssertionError(f"NNLS returned zero weights for row index {row_index}.")
        weights = weights / weights.sum()
        weight_rows.append(weights.astype(float))
    return np.vstack(weight_rows)


def build_neighbor_graph(
    coordinates: np.ndarray,
    *,
    neighbor_count: int,
    self_weight: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    if len(coordinates) < 2:
        raise AssertionError("At least two spots are required to build the neighbor graph.")
    if not 0.0 <= self_weight < 1.0:
        raise AssertionError("self_weight must be in [0, 1).")

    capped_neighbor_count = max(1, min(neighbor_count, len(coordinates) - 1))
    distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)
    graph = np.zeros((len(coordinates), len(coordinates)), dtype=float)
    rows: list[dict[str, Any]] = []

    for source_index in range(len(coordinates)):
        order = np.argsort(distances[source_index])
        neighbors = [target for target in order if target != source_index][:capped_neighbor_count]
        inverse_distance = np.asarray(
            [1.0 / max(float(distances[source_index, target]), 1e-9) for target in neighbors],
            dtype=float,
        )
        inverse_distance = inverse_distance / inverse_distance.sum()

        graph[source_index, source_index] = self_weight
        rows.append(
            {
                "source_spot": source_index,
                "target_spot": source_index,
                "relation": "self",
                "distance": 0.0,
                "weight": round(float(self_weight), 6),
            }
        )
        for target, weight in zip(neighbors, inverse_distance, strict=True):
            edge_weight = float((1.0 - self_weight) * weight)
            graph[source_index, target] = edge_weight
            rows.append(
                {
                    "source_spot": source_index,
                    "target_spot": target,
                    "relation": "neighbor",
                    "distance": round(float(distances[source_index, target]), 6),
                    "weight": round(edge_weight, 6),
                }
            )

    return graph, pd.DataFrame(rows)


def build_neighbor_only_graph(graph: np.ndarray) -> np.ndarray:
    neighbor_only = graph.copy()
    np.fill_diagonal(neighbor_only, 0.0)
    row_sums = neighbor_only.sum(axis=1, keepdims=True)
    return np.divide(
        neighbor_only,
        row_sums,
        out=np.zeros_like(neighbor_only),
        where=row_sums > 0,
    )


def evaluate_method(
    *,
    method_key: str,
    weights: np.ndarray,
    reference_norm: np.ndarray,
    observed_norm: np.ndarray,
    fit_reference: np.ndarray,
    fit_observed: np.ndarray,
    truth_weights: np.ndarray,
    cell_types: list[str],
    spot_ids: list[str],
    neighbor_only_graph: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    predicted_norm = weights @ reference_norm
    predicted_fit = weights @ fit_reference
    rows: list[dict[str, Any]] = []

    for row_index, spot_id in enumerate(spot_ids):
        weight_row = weights[row_index]
        ranking = np.argsort(weight_row)[::-1]
        dominant_index = int(ranking[0])
        secondary_index = int(ranking[1])
        neighbor_average = neighbor_only_graph[row_index] @ weights
        rows.append(
            {
                "method": method_key,
                "spot_id": spot_id,
                "dominant_label": cell_types[dominant_index],
                "dominant_score": round(float(weight_row[dominant_index]), 6),
                "secondary_label": cell_types[secondary_index],
                "secondary_score": round(float(weight_row[secondary_index]), 6),
                "dominant_margin": round(
                    float(weight_row[dominant_index] - weight_row[secondary_index]),
                    6,
                ),
                "abundance_sum": round(float(weight_row.sum()), 6),
                "assignment_entropy": round(shannon_entropy(weight_row), 6),
                "fit_space_rmse": round(
                    float(np.sqrt(np.mean((predicted_fit[row_index] - fit_observed[row_index]) ** 2))),
                    6,
                ),
                "fit_space_correlation": round(
                    safe_pearson(fit_observed[row_index], predicted_fit[row_index]),
                    6,
                ),
                "profile_rmse": round(
                    float(np.sqrt(np.mean((predicted_norm[row_index] - observed_norm[row_index]) ** 2))),
                    6,
                ),
                "profile_correlation": round(
                    safe_pearson(observed_norm[row_index], predicted_norm[row_index]),
                    6,
                ),
                "abundance_rmse": round(
                    float(np.sqrt(np.mean((weight_row - truth_weights[row_index]) ** 2))),
                    6,
                ),
                "spatial_coherence_l1": round(
                    float(np.mean(np.abs(weight_row - neighbor_average))),
                    6,
                ),
            }
        )

    frame = pd.DataFrame(rows)
    summary = {
        "mean_fit_space_rmse": float(frame["fit_space_rmse"].mean()),
        "mean_fit_space_correlation": float(frame["fit_space_correlation"].mean()),
        "mean_profile_rmse": float(frame["profile_rmse"].mean()),
        "mean_profile_correlation": float(frame["profile_correlation"].mean()),
        "mean_abundance_rmse": float(frame["abundance_rmse"].mean()),
        "max_normalization_error": float(np.abs(frame["abundance_sum"] - 1.0).max()),
        "mean_spatial_coherence_l1": float(frame["spatial_coherence_l1"].mean()),
    }
    return rows, summary


def choose_selected_method(
    requested_method: str,
    method_summaries: dict[str, dict[str, float]],
) -> str:
    if requested_method != "best":
        return MODEL_KEYS[requested_method]
    ranked = sorted(
        method_summaries.items(),
        key=lambda item: (
            item[1]["mean_abundance_rmse"],
            item[1]["mean_profile_rmse"],
            -item[1]["mean_profile_correlation"],
        ),
    )
    return ranked[0][0]


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
        elif deliverable["kind"] == "h5ad":
            adata = ad.read_h5ad(path)
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise AssertionError(f"Empty AnnData artifact: {path}")
            if "spatial" not in adata.obsm:
                raise AssertionError(f"Missing spatial coordinates in {path}")
            abundance = np.asarray(adata.X, dtype=float)
            if not np.allclose(abundance.sum(axis=1), 1.0, atol=1e-6):
                raise AssertionError(f"Selected abundances must sum to 1 in {path}")
            for method_key in MODEL_KEYS.values():
                if method_key not in adata.layers:
                    raise AssertionError(f"Missing method layer {method_key} in {path}")
            if "selected_method" not in adata.uns:
                raise AssertionError(f"Missing selected_method in {path}")
        else:
            raise AssertionError(f"Unsupported deliverable kind: {deliverable['kind']}")

    for filename in metadata.get("starter_qc_files", QC_FILES):
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    gene_qc = load_json(outdir / "qc_gene_intersection.json")
    if gene_qc["shared_gene_count"] <= 0:
        raise AssertionError("The starter must keep at least one shared gene.")

    diagnostics = pd.read_csv(outdir / "qc_spot_diagnostics.tsv", sep="\t")
    if diagnostics.empty:
        raise AssertionError("Spot diagnostics QC table is empty.")
    if set(diagnostics["method"]) != set(MODEL_KEYS.values()):
        raise AssertionError("Spot diagnostics must include all starter methods.")
    if not np.allclose(diagnostics["abundance_sum"], 1.0, atol=1e-6):
        raise AssertionError("All method abundances must sum to 1 in the QC table.")

    neighbor_graph = pd.read_csv(outdir / "qc_neighbor_graph.tsv", sep="\t")
    if neighbor_graph.empty:
        raise AssertionError("Neighbor graph QC table is empty.")
    weight_sums = neighbor_graph.groupby("source_spot")["weight"].sum().to_numpy(dtype=float)
    if not np.allclose(weight_sums, 1.0, atol=1e-6):
        raise AssertionError("Neighbor graph rows must sum to 1.")

    comparison = pd.read_csv(outdir / "model_comparison.tsv", sep="\t")
    if not {"cell2location_like", "destvi_like", "uniform_baseline"}.issubset(set(comparison["model"])):
        raise AssertionError("Model comparison is missing one or more methods.")

    summary = load_json(outdir / "run_summary.json")
    if summary["selected_method"] not in MODEL_KEYS.values():
        raise AssertionError("Selected method must be one of the starter methods.")
    if summary["method_steps"] != [
        "gene_intersection",
        "cell2location_like_nnls",
        "destvi_like_smoothed_log_nnls",
        "uniform_baseline_comparison",
    ]:
        raise AssertionError("Unexpected method step summary.")


def build_deconvolution(
    payload: dict[str, Any],
    *,
    input_path: Path,
    outdir: Path,
    requested_method: str,
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)

    reference = payload["reference"]
    spatial = payload["spatial"]
    truth = payload["synthetic_truth"]
    expected = payload["expected_invariants"]

    reference_genes = list(reference["genes"])
    spatial_genes = list(spatial["genes"])
    shared_genes = build_shared_gene_order(spatial_genes, reference_genes)
    if not shared_genes:
        raise AssertionError("No shared genes between reference and spatial inputs.")

    reference_profiles = list(reference["profiles"])
    spots = list(spatial["spots"])
    if len(reference_profiles) < 2:
        raise AssertionError("At least two reference profiles are required.")
    if len(spots) < 2:
        raise AssertionError("At least two spatial spots are required.")

    cell_types = [profile["label"] for profile in reference_profiles]
    spot_ids = [spot["spot_id"] for spot in spots]
    sample_ids = [spot["sample_id"] for spot in spots]
    coordinates = np.asarray([[spot["x_coord"], spot["y_coord"]] for spot in spots], dtype=float)
    normalization_target_sum = float(payload.get("normalization_target_sum", 1000.0))
    neighbor_count = int(payload.get("neighbor_count", 2))
    self_weight = float(payload.get("self_weight", 0.95))

    reference_index = {gene: idx for idx, gene in enumerate(reference_genes)}
    spatial_index = {gene: idx for idx, gene in enumerate(spatial_genes)}
    reference_shared_idx = [reference_index[gene] for gene in shared_genes]
    spatial_shared_idx = [spatial_index[gene] for gene in shared_genes]

    reference_raw = np.asarray([profile["centroid"] for profile in reference_profiles], dtype=float)
    spatial_raw = np.asarray([spot["counts"] for spot in spots], dtype=float)
    if reference_raw.shape[1] != len(reference_genes):
        raise AssertionError("Reference centroid length does not match reference genes.")
    if spatial_raw.shape[1] != len(spatial_genes):
        raise AssertionError("Spot count length does not match spatial genes.")

    reference_shared_raw = reference_raw[:, reference_shared_idx]
    spatial_shared_raw = spatial_raw[:, spatial_shared_idx]
    reference_norm, reference_library_sizes = normalize_rows(reference_shared_raw, normalization_target_sum)
    spatial_norm, spatial_library_sizes = normalize_rows(spatial_shared_raw, normalization_target_sum)

    truth_types = list(truth["cell_types"])
    if truth_types != cell_types:
        raise AssertionError("synthetic_truth.cell_types must match reference profile labels.")
    truth_by_spot = {row["spot_id"]: np.asarray(row["weights"], dtype=float) for row in truth["abundances"]}
    truth_weights = np.vstack([truth_by_spot[spot_id] for spot_id in spot_ids])
    if not np.allclose(truth_weights.sum(axis=1), 1.0, atol=1e-6):
        raise AssertionError("Synthetic truth abundances must sum to 1 for each spot.")

    neighbor_graph, neighbor_frame = build_neighbor_graph(
        coordinates,
        neighbor_count=neighbor_count,
        self_weight=self_weight,
    )
    neighbor_only_graph = build_neighbor_only_graph(neighbor_graph)

    cell2location_weights = compute_nnls_weights(spatial_norm, reference_norm)
    smoothed_spatial_norm = neighbor_graph @ spatial_norm
    destvi_weights = compute_nnls_weights(np.log1p(smoothed_spatial_norm), np.log1p(reference_norm))
    baseline_weights = np.full(
        (len(spots), len(cell_types)),
        1.0 / float(len(cell_types)),
        dtype=float,
    )

    diagnostics_rows: list[dict[str, Any]] = []
    method_summaries: dict[str, dict[str, float]] = {}
    method_weights = {
        "cell2location_like": cell2location_weights,
        "destvi_like": destvi_weights,
        "uniform_baseline": baseline_weights,
    }
    method_fit_inputs = {
        "cell2location_like": (reference_norm, spatial_norm),
        "destvi_like": (np.log1p(reference_norm), np.log1p(smoothed_spatial_norm)),
        "uniform_baseline": (reference_norm, spatial_norm),
    }

    for method_key, weights in method_weights.items():
        fit_reference, fit_observed = method_fit_inputs[method_key]
        rows, summary = evaluate_method(
            method_key=method_key,
            weights=weights,
            reference_norm=reference_norm,
            observed_norm=spatial_norm,
            fit_reference=fit_reference,
            fit_observed=fit_observed,
            truth_weights=truth_weights,
            cell_types=cell_types,
            spot_ids=spot_ids,
            neighbor_only_graph=neighbor_only_graph,
        )
        diagnostics_rows.extend(rows)
        method_summaries[method_key] = summary

    selected_method = choose_selected_method(requested_method, method_summaries)
    selected_weights = method_weights[selected_method]
    selected_diagnostics = pd.DataFrame(diagnostics_rows)
    selected_rows = selected_diagnostics.loc[selected_diagnostics["method"] == selected_method].set_index("spot_id")

    obs_frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "selected_method": selected_method,
            "dominant_label": [selected_rows.loc[spot_id, "dominant_label"] for spot_id in spot_ids],
            "dominant_score": [selected_rows.loc[spot_id, "dominant_score"] for spot_id in spot_ids],
            "secondary_label": [selected_rows.loc[spot_id, "secondary_label"] for spot_id in spot_ids],
            "secondary_score": [selected_rows.loc[spot_id, "secondary_score"] for spot_id in spot_ids],
            "dominant_margin": [selected_rows.loc[spot_id, "dominant_margin"] for spot_id in spot_ids],
            "assignment_entropy": [selected_rows.loc[spot_id, "assignment_entropy"] for spot_id in spot_ids],
            "abundance_sum": [selected_rows.loc[spot_id, "abundance_sum"] for spot_id in spot_ids],
            "target_spot": [spot_id == expected["target_spot_id"] for spot_id in spot_ids],
            "true_dominant_label": [expected["dominant_labels"][spot_id] for spot_id in spot_ids],
        },
        index=spot_ids,
    )
    var_frame = pd.DataFrame({"cell_type": cell_types}, index=cell_types)

    abundance = ad.AnnData(X=selected_weights.astype(float), obs=obs_frame, var=var_frame)
    abundance.obsm["spatial"] = coordinates
    for method_key, weights in method_weights.items():
        abundance.layers[method_key] = weights.astype(float)
    abundance.uns["run_label"] = payload["run_label"]
    abundance.uns["selected_method"] = selected_method
    abundance.uns["available_methods"] = list(method_weights)
    abundance.uns["method_steps"] = [
        "gene_intersection",
        "cell2location_like_nnls",
        "destvi_like_smoothed_log_nnls",
        "uniform_baseline_comparison",
    ]
    abundance.uns["model_metrics"] = {
        key: {metric: round(float(value), 6) for metric, value in metrics.items()}
        for key, metrics in method_summaries.items()
    }
    abundance.uns["starter_scope"] = {
        "computed": [
            "gene intersection",
            "library-size normalization",
            "NNLS abundance estimation for a cell2location-like path",
            "toy spatial neighbor graph with smoothed log-NNLS for a DestVI-like path",
            "uniform-baseline comparison and per-spot QC metrics",
        ],
        "approximated": [
            "cell2location is approximated with deterministic NNLS instead of Bayesian posterior inference with tissue priors.",
            "DestVI is approximated with shallow neighbor smoothing and log-scale NNLS instead of latent variable training and gamma-state inference.",
        ],
        "outside_scope": [
            "uncertainty intervals",
            "reference signature learning from real single-cell data",
            "biological interpretation on real tissue",
        ],
    }
    abundance.write_h5ad(outdir / "cell_type_abundance.h5ad")

    comparison_rows: list[dict[str, Any]] = []
    ranked_methods = sorted(
        method_summaries,
        key=lambda key: (
            method_summaries[key]["mean_abundance_rmse"],
            method_summaries[key]["mean_profile_rmse"],
        ),
    )
    for rank, method_key in enumerate(ranked_methods, start=1):
        metrics = method_summaries[method_key]
        for metric_name in [
            "mean_abundance_rmse",
            "mean_fit_space_rmse",
            "mean_fit_space_correlation",
            "mean_profile_rmse",
            "mean_profile_correlation",
            "mean_spatial_coherence_l1",
            "max_normalization_error",
        ]:
            comparison_rows.append(
                {
                    "model": method_key,
                    "metric": metric_name,
                    "value": round(float(metrics[metric_name]), 6),
                }
            )
        comparison_rows.append(
            {
                "model": method_key,
                "metric": "selected_output",
                "value": 1.0 if method_key == selected_method else 0.0,
            }
        )
        comparison_rows.append(
            {
                "model": method_key,
                "metric": "rank_by_abundance_rmse",
                "value": float(rank),
            }
        )
    comparison_frame = pd.DataFrame(comparison_rows)
    comparison_frame.to_csv(outdir / "model_comparison.tsv", sep="\t", index=False)

    diagnostics_frame = pd.DataFrame(diagnostics_rows)
    diagnostics_frame.to_csv(outdir / "qc_spot_diagnostics.tsv", sep="\t", index=False)

    neighbor_frame = neighbor_frame.copy()
    neighbor_frame["source_spot"] = [spot_ids[index] for index in neighbor_frame["source_spot"]]
    neighbor_frame["target_spot"] = [spot_ids[index] for index in neighbor_frame["target_spot"]]
    neighbor_frame.to_csv(outdir / "qc_neighbor_graph.tsv", sep="\t", index=False)

    gene_qc = {
        "run_label": payload["run_label"],
        "shared_gene_count": len(shared_genes),
        "shared_genes": shared_genes,
        "reference_gene_count": len(reference_genes),
        "spatial_gene_count": len(spatial_genes),
        "cell_types": cell_types,
        "spot_ids": spot_ids,
        "normalization_target_sum": normalization_target_sum,
        "neighbor_count": min(max(1, neighbor_count), len(spots) - 1),
        "self_weight": self_weight,
    }
    (outdir / "qc_gene_intersection.json").write_text(
        json.dumps(gene_qc, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    best_metrics = method_summaries[selected_method]
    target_row = selected_rows.loc[expected["target_spot_id"]]
    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}.",
                f"Processed {len(spots)} synthetic spots, {len(cell_types)} cell types, and {len(shared_genes)} shared genes.",
                f"Selected output method: {selected_method} using the smallest synthetic-truth abundance RMSE.",
            ],
        },
        {
            "name": "Abundance matrix",
            "bullets": [
                f"`cell_type_abundance.h5ad` stores the selected abundance matrix in `X` and all starter methods in AnnData layers.",
                f"Per-spot abundance sums stayed within {best_metrics['max_normalization_error']:.2e} of 1.0.",
                f"Target spot {expected['target_spot_id']} remained {target_row['dominant_label']} with score {target_row['dominant_score']:.3f}.",
            ],
        },
        {
            "name": "Model comparison",
            "bullets": [
                f"cell2location-like mean abundance RMSE: {method_summaries['cell2location_like']['mean_abundance_rmse']:.3f}.",
                f"DestVI-like mean abundance RMSE: {method_summaries['destvi_like']['mean_abundance_rmse']:.3f}.",
                f"Uniform baseline mean abundance RMSE: {method_summaries['uniform_baseline']['mean_abundance_rmse']:.3f}.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter computes deterministic NNLS surrogates, not full cell2location or DestVI posterior training.",
                "The score-based model choice uses synthetic truth available only because the toy example is generated data.",
                "Real datasets still require atlas selection, quality control, and method-specific hyperparameter tuning outside this starter.",
            ],
        },
    ]
    write_markdown(outdir / "deconvolution_summary.md", "Deconvolution Summary", report_sections)

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path.resolve()),
        "outdir": str(outdir.resolve()),
        "selected_method": selected_method,
        "available_methods": list(method_weights),
        "method_steps": [
            "gene_intersection",
            "cell2location_like_nnls",
            "destvi_like_smoothed_log_nnls",
            "uniform_baseline_comparison",
        ],
        "spot_count": len(spots),
        "cell_type_count": len(cell_types),
        "shared_gene_count": len(shared_genes),
        "reference_library_sizes": [float(value) for value in reference_library_sizes.tolist()],
        "spatial_library_sizes": [float(value) for value in spatial_library_sizes.tolist()],
        "target_spot_id": expected["target_spot_id"],
        "target_spot_dominant_label": str(target_row["dominant_label"]),
        "model_metrics": {
            key: {metric: round(float(value), 6) for metric, value in metrics.items()}
            for key, metrics in method_summaries.items()
        },
    }
    (outdir / "run_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--method",
        choices=["best", "cell2location", "destvi", "baseline"],
        default="best",
    )
    args = parser.parse_args()
    payload = load_json(args.input)
    build_deconvolution(
        payload,
        input_path=args.input,
        outdir=args.outdir,
        requested_method=args.method,
    )
    validate_outputs(SKILL_DIR, args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
