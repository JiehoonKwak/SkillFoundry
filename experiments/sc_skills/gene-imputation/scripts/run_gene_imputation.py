#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SKILL_DIR = Path(__file__).resolve().parents[1]


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


def write_markdown_report(path: Path, sections: list[tuple[str, list[str]]]) -> None:
    lines = ["# Gene Imputation Report", ""]
    for name, bullets in sections:
        lines.append(f"## {name}")
        lines.append("")
        for bullet in bullets:
            lines.append(f"- {bullet}")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


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
        elif deliverable["kind"] == "md":
            text = path.read_text(encoding="utf-8")
            for section in deliverable.get("required_sections", []):
                if f"## {section}" not in text:
                    raise AssertionError(f"Missing markdown section {section} in {path}")
        elif deliverable["kind"] == "h5ad":
            adata = ad.read_h5ad(path)
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise AssertionError(f"Empty AnnData output: {path}")
            if "mapping_weights" not in adata.obsm:
                raise AssertionError("Expected `mapping_weights` in `imputed_expression.h5ad`.")
            weights = np.asarray(adata.obsm["mapping_weights"], dtype=float)
            if not np.allclose(weights.sum(axis=1), 1.0, atol=1e-6):
                raise AssertionError("Mapping weights must sum to one for every spot.")

    qc_path = outdir / "intermediate_qc.json"
    if not qc_path.exists():
        raise AssertionError(f"Missing intermediate QC file: {qc_path}")
    qc = load_json(qc_path)
    required_qc_keys = [
        "cosine_similarity",
        "mapping_weights",
        "reference_cell_types",
        "shared_reconstruction_rmse",
        "spot_assignments",
    ]
    missing_qc = [key for key in required_qc_keys if key not in qc]
    if missing_qc:
        raise AssertionError(f"Missing QC keys: {missing_qc}")

    metrics_frame = pd.read_csv(outdir / "heldout_metrics.tsv", sep="\t")
    metrics = dict(zip(metrics_frame["metric"], metrics_frame["value"]))

    if expected:
        if qc["gene_intersection"]["count"] != expected["gene_intersection_size"]:
            raise AssertionError("Unexpected shared gene count.")
        expected_assignments = expected["top_assignments"]
        observed_assignments = {item["spot_id"]: item["top_cell_type"] for item in qc["spot_assignments"]}
        if observed_assignments != expected_assignments:
            raise AssertionError(f"Unexpected top assignments: {observed_assignments}")
        floors = expected["metrics"]
        if float(metrics["heldout_pearson_r"]) < floors["heldout_pearson_r_min"]:
            raise AssertionError("Held-out Pearson correlation fell below the expected floor.")
        if float(metrics["heldout_rmse"]) > floors["heldout_rmse_max"]:
            raise AssertionError("Held-out RMSE exceeded the expected ceiling.")
        if float(metrics["marker_recovery_rate"]) < floors["marker_recovery_rate_min"]:
            raise AssertionError("Marker recovery rate fell below the expected floor.")
        if float(metrics["mean_assignment_entropy"]) > floors["mean_assignment_entropy_max"]:
            raise AssertionError("Assignment entropy exceeded the expected ceiling.")

    return {"metrics": metrics, "qc_path": str(qc_path)}


def run_skill(payload: dict, outdir: Path, *, input_path: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)

    shared_genes = payload["shared_genes"]
    heldout_genes = payload["heldout_genes"]
    all_genes = [*shared_genes, *heldout_genes]
    cells = payload["reference"]["cells"]
    spots = payload["spatial"]["spots"]
    marker_genes = payload["marker_genes"]
    solver = payload.get("solver", {})

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
                max_iter=int(solver.get("max_iter", 400)),
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
    max_weights = weights.max(axis=1)

    gene_index = {gene: index for index, gene in enumerate(all_genes)}
    marker_status = []
    for spot_index, expected_label in enumerate(expected_labels):
        expected_score = float(np.mean([predicted_all[spot_index, gene_index[gene]] for gene in marker_genes[expected_label]]))
        competing_labels = [label for label in cell_types if label != expected_label]
        competing_score = max(
            float(np.mean([predicted_all[spot_index, gene_index[gene]] for gene in marker_genes[label]]))
            for label in competing_labels
        )
        marker_status.append(
            {
                "spot_id": spot_ids[spot_index],
                "expected_cell_type": expected_label,
                "expected_marker_score": expected_score,
                "competing_marker_score": competing_score,
                "marker_gap": expected_score - competing_score,
                "marker_recovered": expected_score > competing_score,
            }
        )

    heldout_pearson = pearson_correlation(heldout_truth.ravel(), predicted_heldout.ravel())
    heldout_rmse = float(np.sqrt(np.mean((predicted_heldout - heldout_truth) ** 2)))
    shared_rmse = float(np.sqrt(np.mean((reconstructed_shared - observed_shared) ** 2)))
    marker_recovery_rate = float(np.mean([item["marker_recovered"] for item in marker_status]))
    mean_entropy = float(entropies.mean())
    mean_top_cosine = float(np.mean(np.take_along_axis(cosine_scores, top_indices[:, None], axis=1)))

    metrics_frame = pd.DataFrame(
        [
            {"metric": "heldout_pearson_r", "value": heldout_pearson, "higher_is_better": True},
            {"metric": "heldout_rmse", "value": heldout_rmse, "higher_is_better": False},
            {"metric": "shared_reconstruction_rmse", "value": shared_rmse, "higher_is_better": False},
            {"metric": "marker_recovery_rate", "value": marker_recovery_rate, "higher_is_better": True},
            {"metric": "mean_assignment_entropy", "value": mean_entropy, "higher_is_better": False},
            {"metric": "mean_top_cosine_similarity", "value": mean_top_cosine, "higher_is_better": True},
        ]
    )
    metrics_path = outdir / "heldout_metrics.tsv"
    metrics_frame.to_csv(metrics_path, sep="\t", index=False)

    obs_frame = pd.DataFrame(
        {
            "expected_cell_type": expected_labels,
            "top_cell_type": top_labels,
            "assignment_entropy": entropies,
            "max_mapping_weight": max_weights,
            "x_coord": coordinates[:, 0],
            "y_coord": coordinates[:, 1],
        },
        index=spot_ids,
    )
    var_frame = pd.DataFrame(
        {
            "gene_role": ["shared"] * len(shared_genes) + ["heldout"] * len(heldout_genes),
            "is_heldout": [False] * len(shared_genes) + [True] * len(heldout_genes),
        },
        index=all_genes,
    )
    adata = ad.AnnData(X=predicted_all, obs=obs_frame, var=var_frame)
    adata.obsm["spatial"] = coordinates
    adata.obsm["mapping_weights"] = weights
    adata.uns["mapping_cell_types"] = np.asarray(cell_types, dtype="U")
    adata.uns["cosine_similarity"] = cosine_scores
    adata.uns["starter_scope"] = {
        "computed": "Reference centroids, simplex-constrained NNLS mapping weights, projected genes, held-out metrics.",
        "approximated": "Cell-level Tangram optimization is approximated by centroid-level deterministic weights.",
        "surrogate": "Synthetic held-out truth replaces a real paired benchmark slice.",
    }
    h5ad_path = outdir / "imputed_expression.h5ad"
    adata.write_h5ad(h5ad_path)

    report_sections = [
        (
            "Run context",
            [
                f"Computed centroid-level mapping for {len(spots)} spatial spots against {len(cell_types)} reference cell types.",
                f"Used {len(shared_genes)} shared genes for fitting and {len(heldout_genes)} held-out genes for validation.",
                "The starter replaces full Tangram optimization with simplex-constrained NNLS on reference centroids.",
            ],
        ),
        (
            "Held-out validation",
            [
                f"Held-out Pearson r = {heldout_pearson:.3f}.",
                f"Held-out RMSE = {heldout_rmse:.4f}.",
                f"Shared-gene reconstruction RMSE = {shared_rmse:.4f}.",
            ],
        ),
        (
            "Interpretation",
            [
                f"Top assignments were {', '.join(f'{spot_ids[i]}->{top_labels[i]}' for i in range(len(spot_ids)))}.",
                f"Marker recovery rate = {marker_recovery_rate:.3f}; lower entropy indicates cleaner spot-to-centroid assignments.",
                "Projected held-out genes preserve the epithelial versus immune contrast encoded in the toy reference.",
            ],
        ),
        (
            "Caveats",
            [
                "This is a deterministic synthetic starter, not a biological benchmark on real paired scRNA-seq and spatial data.",
                "Weights are estimated from cell-type centroids rather than a full Tangram cell-to-spot optimization.",
                "Held-out validation uses synthetic truth carried in the toy input so the contract stays fast and portable.",
            ],
        ),
    ]
    report_path = outdir / "imputation_report.md"
    write_markdown_report(report_path, report_sections)

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
        "shared_reconstruction_rmse": shared_rmse,
        "spot_assignments": [
            {
                "spot_id": spot_ids[index],
                "expected_cell_type": expected_labels[index],
                "top_cell_type": top_labels[index],
                "assignment_entropy": float(entropies[index]),
                "max_mapping_weight": float(max_weights[index]),
            }
            for index in range(len(spot_ids))
        ],
        "marker_recovery": marker_status,
        "heldout_gene_metrics": [
            {
                "gene": gene,
                "pearson_r": pearson_correlation(heldout_truth[:, index], predicted_heldout[:, index]),
                "rmse": float(np.sqrt(np.mean((predicted_heldout[:, index] - heldout_truth[:, index]) ** 2))),
            }
            for index, gene in enumerate(heldout_genes)
        ],
    }
    qc_path = outdir / "intermediate_qc.json"
    qc_path.write_text(json.dumps(qc_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    validation = validate_outputs(outdir, expected=payload.get("expected_invariants"))
    result = {
        "input_path": str(input_path.resolve()),
        "metrics": validation["metrics"],
        "outdir": str(outdir.resolve()),
        "written_files": [],
    }
    (outdir / "run_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result["written_files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    (outdir / "run_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the spatial gene imputation starter.")
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = load_json(args.input)
    result = run_skill(payload, args.outdir, input_path=args.input)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
