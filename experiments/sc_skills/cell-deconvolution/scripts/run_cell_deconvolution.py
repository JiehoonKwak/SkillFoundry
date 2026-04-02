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
ASSIGNMENT_COLUMNS = [
    "spot_id",
    "dominant_label",
    "dominant_score",
    "secondary_label",
    "secondary_score",
]
QC_FILES = [
    "qc_gene_intersection.json",
    "qc_spot_metrics.tsv",
    "qc_segment_projection.tsv",
    "qc_holdout_predictions.tsv",
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


def normalize_log1p(matrix: np.ndarray, target_sum: float) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = matrix.sum(axis=1)
    if np.any(library_sizes <= 0):
        raise AssertionError("All profiles and spots must have positive total counts.")
    normalized = (matrix / library_sizes[:, None]) * float(target_sum)
    return np.log1p(normalized), library_sizes


def shannon_entropy(probabilities: np.ndarray) -> float:
    clipped = np.clip(probabilities, 1e-12, None)
    return float(-(clipped * np.log2(clipped)).sum())


def safe_pearson(observed: np.ndarray, predicted: np.ndarray) -> float:
    if np.allclose(observed, observed[0]) or np.allclose(predicted, predicted[0]):
        return 1.0 if np.allclose(observed, predicted) else 0.0
    return float(np.corrcoef(observed, predicted)[0, 1])


def integerize_expected_counts(weights: np.ndarray, total_count: int) -> tuple[np.ndarray, np.ndarray]:
    expected = weights * float(total_count)
    assigned = np.floor(expected).astype(int)
    remaining = int(total_count - int(assigned.sum()))
    remainders = expected - assigned
    order = sorted(range(len(weights)), key=lambda idx: (-float(remainders[idx]), idx))
    for idx in order[:remaining]:
        assigned[idx] += 1
    return expected, assigned


def compute_nnls_weights(spot_matrix: np.ndarray, reference_matrix: np.ndarray) -> np.ndarray:
    weight_rows: list[np.ndarray] = []
    for spot_index, spot_profile in enumerate(spot_matrix):
        weights, residual = nnls(reference_matrix.T, spot_profile)
        if residual < 0:
            raise AssertionError(f"Unexpected negative residual for spot index {spot_index}.")
        if weights.sum() <= 0:
            raise AssertionError(f"NNLS returned zero weights for spot index {spot_index}.")
        weights = weights / weights.sum()
        weight_rows.append(weights.astype(float))
    return np.vstack(weight_rows)


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
            weights = np.asarray(adata.X, dtype=float)
            if not np.allclose(weights.sum(axis=1), 1.0, atol=1e-6):
                raise AssertionError(f"Mapping weights must sum to 1 for each spot in {path}")
            required_obs = [
                "dominant_label",
                "secondary_label",
                "dominant_margin",
                "assignment_entropy",
                "segmentation_cell_count",
            ]
            missing_obs = [column for column in required_obs if column not in adata.obs.columns]
            if missing_obs:
                raise AssertionError(f"Missing mapping obs columns in {path}: {missing_obs}")
        else:
            raise AssertionError(f"Unsupported deliverable kind: {deliverable['kind']}")

    for filename in metadata.get("starter_qc_files", QC_FILES):
        qc_path = outdir / filename
        if not qc_path.exists():
            raise AssertionError(f"Missing QC artifact: {qc_path}")

    gene_qc = load_json(outdir / "qc_gene_intersection.json")
    if gene_qc["fit_gene_count"] <= 0 or gene_qc["holdout_gene_count"] <= 0:
        raise AssertionError("The starter must keep both fit and holdout genes.")

    spot_metrics = pd.read_csv(outdir / "qc_spot_metrics.tsv", sep="\t")
    if spot_metrics.empty:
        raise AssertionError("Spot metrics QC table is empty.")

    segment_projection = pd.read_csv(outdir / "qc_segment_projection.tsv", sep="\t")
    if segment_projection.empty:
        raise AssertionError("Segment projection QC table is empty.")
    projected_by_spot = segment_projection.groupby("spot_id")["assigned_segment_count"].sum().to_dict()
    expected_by_spot = spot_metrics.set_index("spot_id")["segmentation_cell_count"].to_dict()
    if projected_by_spot != expected_by_spot:
        raise AssertionError("Projected segment counts do not match segmentation cell counts.")

    holdout_predictions = pd.read_csv(outdir / "qc_holdout_predictions.tsv", sep="\t")
    required_holdout_columns = {
        "spot_id",
        "gene",
        "observed_value",
        "predicted_value",
        "absolute_error",
    }
    if not required_holdout_columns.issubset(holdout_predictions.columns):
        raise AssertionError("Holdout prediction QC table is missing required columns.")
    if holdout_predictions.empty:
        raise AssertionError("Holdout prediction QC table is empty.")

    assignments = pd.read_csv(outdir / "cell_assignment.tsv", sep="\t").set_index("spot_id")
    mapping = ad.read_h5ad(outdir / "mapping_matrix.h5ad")
    mapping_frame = pd.DataFrame(np.asarray(mapping.X), index=mapping.obs_names, columns=mapping.var_names)
    for spot_id, row in assignments.iterrows():
        ranked = mapping_frame.loc[spot_id].sort_values(ascending=False)
        if row["dominant_label"] != ranked.index[0]:
            raise AssertionError(f"Dominant label mismatch for {spot_id}")
        if row["secondary_label"] != ranked.index[1]:
            raise AssertionError(f"Secondary label mismatch for {spot_id}")


def build_deconvolution(payload: dict[str, Any], *, input_path: Path, outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)

    reference = payload["reference"]
    spatial = payload["spatial"]
    marker_table = payload["marker_table"]
    fit_genes = list(payload["fit_genes"])
    holdout_genes = list(payload["holdout_genes"])
    target_sum = float(payload.get("normalization_target_sum", 1000.0))

    reference_genes = list(reference["genes"])
    spatial_genes = list(spatial["genes"])
    shared_genes = build_shared_gene_order(spatial_genes, reference_genes)
    if not shared_genes:
        raise AssertionError("No shared genes between reference and spatial inputs.")
    if set(fit_genes) & set(holdout_genes):
        raise AssertionError("Fit and holdout genes must be disjoint.")
    if not set(fit_genes).issubset(shared_genes):
        raise AssertionError("Fit genes must be shared between reference and spatial inputs.")
    if not set(holdout_genes).issubset(shared_genes):
        raise AssertionError("Holdout genes must be shared between reference and spatial inputs.")
    if len(fit_genes) < 2:
        raise AssertionError("At least two fit genes are required.")
    if len(holdout_genes) < 1:
        raise AssertionError("At least one holdout gene is required.")

    profiles = list(reference["profiles"])
    spots = list(spatial["spots"])
    if len(profiles) < 2:
        raise AssertionError("At least two reference profiles are required.")
    if len(spots) < 2:
        raise AssertionError("At least two spatial spots are required.")

    reference_labels = [profile["label"] for profile in profiles]
    broad_labels = [profile.get("broad_label", profile["label"]) for profile in profiles]
    marker_table_by_label = {
        entry["label"]: [gene for gene in entry["markers"] if gene in fit_genes]
        for entry in marker_table
    }
    missing_marker_labels = [label for label in reference_labels if not marker_table_by_label.get(label)]
    if missing_marker_labels:
        raise AssertionError(f"Missing fit-gene markers for labels: {missing_marker_labels}")

    reference_index = {gene: idx for idx, gene in enumerate(reference_genes)}
    spatial_index = {gene: idx for idx, gene in enumerate(spatial_genes)}
    shared_reference_idx = [reference_index[gene] for gene in shared_genes]
    shared_spatial_idx = [spatial_index[gene] for gene in shared_genes]
    shared_gene_to_idx = {gene: idx for idx, gene in enumerate(shared_genes)}

    reference_raw = np.asarray([profile["centroid"] for profile in profiles], dtype=float)
    spatial_raw = np.asarray([spot["counts"] for spot in spots], dtype=float)
    if reference_raw.shape[1] != len(reference_genes):
        raise AssertionError("Reference centroid length does not match reference genes.")
    if spatial_raw.shape[1] != len(spatial_genes):
        raise AssertionError("Spatial spot length does not match spatial genes.")

    reference_shared_raw = reference_raw[:, shared_reference_idx]
    spatial_shared_raw = spatial_raw[:, shared_spatial_idx]
    reference_shared_norm, reference_library_sizes = normalize_log1p(reference_shared_raw, target_sum)
    spatial_shared_norm, spatial_library_sizes = normalize_log1p(spatial_shared_raw, target_sum)

    fit_indices = [shared_gene_to_idx[gene] for gene in fit_genes]
    holdout_indices = [shared_gene_to_idx[gene] for gene in holdout_genes]
    reference_fit = reference_shared_norm[:, fit_indices]
    spatial_fit = spatial_shared_norm[:, fit_indices]
    reference_holdout = reference_shared_norm[:, holdout_indices]
    spatial_holdout = spatial_shared_norm[:, holdout_indices]

    weights = compute_nnls_weights(spatial_fit, reference_fit)
    predicted_fit = weights @ reference_fit
    predicted_holdout = weights @ reference_holdout

    spot_ids = [spot["spot_id"] for spot in spots]
    sample_ids = [spot["sample_id"] for spot in spots]
    segmentation_counts = [int(spot["segmentation_cell_count"]) for spot in spots]
    coordinates = np.asarray([[spot["x_coord"], spot["y_coord"]] for spot in spots], dtype=float)

    holdout_rows: list[dict[str, Any]] = []
    segment_rows: list[dict[str, Any]] = []
    assignment_rows: list[dict[str, Any]] = []
    spot_metric_rows: list[dict[str, Any]] = []

    for row_index, spot_id in enumerate(spot_ids):
        weight_row = weights[row_index]
        ranking = np.argsort(weight_row)[::-1]
        dominant_idx = int(ranking[0])
        secondary_idx = int(ranking[1])
        dominant_label = reference_labels[dominant_idx]
        secondary_label = reference_labels[secondary_idx]
        dominant_score = float(weight_row[dominant_idx])
        secondary_score = float(weight_row[secondary_idx])
        dominant_margin = dominant_score - secondary_score
        entropy = shannon_entropy(weight_row)
        fit_reconstruction_mae = float(np.abs(predicted_fit[row_index] - spatial_fit[row_index]).mean())

        dominant_marker_values = [
            spatial_raw[row_index, spatial_index[gene]]
            for gene in marker_table_by_label[dominant_label]
        ]
        secondary_marker_values = [
            spatial_raw[row_index, spatial_index[gene]]
            for gene in marker_table_by_label[secondary_label]
        ]
        dominant_marker_mean = float(np.mean(dominant_marker_values))
        secondary_marker_mean = float(np.mean(secondary_marker_values))
        dominant_marker_recovery = dominant_marker_mean / (
            dominant_marker_mean + secondary_marker_mean + 1e-12
        )

        expected_segment_counts, assigned_segment_counts = integerize_expected_counts(
            weight_row,
            segmentation_counts[row_index],
        )
        segment_slot = 1
        for label_idx, label in enumerate(reference_labels):
            segment_ids = []
            for _ in range(int(assigned_segment_counts[label_idx])):
                segment_ids.append(f"{spot_id}_cell_{segment_slot:02d}")
                segment_slot += 1
            segment_rows.append(
                {
                    "spot_id": spot_id,
                    "label": label,
                    "weight": round(float(weight_row[label_idx]), 6),
                    "expected_segment_count": round(float(expected_segment_counts[label_idx]), 6),
                    "assigned_segment_count": int(assigned_segment_counts[label_idx]),
                    "assigned_segment_ids": ";".join(segment_ids),
                }
            )

        assignment_rows.append(
            {
                "spot_id": spot_id,
                "dominant_label": dominant_label,
                "dominant_score": round(dominant_score, 6),
                "secondary_label": secondary_label,
                "secondary_score": round(secondary_score, 6),
            }
        )
        spot_metric_rows.append(
            {
                "spot_id": spot_id,
                "segmentation_cell_count": segmentation_counts[row_index],
                "dominant_label": dominant_label,
                "dominant_score": round(dominant_score, 6),
                "secondary_label": secondary_label,
                "secondary_score": round(secondary_score, 6),
                "dominant_margin": round(dominant_margin, 6),
                "assignment_entropy": round(entropy, 6),
                "fit_reconstruction_mae": round(fit_reconstruction_mae, 6),
                "dominant_marker_recovery": round(dominant_marker_recovery, 6),
            }
        )

        for gene_index, gene in enumerate(holdout_genes):
            observed_value = float(spatial_holdout[row_index, gene_index])
            predicted_value = float(predicted_holdout[row_index, gene_index])
            holdout_rows.append(
                {
                    "spot_id": spot_id,
                    "gene": gene,
                    "observed_value": round(observed_value, 6),
                    "predicted_value": round(predicted_value, 6),
                    "absolute_error": round(abs(observed_value - predicted_value), 6),
                }
            )

    holdout_frame = pd.DataFrame(holdout_rows)
    holdout_gene_metrics: dict[str, dict[str, float]] = {}
    for gene in holdout_genes:
        subset = holdout_frame.loc[holdout_frame["gene"] == gene]
        observed = subset["observed_value"].to_numpy(dtype=float)
        predicted = subset["predicted_value"].to_numpy(dtype=float)
        holdout_gene_metrics[gene] = {
            "pearson_r": round(safe_pearson(observed, predicted), 6),
            "mae": round(float(np.abs(observed - predicted).mean()), 6),
            "rmse": round(float(np.sqrt(np.mean((observed - predicted) ** 2))), 6),
        }

    assignment_frame = pd.DataFrame(assignment_rows, columns=ASSIGNMENT_COLUMNS)
    spot_metric_frame = pd.DataFrame(spot_metric_rows)
    segment_frame = pd.DataFrame(segment_rows)

    obs_frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "segmentation_cell_count": segmentation_counts,
            "dominant_label": spot_metric_frame["dominant_label"].tolist(),
            "secondary_label": spot_metric_frame["secondary_label"].tolist(),
            "dominant_score": spot_metric_frame["dominant_score"].tolist(),
            "secondary_score": spot_metric_frame["secondary_score"].tolist(),
            "dominant_margin": spot_metric_frame["dominant_margin"].tolist(),
            "assignment_entropy": spot_metric_frame["assignment_entropy"].tolist(),
            "fit_reconstruction_mae": spot_metric_frame["fit_reconstruction_mae"].tolist(),
            "dominant_marker_recovery": spot_metric_frame["dominant_marker_recovery"].tolist(),
        },
        index=spot_ids,
    )
    var_frame = pd.DataFrame(
        {
            "cell_type": reference_labels,
            "broad_label": broad_labels,
        },
        index=reference_labels,
    )
    mapping = ad.AnnData(X=weights.astype(float), obs=obs_frame, var=var_frame)
    mapping.obsm["spatial"] = coordinates
    mapping.uns["run_label"] = payload["run_label"]
    mapping.uns["method_steps"] = [
        "gene_intersection",
        "nnls_mapping",
        "segmentation_projection",
        "holdout_gene_validation",
    ]
    mapping.uns["fit_genes"] = fit_genes
    mapping.uns["holdout_genes"] = holdout_genes
    mapping.uns["holdout_gene_metrics"] = holdout_gene_metrics
    mapping.uns["starter_scope"] = {
        "computed": [
            "gene intersection",
            "library-size normalization with log1p",
            "NNLS mapping weights from reference centroids to spots",
            "integerized segment-count projection from per-spot weights",
            "held-out gene prediction and entropy-based QC",
        ],
        "approximated": [
            "Tangram constrained mapping is approximated with centroid-level NNLS rather than GPU-trained cell-to-space optimization.",
            "Segmentation-aware assignment uses provided toy spot cell counts instead of image-derived segmentation objects.",
        ],
        "outside_scope": [
            "real histology segmentation",
            "true cell-level Tangram optimization",
            "atlas selection and biological interpretation on real tissue",
        ],
    }
    mapping.write_h5ad(outdir / "mapping_matrix.h5ad")

    prediction_counts = assignment_frame["dominant_label"].value_counts().to_dict()
    max_entropy_spot = spot_metric_frame.sort_values("assignment_entropy", ascending=False).iloc[0]

    report_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}.",
                f"Computed NNLS mapping weights for {len(spots)} spots against {len(reference_labels)} reference centroids.",
                "Segmentation-aware projection converted per-spot weights into integerized cell-slot assignments using provided spot cell counts.",
            ],
        },
        {
            "name": "Input QC",
            "bullets": [
                f"Shared genes: {len(shared_genes)} total, with {len(fit_genes)} fit genes and {len(holdout_genes)} held-out genes.",
                f"Reference labels: {', '.join(reference_labels)}.",
                f"Segmentation counts per spot: {', '.join(f'{spot_id}={count}' for spot_id, count in zip(spot_ids, segmentation_counts))}.",
            ],
        },
        {
            "name": "Assignments",
            "bullets": [
                f"Dominant label counts: {', '.join(f'{label}={count}' for label, count in sorted(prediction_counts.items()))}.",
                f"Most ambiguous spot: {max_entropy_spot['spot_id']} with entropy {max_entropy_spot['assignment_entropy']:.3f}.",
                f"{assignment_frame.iloc[0]['spot_id']}={assignment_frame.iloc[0]['dominant_label']} ({assignment_frame.iloc[0]['dominant_score']:.3f}); "
                f"{assignment_frame.iloc[1]['spot_id']}={assignment_frame.iloc[1]['dominant_label']} ({assignment_frame.iloc[1]['dominant_score']:.3f}); "
                f"{assignment_frame.iloc[2]['spot_id']} remains mixed with secondary support {assignment_frame.iloc[2]['secondary_score']:.3f}.",
            ],
        },
        {
            "name": "Validation",
            "bullets": [
                "Held-out genes were predicted from mapping weights and reference centroids rather than copied from the toy input.",
                "; ".join(
                    f"{gene}: r={metrics['pearson_r']:.3f}, rmse={metrics['rmse']:.3f}"
                    for gene, metrics in holdout_gene_metrics.items()
                ),
                f"Median fit-gene reconstruction MAE: {spot_metric_frame['fit_reconstruction_mae'].median():.3f}.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter approximates Tangram constrained mode with centroid-level NNLS and tiny synthetic data.",
                "Spot cell counts are supplied directly in the toy input instead of being extracted from histology segmentation.",
                "The outputs validate a portable shell workflow and QC contract, not a production biological analysis.",
            ],
        },
    ]
    write_markdown(outdir / "deconvolution_report.md", "Deconvolution Report", report_sections)

    assignment_frame.to_csv(outdir / "cell_assignment.tsv", sep="\t", index=False)
    spot_metric_frame.to_csv(outdir / "qc_spot_metrics.tsv", sep="\t", index=False)
    segment_frame.to_csv(outdir / "qc_segment_projection.tsv", sep="\t", index=False)
    holdout_frame.to_csv(outdir / "qc_holdout_predictions.tsv", sep="\t", index=False)

    gene_qc = {
        "run_label": payload["run_label"],
        "shared_gene_count": len(shared_genes),
        "fit_gene_count": len(fit_genes),
        "holdout_gene_count": len(holdout_genes),
        "shared_genes": shared_genes,
        "fit_genes": fit_genes,
        "holdout_genes": holdout_genes,
        "reference_gene_count": len(reference_genes),
        "spatial_gene_count": len(spatial_genes),
        "reference_labels": reference_labels,
        "spot_ids": spot_ids,
    }
    (outdir / "qc_gene_intersection.json").write_text(
        json.dumps(gene_qc, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path.resolve()),
        "outdir": str(outdir.resolve()),
        "method_steps": [
            "gene_intersection",
            "nnls_mapping",
            "segmentation_projection",
            "holdout_gene_validation",
        ],
        "spot_count": len(spots),
        "reference_profile_count": len(reference_labels),
        "shared_gene_count": len(shared_genes),
        "fit_gene_count": len(fit_genes),
        "holdout_gene_count": len(holdout_genes),
        "spatial_library_sizes": [float(value) for value in spatial_library_sizes.tolist()],
        "reference_library_sizes": [float(value) for value in reference_library_sizes.tolist()],
        "prediction_counts": {key: int(value) for key, value in prediction_counts.items()},
        "holdout_gene_metrics": holdout_gene_metrics,
        "written_files": sorted(
            {path.name for path in outdir.iterdir() if path.is_file()} | {"run_summary.json"}
        ),
    }
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    validate_outputs(SKILL_DIR, outdir)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the Tangram-style cell deconvolution starter on tiny local toy data."
    )
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    result = build_deconvolution(payload, input_path=args.input, outdir=args.outdir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
