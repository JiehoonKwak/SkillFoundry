#!/usr/bin/env python3
"""Materialize runnable toy exercise scaffolds for experiment-only sc_skills."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_ROOT = ROOT / "experiments" / "sc_skills"
MANIFEST_PATH = EXPERIMENT_ROOT / "batch_design_manifest.json"
SHARED_ROOT = EXPERIMENT_ROOT / "_shared"


PRIMARY_RUNNERS = {
    "cellphonedb-analysis": "run_exercise.py",
    "database-query": "run_database_query.py",
    "liana-analysis": "rank_and_report.py",
    "ligand-receptor-discovery": "run_exercise.py",
}


EXTRA_RUNNERS = {
    "cellphonedb-analysis": ["rank_and_report.py"],
    "ligand-receptor-discovery": ["rank_and_report.py"],
    "liana-analysis": ["run_exercise.py"],
}


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def runtime_module_text() -> str:
    return textwrap.dedent(
        """
        from __future__ import annotations

        import argparse
        import json
        from pathlib import Path

        import anndata as ad
        import h5py
        import numpy as np
        import pandas as pd


        def _write_table(path: Path, spec: dict) -> None:
            frame = pd.DataFrame(spec["rows"], columns=spec["columns"])
            frame.to_csv(path, sep="\\t", index=False)


        def _write_markdown(path: Path, spec: dict) -> None:
            lines: list[str] = []
            title = spec.get("title")
            if title:
                lines.append(title)
                lines.append("")
            for section in spec.get("sections", []):
                lines.append(f"## {section['heading']}")
                lines.append("")
                lines.append(section["body"])
                lines.append("")
            path.write_text("\\n".join(lines).rstrip() + "\\n", encoding="utf-8")


        def _write_h5ad(path: Path, spec: dict) -> None:
            obs = pd.DataFrame(index=spec["obs_names"])
            for column, values in spec.get("obs_columns", {}).items():
                obs[column] = values
            var = pd.DataFrame(index=spec["var_names"])
            for column, values in spec.get("var_columns", {}).items():
                var[column] = values
            adata = ad.AnnData(X=np.asarray(spec["X"], dtype=float), obs=obs, var=var)
            adata.write_h5ad(path)


        def _write_h5mu(path: Path, spec: dict) -> None:
            with h5py.File(path, "w") as handle:
                handle.attrs["encoding-type"] = "MuData"
                handle.attrs["encoding-version"] = "0.1.0"
                handle.create_dataset("latent", data=np.asarray(spec["latent"], dtype=float))
                handle.create_dataset("obs_names", data=np.asarray(spec["obs_names"], dtype="S"))
                handle.create_dataset("latent_names", data=np.asarray(spec["latent_names"], dtype="S"))
                mods = handle.create_group("mod")
                for modality_name, modality in spec.get("modalities", {}).items():
                    group = mods.create_group(modality_name)
                    group.create_dataset("shape", data=np.asarray(modality["shape"], dtype=int))
                    group.attrs["summary"] = modality["summary"]


        def _validate_table(path: Path, spec: dict) -> None:
            frame = pd.read_csv(path, sep="\\t")
            expected = spec["columns"]
            if list(frame.columns) != expected:
                raise AssertionError(f"{path} columns {list(frame.columns)} != {expected}")
            if len(frame) < spec.get("min_rows", 1):
                raise AssertionError(f"{path} has too few rows: {len(frame)}")


        def _validate_markdown(path: Path, spec: dict) -> None:
            text = path.read_text(encoding="utf-8")
            title = spec.get("title")
            if title and title not in text:
                raise AssertionError(f"{path} missing title {title!r}")
            for section in spec.get("sections", []):
                heading = f"## {section['heading']}"
                if heading not in text:
                    raise AssertionError(f"{path} missing heading {heading!r}")


        def _validate_h5ad(path: Path, spec: dict) -> None:
            adata = ad.read_h5ad(path)
            expected_shape = (len(spec["obs_names"]), len(spec["var_names"]))
            if adata.shape != expected_shape:
                raise AssertionError(f"{path} shape {adata.shape} != {expected_shape}")


        def _validate_h5mu(path: Path, spec: dict) -> None:
            with h5py.File(path, "r") as handle:
                if "latent" not in handle:
                    raise AssertionError(f"{path} missing latent dataset")
                if "mod" not in handle:
                    raise AssertionError(f"{path} missing modality group")
                latent_shape = tuple(handle["latent"].shape)
                expected = (len(spec["obs_names"]), len(spec["latent_names"]))
                if latent_shape != expected:
                    raise AssertionError(f"{path} latent shape {latent_shape} != {expected}")


        def load_example(task_root: Path) -> dict:
            return json.loads((task_root / "examples" / "toy_input.json").read_text(encoding="utf-8"))


        def run_task(task_root: Path, output_dir: Path) -> None:
            spec = load_example(task_root)
            output_dir.mkdir(parents=True, exist_ok=True)
            for name, table_spec in spec.get("tables", {}).items():
                _write_table(output_dir / name, table_spec)
            for name, markdown_spec in spec.get("markdown", {}).items():
                _write_markdown(output_dir / name, markdown_spec)
            for name, h5ad_spec in spec.get("h5ad", {}).items():
                _write_h5ad(output_dir / name, h5ad_spec)
            for name, h5mu_spec in spec.get("h5mu", {}).items():
                _write_h5mu(output_dir / name, h5mu_spec)


        def validate_task(task_root: Path, output_dir: Path) -> None:
            spec = load_example(task_root)
            for name, table_spec in spec.get("tables", {}).items():
                _validate_table(output_dir / name, table_spec)
            for name, markdown_spec in spec.get("markdown", {}).items():
                _validate_markdown(output_dir / name, markdown_spec)
            for name, h5ad_spec in spec.get("h5ad", {}).items():
                _validate_h5ad(output_dir / name, h5ad_spec)
            for name, h5mu_spec in spec.get("h5mu", {}).items():
                _validate_h5mu(output_dir / name, h5mu_spec)


        def cli(task_root: Path, mode: str) -> int:
            parser = argparse.ArgumentParser()
            parser.add_argument("--output-dir", type=Path, required=True)
            args = parser.parse_args()
            if mode == "run":
                run_task(task_root, args.output_dir)
            else:
                validate_task(task_root, args.output_dir)
            return 0
        """
    ).strip() + "\n"


def runner_wrapper_text(task_slug: str, mode: str) -> str:
    return textwrap.dedent(
        f"""
        #!/usr/bin/env python3
        from __future__ import annotations

        import importlib.util
        import sys
        from pathlib import Path

        TASK_ROOT = Path(__file__).resolve().parents[1]
        RUNTIME_PATH = TASK_ROOT.parent / "_shared" / "runtime.py"

        spec = importlib.util.spec_from_file_location("sc_skill_runtime", RUNTIME_PATH)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)

        if __name__ == "__main__":
            raise SystemExit(module.cli(TASK_ROOT, "{mode}"))
        """
    ).strip() + "\n"


def local_test_text(primary_runner: str) -> str:
    return textwrap.dedent(
        f"""
        from __future__ import annotations

        import subprocess
        import tempfile
        import unittest
        from pathlib import Path


        SKILL_ROOT = Path(__file__).resolve().parents[1]


        class ExperimentContractTest(unittest.TestCase):
            def test_toy_contract(self) -> None:
                with tempfile.TemporaryDirectory() as tmpdir:
                    output_dir = Path(tmpdir)
                    run_completed = subprocess.run(
                        ["python3", str(SKILL_ROOT / "scripts" / "{primary_runner}"), "--output-dir", str(output_dir)],
                        cwd=SKILL_ROOT.parents[4],
                        check=False,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(run_completed.returncode, 0, run_completed.stderr)
                    validate_completed = subprocess.run(
                        ["python3", str(SKILL_ROOT / "scripts" / "validate_outputs.py"), "--output-dir", str(output_dir)],
                        cwd=SKILL_ROOT.parents[4],
                        check=False,
                        capture_output=True,
                        text=True,
                    )
                    self.assertEqual(validate_completed.returncode, 0, validate_completed.stderr)


        if __name__ == "__main__":
            unittest.main()
        """
    ).strip() + "\n"


def sections_to_markdown(title: str, headings: list[str]) -> dict:
    return {
        "title": title,
        "sections": [{"heading": heading, "body": f"Toy summary for {heading.lower()}."} for heading in headings],
    }


def table(columns: list[str], rows: list[list[object]]) -> dict:
    return {"columns": columns, "rows": rows, "min_rows": len(rows)}


def h5ad(X: list[list[float]], obs_names: list[str], var_names: list[str], obs_columns: dict | None = None) -> dict:
    return {
        "X": X,
        "obs_names": obs_names,
        "var_names": var_names,
        "obs_columns": obs_columns or {},
    }


def h5mu(latent: list[list[float]], obs_names: list[str], latent_names: list[str], modalities: dict[str, dict]) -> dict:
    return {
        "latent": latent,
        "obs_names": obs_names,
        "latent_names": latent_names,
        "modalities": modalities,
    }


TASK_EXAMPLES = {
    "annotation": {
        "tables": {
            "annotation_table.tsv": table(
                [
                    "cell_id",
                    "sample_id",
                    "x_coord",
                    "y_coord",
                    "predicted_label",
                    "label_source",
                    "label_confidence",
                    "broad_label",
                    "marker_status",
                    "niche_label",
                    "niche_evidence",
                    "notes",
                ],
                [
                    ["cell_1", "slide_a", 10.0, 5.0, "T cell", "CellTypist", 0.92, "immune", "supported", "immune_border", "neighbor_enrichment", "toy"],
                    ["cell_2", "slide_a", 12.5, 6.0, "Endothelial", "marker_review", 0.88, "stromal", "supported", "vascular_stroma", "marker_plus_graph", "toy"],
                ],
            ),
        },
        "markdown": {
            "marker_evidence.md": sections_to_markdown("# Marker Evidence", ["Run context", "Reference source", "Marker review", "Conflicts and resolutions", "Remaining uncertainty"]),
            "niche_summary.md": sections_to_markdown("# Niche Summary", ["Run context", "Spatial graph", "Neighborhood validation", "Niche labels", "Caveats"]),
        },
    },
    "cell-cell-communication": {
        "tables": {
            "communication_results.tsv": table(
                ["source_cell", "target_cell", "interacting_pair", "consensus_score", "spatial_support", "notes"],
                [["T_cell", "Fibroblast", "CXCL12-CXCR4", 0.81, "adjacent", "toy"], ["Myeloid", "Endothelial", "TNF-TNFRSF1A", 0.74, "same_niche", "toy"]],
            ),
            "priority_pairs.tsv": table(
                ["source_cell", "target_cell", "interacting_pair", "priority_tier", "evidence_summary"],
                [["T_cell", "Fibroblast", "CXCL12-CXCR4", "high", "consensus+spatial"], ["Myeloid", "Endothelial", "TNF-TNFRSF1A", "medium", "ranked"]],
            ),
        },
        "markdown": {
            "interpretation_report.md": sections_to_markdown("# Communication Interpretation", ["Run context", "Primary method", "Filtering", "Prioritized biology", "Caveats"]),
        },
    },
    "cell-deconvolution": {
        "h5ad": {
            "mapping_matrix.h5ad": h5ad([[0.8, 0.2], [0.1, 0.9]], ["spot_1", "spot_2"], ["T_cell", "Endothelial"], {"sample_id": ["slide_a", "slide_a"]}),
        },
        "tables": {
            "cell_assignment.tsv": table(
                ["spot_id", "assigned_label", "assignment_score", "qc_status"],
                [["spot_1", "T_cell", 0.8, "pass"], ["spot_2", "Endothelial", 0.9, "pass"]],
            ),
        },
        "markdown": {
            "deconvolution_report.md": sections_to_markdown("# Deconvolution Report", ["Run context", "Inputs", "Mapping setup", "Assignment QC", "Caveats"]),
        },
    },
    "cellphonedb-analysis": {
        "tables": {
            "cellphonedb_significant_means.tsv": table(
                ["source_cell", "target_cell", "interacting_pair", "significant_mean", "pvalue"],
                [["T_cell", "Fibroblast", "CXCL12-CXCR4", 0.42, 0.004], ["Myeloid", "Endothelial", "TNF-TNFRSF1A", 0.38, 0.01]],
            ),
            "interaction_ranked.tsv": table(
                ["interacting_pair", "priority_score", "priority_tier", "notes"],
                [["CXCL12-CXCR4", 0.91, "high", "toy"], ["TNF-TNFRSF1A", 0.72, "medium", "toy"]],
            ),
        },
        "markdown": {
            "cellphonedb_report.md": sections_to_markdown("# CellPhoneDB Report", ["Run mode", "Input summary", "Evidence summary", "Top interaction programs", "Caveats and blockers"]),
        },
    },
    "database-query": {
        "tables": {
            "query_results.tsv": table(
                ["query_id", "resolved_id", "source", "entity_type", "label"],
                [["gene_1", "ENSG000001", "Ensembl", "gene", "BRAF"], ["protein_1", "P15056", "UniProt", "protein", "BRAF_HUMAN"]],
            ),
        },
        "markdown": {
            "source_provenance.md": sections_to_markdown("# Source Provenance", ["Run context", "Sources queried", "Resolution choices"]),
            "resolution_notes.md": sections_to_markdown("# Resolution Notes", ["Run context", "Alias normalization", "Ambiguities", "Caveats"]),
        },
    },
    "gene-imputation": {
        "h5ad": {
            "imputed_expression.h5ad": h5ad([[1.2, 0.5, 0.8], [0.7, 1.5, 0.4]], ["spot_1", "spot_2"], ["GENE1", "GENE2", "GENE3"], {"sample_id": ["slide_a", "slide_a"]}),
        },
        "tables": {
            "heldout_metrics.tsv": table(
                ["metric_name", "value", "higher_is_better", "qc_status"],
                [["pearson_r", 0.82, True, "pass"], ["rmse", 0.18, False, "pass"]],
            ),
        },
        "markdown": {
            "imputation_report.md": sections_to_markdown("# Imputation Report", ["Run context", "Inputs", "Held-out evaluation", "Caveats"]),
        },
    },
    "liana-analysis": {
        "tables": {
            "liana_rankings.tsv": table(
                ["interacting_pair", "method", "rank", "consensus_score"],
                [["CXCL12-CXCR4", "liana", 1, 0.93], ["TNF-TNFRSF1A", "liana", 2, 0.78]],
            ),
            "resource_overlap.tsv": table(
                ["interacting_pair", "resource_count", "resources"],
                [["CXCL12-CXCR4", 3, "CellPhoneDB;OmniPath;LIANA"], ["TNF-TNFRSF1A", 2, "CellPhoneDB;OmniPath"]],
            ),
        },
        "markdown": {
            "liana_report.md": sections_to_markdown("# LIANA Report", ["Run context", "Input summary", "Ranking logic", "Top interactions", "Caveats"]),
        },
    },
    "ligand-receptor-discovery": {
        "tables": {
            "candidate_pairs.tsv": table(
                ["source_cell", "target_cell", "ligand", "receptor", "priority_score"],
                [["T_cell", "Fibroblast", "CXCL12", "CXCR4", 0.94], ["Myeloid", "Endothelial", "TNF", "TNFRSF1A", 0.71]],
            ),
            "evidence_table.tsv": table(
                ["pair_id", "resource_support", "expression_support", "spatial_support", "notes"],
                [["CXCL12-CXCR4", "high", "high", "adjacent", "toy"], ["TNF-TNFRSF1A", "medium", "medium", "same_region", "toy"]],
            ),
        },
        "markdown": {
            "discovery_summary.md": sections_to_markdown("# Discovery Summary", ["Run context", "Search space", "Prioritization", "Caveats"]),
        },
    },
    "mapping-validation": {
        "tables": {
            "validation_metrics.tsv": table(
                ["metric_name", "value", "qc_status", "notes"],
                [["shared_gene_fraction", 0.67, "pass", "toy"], ["neighbor_label_agreement", 0.81, "pass", "toy"]],
            ),
        },
        "markdown": {
            "holdout_summary.md": sections_to_markdown("# Holdout Summary", ["Run context", "Held-out genes", "Performance summary"]),
            "mapping_validation_report.md": sections_to_markdown("# Mapping Validation Report", ["Run context", "Input QC", "Validation metrics", "Caveats"]),
        },
    },
    "multimodal-integration": {
        "h5mu": {
            "integrated_latent.h5mu": h5mu(
                [[0.1, 0.8], [0.4, 0.6], [0.9, 0.2]],
                ["cell_1", "cell_2", "cell_3"],
                ["LV1", "LV2"],
                {"rna": {"shape": [3, 4], "summary": "toy RNA counts"}, "atac": {"shape": [3, 5], "summary": "toy chromatin accessibility"}},
            ),
        },
        "tables": {
            "modality_qc.tsv": table(
                ["metric_name", "value", "qc_status"],
                [["batch_mixing", 0.79, "pass"], ["modality_balance", 0.83, "pass"]],
            ),
        },
        "markdown": {
            "integration_report.md": sections_to_markdown("# Integration Report", ["Run context", "Modalities", "Latent QC", "Interpretation", "Caveats"]),
        },
    },
    "panel-design": {
        "tables": {
            "panel_candidates.tsv": table(
                ["gene", "target_label", "specificity_score", "platform_fit"],
                [["CXCL13", "B_cell", 0.91, "Xenium"], ["PECAM1", "Endothelial", 0.88, "MERFISH"]],
            ),
        },
        "markdown": {
            "panel_rationale.md": sections_to_markdown("# Panel Rationale", ["Run context", "Selection rules", "Chosen genes", "Caveats"]),
            "platform_notes.md": sections_to_markdown("# Platform Notes", ["Run context", "Platform constraints", "Design follow-up"]),
        },
    },
    "sequence-analysis": {
        "tables": {
            "sequence_summary.tsv": table(
                ["query_id", "molecule_type", "matched_symbol", "matched_accession", "best_hit_label"],
                [["query_1", "DNA", "BRAF", "NM_004333.6", "BRAF transcript"], ["query_1_protein", "protein", "BRAF", "P15056", "BRAF_HUMAN"]],
            ),
            "feature_annotations.tsv": table(
                ["query_id", "feature_source", "feature_type", "feature_id", "start", "end", "label"],
                [["query_1", "Ensembl", "exon", "ENSE0001", 100, 210, "toy exon"], ["query_1_protein", "UniProt", "domain", "DOMAIN1", 50, 120, "kinase"]],
            ),
            "primer_candidates.tsv": table(
                ["query_id", "pair_rank", "left_primer", "right_primer", "product_size", "left_tm_c", "right_tm_c"],
                [["query_1", 1, "ACGTACGTACGT", "TGCATGCATGCA", 180, 60.1, 60.4], ["query_1", 2, "CGTACGTACGTA", "CATGCATGCATG", 210, 59.8, 60.0]],
            ),
        },
    },
    "spatial-deconvolution": {
        "h5ad": {
            "cell_type_abundance.h5ad": h5ad([[0.7, 0.3], [0.2, 0.8]], ["spot_1", "spot_2"], ["T_cell", "Fibroblast"], {"sample_id": ["slide_a", "slide_a"]}),
        },
        "tables": {
            "model_comparison.tsv": table(
                ["model_name", "metric_name", "value", "qc_status"],
                [["cell2location", "marker_agreement", 0.84, "pass"], ["DestVI", "marker_agreement", 0.79, "pass"]],
            ),
        },
        "markdown": {
            "deconvolution_summary.md": sections_to_markdown("# Spatial Deconvolution Summary", ["Run context", "Models compared", "Abundance summary", "Caveats"]),
        },
    },
    "spatial-domain-detection": {
        "tables": {
            "domain_labels.tsv": table(
                ["cell_id", "spagcn_domain", "graphst_domain", "consensus_domain"],
                [["cell_1", "tumor_core", "tumor_core", "tumor_core"], ["cell_2", "immune_border", "immune_border", "immune_border"]],
            ),
            "domain_markers.tsv": table(
                ["domain_label", "marker_gene", "effect_size", "notes"],
                [["tumor_core", "EPCAM", 1.8, "toy"], ["immune_border", "CXCL13", 1.5, "toy"]],
            ),
        },
        "markdown": {
            "domain_detection_report.md": sections_to_markdown("# Domain Detection Report", ["Run context", "Spatial graph", "Domain agreement", "Marker evidence", "Caveats"]),
        },
    },
    "spatial-mapping": {
        "tables": {
            "mapped_labels.tsv": table(
                ["observation_id", "x_coord", "y_coord", "primary_label", "primary_score", "score_margin", "secondary_label", "mapping_method", "gene_overlap_count", "gene_overlap_fraction", "marker_review_status", "review_notes"],
                [["spot_1", 10.0, 5.0, "T_cell", 0.91, 0.25, "Myeloid", "Tangram", 4, 0.8, "accepted", "toy"], ["spot_2", 13.0, 6.5, "Endothelial", 0.88, 0.2, "Fibroblast", "Tangram", 4, 0.8, "accepted", "toy"]],
            ),
            "mapping_scores.tsv": table(
                ["metric_category", "metric_name", "scope", "target_id", "method", "value", "higher_is_better", "suggested_pass_floor", "qc_status", "notes"],
                [["gene_intersection", "shared_gene_fraction", "run", "all", "Tangram", 0.8, True, 0.5, "pass", "toy"], ["mapping", "mean_primary_score", "run", "all", "Tangram", 0.895, True, 0.6, "pass", "toy"]],
            ),
        },
        "markdown": {
            "spatial_mapping_report.md": sections_to_markdown("# Spatial Mapping Report", ["Run context", "Inputs and preprocessing", "Gene intersection QC", "Tangram mapping", "Label transfer review", "Caveats and next actions"]),
        },
    },
    "squidpy-analysis": {
        "tables": {
            "spatial_stats.tsv": table(
                ["metric_name", "value", "qc_status"],
                [["moran_i", 0.41, "pass"], ["co_occurrence_score", 0.66, "pass"]],
            ),
            "neighbor_enrichment.tsv": table(
                ["label_a", "label_b", "enrichment_score", "notes"],
                [["T_cell", "Fibroblast", 1.9, "toy"], ["Myeloid", "Endothelial", 1.4, "toy"]],
            ),
        },
        "markdown": {
            "squidpy_report.md": sections_to_markdown("# Squidpy Report", ["Run context", "Spatial graph", "Statistics", "Caveats"]),
        },
    },
    "trajectory-inference": {
        "tables": {
            "trajectory_coordinates.tsv": table(
                ["cell_id", "pseudotime", "lineage_label", "embedding_x", "embedding_y"],
                [["cell_1", 0.1, "lineage_a", 0.2, 1.1], ["cell_2", 0.8, "lineage_a", 1.0, 0.3]],
            ),
            "fate_probabilities.tsv": table(
                ["cell_id", "terminal_state", "probability"],
                [["cell_1", "state_1", 0.72], ["cell_2", "state_2", 0.81]],
            ),
        },
        "markdown": {
            "trajectory_report.md": sections_to_markdown("# Trajectory Report", ["Run context", "Preprocessing", "Trajectory outputs", "Caveats"]),
        },
    },
}


def build_example_spec(task_slug: str) -> dict:
    return TASK_EXAMPLES[task_slug]


def update_skill_text(task_dir: Path, task_slug: str, primary_runner: str) -> None:
    skill_path = task_dir / "SKILL.md"
    text = skill_path.read_text(encoding="utf-8")
    text = text.replace("ref-skill", "sc-skill")
    if "## Exercise path" not in text:
        exercise_block = textwrap.dedent(
            f"""

            ## Exercise path

            Use the bundled deterministic exercise to verify the contract without requiring the full external runtime:

            ```bash
            python3 experiments/sc_skills/{task_slug}/scripts/{primary_runner} --output-dir /tmp/{task_slug}-exercise
            python3 experiments/sc_skills/{task_slug}/scripts/validate_outputs.py --output-dir /tmp/{task_slug}-exercise
            python3 -m unittest discover -s experiments/sc_skills/{task_slug}/tests -p 'test_*.py'
            ```
            """
        )
        text = text.rstrip() + "\n" + exercise_block
    skill_path.write_text(text.rstrip() + "\n", encoding="utf-8")


def materialize_task(task: dict) -> None:
    task_slug = task["task_slug"]
    task_dir = ROOT / task["target_dir"]
    primary_runner = PRIMARY_RUNNERS.get(task_slug, "run_exercise.py")
    example_spec = build_example_spec(task_slug)

    write_text(task_dir / "examples" / "toy_input.json", json.dumps(example_spec, indent=2) + "\n")
    write_text(task_dir / "scripts" / primary_runner, runner_wrapper_text(task_slug, "run"))
    for extra_runner in EXTRA_RUNNERS.get(task_slug, []):
        write_text(task_dir / "scripts" / extra_runner, runner_wrapper_text(task_slug, "run"))
    write_text(task_dir / "scripts" / "validate_outputs.py", runner_wrapper_text(task_slug, "validate"))
    write_text(task_dir / "tests" / "test_contract.py", local_test_text(primary_runner))
    update_skill_text(task_dir, task_slug, primary_runner)


def main() -> int:
    manifest = load_manifest()
    write_text(SHARED_ROOT / "runtime.py", runtime_module_text())
    for task in manifest["tasks"]:
        materialize_task(task)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
