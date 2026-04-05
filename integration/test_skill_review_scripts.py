from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRATCH = ROOT / "scratch" / "tests" / "skill-review"
REGISTRY = ROOT / "registry" / "skills.jsonl"


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=1800,
    )


class SkillReviewScriptsIntegrationTest(unittest.TestCase):
    def setUp(self) -> None:
        SCRATCH.mkdir(parents=True, exist_ok=True)

    def skill_count(self) -> int:
        return sum(1 for line in REGISTRY.read_text(encoding="utf-8").splitlines() if line.strip())

    def test_audit_skill_suite_reports_full_coverage(self) -> None:
        audit_json = SCRATCH / "audit.json"
        run_command(
            [
                "python3",
                "scripts/audit_skill_suite.py",
                "--json-out",
                str(audit_json),
            ]
        )
        payload = json.loads(audit_json.read_text(encoding="utf-8"))
        summary = payload["summary"]
        self.assertEqual(summary["skill_count"], self.skill_count())
        self.assertEqual(summary["hard_failure_count"], 0)
        self.assertEqual(summary["skills_with_smoke_target"], summary["skill_count"])
        self.assertEqual(summary["skills_with_repo_test_refs"], summary["skill_count"])

    def test_skill_smoke_matrix_dry_run_resolves_all_targets(self) -> None:
        matrix_json = SCRATCH / "smoke_matrix_dry_run.json"
        run_command(
            [
                "python3",
                "scripts/run_skill_smoke_matrix.py",
                "--dry-run",
                "--json-out",
                str(matrix_json),
            ]
        )
        payload = json.loads(matrix_json.read_text(encoding="utf-8"))
        summary = payload["summary"]
        self.assertEqual(summary["skill_count"], self.skill_count())
        self.assertEqual(summary["targets_resolved"], summary["skill_count"])
        self.assertEqual(summary["missing_targets"], [])

    def test_skill_advantage_subset_shows_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "skill_advantage_subset.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "geopandas-spatial-join",
                "--case",
                "papermill-parameterized-notebook",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_workflow_generation_agents_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "workflow_generation_agents_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "workflow-generation-agents-starter-summary",
                "--case",
                "workflow-generation-agents-starter-checklist",
                "--case",
                "workflow-generation-agents-starter-nested-output",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertEqual(payload["aggregate"]["skill"]["perfect_case_count"], 3)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], 3)

    def test_uniprot_sequence_feature_annotation_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "uniprot_sequence_feature_annotation_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "uniprot-sequence-feature-annotation-starter-p04637",
                "--case",
                "uniprot-sequence-feature-annotation-starter-p38398",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertEqual(payload["aggregate"]["skill"]["perfect_case_count"], 2)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], 2)

    def test_somatic_pipelines_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "somatic_pipelines_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "somatic-pipelines-starter-summary",
                "--case",
                "somatic-pipelines-starter-augmented",
                "--case",
                "somatic-pipelines-starter-nested-output",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertEqual(payload["aggregate"]["skill"]["perfect_case_count"], 3)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], 3)

    def test_skill_browser_mindmap_generation_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "skill_browser_mindmap_generation_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "skill-browser-mindmap-generation-starter-summary",
                "--case",
                "skill-browser-mindmap-generation-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertEqual(payload["aggregate"]["skill"]["perfect_case_count"], 2)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], 2)

    def test_reactome_event_summary_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "reactome_event_summary_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "reactome-event-summary-canonical",
                "--case",
                "reactome-event-summary-noisy-input",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertEqual(payload["aggregate"]["skill"]["perfect_case_count"], 2)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], 2)

    def test_rapidfuzz_skill_deduplication_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "rapidfuzz_skill_deduplication_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "rapidfuzz-skill-deduplication-registry-slice",
                "--case",
                "rapidfuzz-skill-deduplication-mixed-registry",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertGreater(payload["aggregate"]["skill"]["perfect_case_count"], 0)
        self.assertLess(payload["aggregate"]["baseline"]["perfect_case_count"], payload["aggregate"]["skill"]["perfect_case_count"])

    def test_nf_core_pipeline_list_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "nf_core_pipeline_list_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "nf-core-pipeline-list-pulled-three",
                "--case",
                "nf-core-pipeline-list-release-five",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_polygenic_risk_scoring_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "polygenic_risk_scoring_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "polygenic-risk-scoring-starter-summary",
                "--case",
                "polygenic-risk-scoring-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_rna_velocity_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "rna_velocity_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "rna-velocity-starter-canonical-summary",
                "--case",
                "rna-velocity-starter-nested-output",
                "--case",
                "rna-velocity-starter-promotion-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )
        self.assertGreater(
            payload["aggregate"]["skill"]["success_rate"],
            payload["aggregate"]["baseline"]["success_rate"],
        )

    def test_pseudobulk_analysis_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "pseudobulk_analysis_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "pseudobulk-analysis-starter-summary",
                "--case",
                "pseudobulk-analysis-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_qsar_property_prediction_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "qsar_property_prediction_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "qsar-property-prediction-starter-summary",
                "--case",
                "qsar-property-prediction-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_pydoe3_experimental_design_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "pydoe3_experimental_design_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "pydoe3-experimental-design-starter-summary",
                "--case",
                "pydoe3-experimental-design-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_numcodecs_compression_decompression_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "numcodecs_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "numcodecs-compression-decompression-canonical",
                "--case",
                "numcodecs-compression-decompression-nested-output",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_microscopy_pipelines_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "microscopy_pipelines_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "microscopy-pipelines-starter-summary",
                "--case",
                "microscopy-pipelines-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_pde_cfd_simulation_workflows_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "pde_cfd_simulation_workflows_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "pde-cfd-simulation-workflows-starter-summary",
                "--case",
                "pde-cfd-simulation-workflows-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_multi_modal_image_omics_integration_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "multi_modal_image_omics_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "multi-modal-image-omics-integration-starter-summary",
                "--case",
                "multi-modal-image-omics-integration-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_pathology_histology_workflows_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "pathology_histology_workflows_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "pathology-histology-workflows-starter-summary",
                "--case",
                "pathology-histology-workflows-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_langgraph_planning_execution_agent_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "langgraph_planning_execution_agent_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "langgraph-planning-execution-agent-single-cell-report",
                "--case",
                "langgraph-planning-execution-agent-literature-single-cell-report",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_openalex_literature_search_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "openalex_literature_search_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "openalex-literature-search-single-cell",
                "--case",
                "openalex-literature-search-spatial-transcriptomics",
                "--case",
                "openalex-literature-search-gwas-methods",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_scientific_summarization_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "scientific_summarization_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "scientific-summarization-starter-summary",
                "--case",
                "scientific-summarization-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_dash_scientific_dashboard_benchmark_case_shows_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "dash_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "dash-scientific-dashboard",
                "--case",
                "dash-scientific-dashboard-extended",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_mkdocs_summary_catalog_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "mkdocs_summary_catalog_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "mkdocs-summary-catalog-canonical",
                "--case",
                "mkdocs-summary-catalog-stale-rebuild",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_ebi_proteins_entry_summary_benchmark_case_shows_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "ebi_proteins_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "ebi-proteins-entry-summary-canonical",
                "--case",
                "ebi-proteins-entry-summary-normalized-input",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_ms_proteomics_preprocessing_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "ms_proteomics_preprocessing_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "ms-proteomics-preprocessing-starter-summary",
                "--case",
                "ms-proteomics-preprocessing-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_multimodal_fusion_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "multimodal_fusion_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "multimodal-fusion-starter-summary",
                "--case",
                "multimodal-fusion-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_metabolights_study_search_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "metabolights_study_search_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "metabolights-study-search-canonical",
                "--case",
                "metabolights-study-search-normalized-multi",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_ensembl_gene_lookup_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "ensembl_gene_lookup_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "ensembl-gene-lookup-canonical",
                "--case",
                "ensembl-gene-lookup-fallback",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_germline_pipelines_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "germline_pipelines_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "germline-pipelines-starter-canonical",
                "--case",
                "germline-pipelines-starter-mutated",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_gpu_jobs_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "gpu_jobs_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "gpu-jobs-starter-summary",
                "--case",
                "gpu-jobs-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_gwas_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "gwas_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "gwas-starter-summary",
                "--case",
                "gwas-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_protocol_and_workflow_extraction_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "protocol_and_workflow_extraction_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "protocol-and-workflow-extraction-starter-summary",
                "--case",
                "protocol-and-workflow-extraction-starter-nested-output",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_missing_data_handling_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "missing_data_handling_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "missing-data-handling-starter-summary",
                "--case",
                "missing-data-handling-starter-checklist",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_long_read_genomics_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "long_read_genomics_starter_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "long-read-genomics-starter-summary",
                "--case",
                "long-read-genomics-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_materials_benchmark_datasets_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "materials_benchmark_datasets_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "materials-benchmark-datasets-starter-summary",
                "--case",
                "materials-benchmark-datasets-starter-augmented",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_isoform_transcript_level_analysis_starter_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "isoform_transcript_level_analysis_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "isoform-transcript-level-analysis-starter-summary",
                "--case",
                "isoform-transcript-level-analysis-starter-nested-output",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_ncbi_pubmed_search_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "ncbi_pubmed_search_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "ncbi-pubmed-search-single-cell-top-hit",
                "--case",
                "ncbi-pubmed-search-single-cell-top-three",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )

    def test_quickgo_term_search_benchmark_cases_show_better_deliverable_rate(self) -> None:
        benchmark_json = SCRATCH / "quickgo_term_search_skill_advantage.json"
        run_command(
            [
                "python3",
                "scripts/benchmark_skill_advantage.py",
                "--case",
                "quickgo-term-search-apoptosis-local-fixture",
                "--case",
                "quickgo-term-search-cell-cycle-local-fixture",
                "--json-out",
                str(benchmark_json),
            ]
        )
        payload = json.loads(benchmark_json.read_text(encoding="utf-8"))
        self.assertTrue(payload["summary"]["skill_better_on_deliverable_rate"])
        self.assertGreater(
            payload["aggregate"]["skill"]["average_deliverable_rate"],
            payload["aggregate"]["baseline"]["average_deliverable_rate"],
        )


if __name__ == "__main__":
    unittest.main()
