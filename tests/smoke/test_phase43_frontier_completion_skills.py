from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCIENTIFIC_PYTHON = ROOT / "slurm" / "envs" / "scientific-python" / "bin" / "python"
MATERIALS_PYTHON = ROOT / "slurm" / "envs" / "materials" / "bin" / "python"
CHEMTOOLS_PYTHON = ROOT / "slurm" / "envs" / "chemtools" / "bin" / "python"


class Phase43FrontierCompletionSkillTests(unittest.TestCase):
    def test_scientific_knowledge_frontier_skills(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            triage_out = tmp / "triage.json"
            review_out = tmp / "review.json"
            captions_out = tmp / "captions.json"
            links_out = tmp / "links.json"
            benchmark_out = tmp / "benchmark.json"
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/semantic-scholar-paper-triage-starter/scripts/run_semantic_scholar_paper_triage.py",
                    "--input",
                    "skills/scientific-knowledge/semantic-scholar-paper-triage-starter/examples/candidate_papers.json",
                    "--query",
                    "single-cell RNA-seq atlas integration",
                    "--out",
                    str(triage_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/semantic-scholar-review-paper-mining-starter/scripts/run_semantic_scholar_review_mining.py",
                    "--input",
                    "skills/scientific-knowledge/semantic-scholar-review-paper-mining-starter/examples/paper_metadata.json",
                    "--out",
                    str(review_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/figure-table-caption-extraction-starter/scripts/run_figure_table_caption_extraction.py",
                    "--input",
                    "skills/scientific-knowledge/figure-table-caption-extraction-starter/examples/paper_excerpt.txt",
                    "--out",
                    str(captions_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/dataset-code-link-extraction-starter/scripts/run_dataset_code_link_extraction.py",
                    "--input",
                    "skills/scientific-knowledge/dataset-code-link-extraction-starter/examples/paper_text.md",
                    "--out",
                    str(links_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/benchmark-table-mining-starter/scripts/run_benchmark_table_mining.py",
                    "--input",
                    "skills/scientific-knowledge/benchmark-table-mining-starter/examples/benchmark_note.md",
                    "--out",
                    str(benchmark_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            triage = json.loads(triage_out.read_text(encoding="utf-8"))
            review = json.loads(review_out.read_text(encoding="utf-8"))
            captions = json.loads(captions_out.read_text(encoding="utf-8"))
            links = json.loads(links_out.read_text(encoding="utf-8"))
            benchmark = json.loads(benchmark_out.read_text(encoding="utf-8"))
            self.assertEqual(triage["top_candidates"][0]["paper_id"], "P1")
            self.assertEqual(review["review_paper_count"], 3)
            self.assertEqual(captions["table_caption_count"], 2)
            self.assertEqual(links["url_count"], 4)
            self.assertEqual(benchmark["best_method"]["method"], "scVI")

    def test_data_and_genomics_frontier_skills(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            harmonized_tsv = tmp / "harmonized.tsv"
            harmonized_out = tmp / "harmonized.json"
            toy_bundle = tmp / "toy_bundle"
            toy_bundle_out = tmp / "toy_bundle.json"
            cooler_out = tmp / "cooler.json"
            cooler_store = tmp / "toy.cool"
            bcftools_out = tmp / "bcftools.json"
            bcftools_vcf = tmp / "toy_variants.filtered.vcf.gz"
            subprocess.run(
                [
                    "python3",
                    "skills/data-acquisition-and-dataset-handling/metadata-harmonization-starter/scripts/run_metadata_harmonization.py",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/metadata-harmonization-starter/examples/cohort_a.tsv",
                    "--input",
                    "skills/data-acquisition-and-dataset-handling/metadata-harmonization-starter/examples/cohort_b.tsv",
                    "--mapping",
                    "skills/data-acquisition-and-dataset-handling/metadata-harmonization-starter/examples/column_mapping.json",
                    "--out-tsv",
                    str(harmonized_tsv),
                    "--summary-out",
                    str(harmonized_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/data-acquisition-and-dataset-handling/synthetic-toy-dataset-generator-starter/scripts/generate_synthetic_toy_dataset.py",
                    "--sample-count",
                    "6",
                    "--feature-count",
                    "4",
                    "--seed",
                    "17",
                    "--out-dir",
                    str(toy_bundle),
                    "--summary-out",
                    str(toy_bundle_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    str(SCIENTIFIC_PYTHON),
                    "skills/epigenomics-and-chromatin/cooler-hic-matrix-summary-starter/scripts/run_cooler_hic_matrix_summary.py",
                    "--out",
                    str(cooler_out),
                    "--cooler-out",
                    str(cooler_store),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            subprocess.run(
                [
                    "python3",
                    "skills/genomics/bcftools-variant-filtering-starter/scripts/run_bcftools_variant_filtering.py",
                    "--input",
                    "skills/genomics/bcftools-variant-filtering-starter/examples/toy_variants.vcf",
                    "--filtered-vcf-out",
                    str(bcftools_vcf),
                    "--out",
                    str(bcftools_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            harmonized = json.loads(harmonized_out.read_text(encoding="utf-8"))
            toy_bundle_summary = json.loads(toy_bundle_out.read_text(encoding="utf-8"))
            cooler = json.loads(cooler_out.read_text(encoding="utf-8"))
            bcftools = json.loads(bcftools_out.read_text(encoding="utf-8"))
            self.assertEqual(harmonized["row_count"], 4)
            self.assertEqual(toy_bundle_summary["matrix_shape"], [4, 6])
            self.assertEqual(cooler["total_contact_count"], 38)
            self.assertEqual(bcftools["kept_ids"], ["varA", "varC"])

    def test_matminer_toy_property_prediction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "matminer.json"
            subprocess.run(
                [
                    str(MATERIALS_PYTHON),
                    "skills/materials-science-and-engineering/matminer-toy-property-prediction-starter/scripts/run_matminer_toy_property_prediction.py",
                    "--input",
                    "skills/materials-science-and-engineering/matminer-toy-property-prediction-starter/examples/toy_materials.tsv",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_count"], 6)
            self.assertEqual(payload["training_mae"], 0.0)

    def test_deepchem_molgraph_featurization(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "deepchem_molgraph.json"
            subprocess.run(
                [
                    str(CHEMTOOLS_PYTHON),
                    "skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/scripts/featurize_molecules.py",
                    "--input",
                    "skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/examples/molecules.tsv",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["molecule_count"], 2)
            self.assertGreater(payload["graphs"][0]["node_feature_count"], 0)

    def test_generated_frontier_starters_emit_summaries(self) -> None:
        skills = [
            json.loads(line)
            for line in (ROOT / "registry" / "skills.jsonl").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        generated = [skill for skill in skills if "frontier-closure" in skill.get("tags", [])]
        self.assertGreaterEqual(len(generated), 100)
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            for skill in generated:
                with self.subTest(skill=skill["slug"]):
                    out_path = tmp / f"{skill['slug']}.json"
                    subprocess.run(
                        [
                            "python3",
                            f"{skill['path']}/scripts/run_frontier_starter.py",
                            "--out",
                            str(out_path),
                        ],
                        cwd=ROOT,
                        check=True,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    payload = json.loads(out_path.read_text(encoding="utf-8"))
                    self.assertEqual(payload["skill_slug"], skill["slug"])
                    self.assertEqual(payload["leaf_slug"], skill["topic_path"][1])
                    self.assertEqual(payload["source_resource_ids"], skill["source_resource_ids"])


if __name__ == "__main__":
    unittest.main()
