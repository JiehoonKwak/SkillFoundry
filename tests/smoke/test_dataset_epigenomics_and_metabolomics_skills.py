from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run_json_command(*args: str) -> dict:
    completed = subprocess.run(
        ["python3", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=90,
    )
    return json.loads(completed.stdout)


class DatasetEpigenomicsAndMetabolomicsSkillSmokeTests(unittest.TestCase):
    def test_encode_experiment_search_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/epigenomics-and-chromatin/encode-experiment-search/scripts/search_encode_experiments.py",
            "--search-term",
            "ATAC-seq",
            "--limit",
            "1",
        )
        self.assertEqual(payload["search_term"], "ATAC-seq")
        self.assertEqual(payload["experiment_count"], 1)
        self.assertGreaterEqual(payload["total"], 1)
        self.assertIn("ATAC-seq", payload["experiments"][0]["assay_title"])

    def test_biosamples_sample_search_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/data-acquisition-and-dataset-handling/biosamples-sample-search/scripts/search_biosamples_samples.py",
            "--text",
            "breast cancer",
            "--size",
            "1",
        )
        self.assertEqual(payload["text"], "breast cancer")
        self.assertEqual(payload["sample_count"], 1)
        self.assertGreaterEqual(payload["total_elements"], 1)
        self.assertTrue(payload["samples"][0]["accession"].startswith("SAME"))

    def test_metabolights_study_search_returns_expected_fields(self) -> None:
        payload = run_json_command(
            "skills/metabolomics-and-other-omics/metabolights-study-search/scripts/search_metabolights_studies.py",
            "--query",
            "diabetes",
            "--rows",
            "1",
        )
        self.assertEqual(payload["query"], "diabetes")
        self.assertEqual(payload["study_count"], 1)
        self.assertGreaterEqual(payload["total_studies"], 1)
        self.assertTrue(payload["studies"][0]["accession"].startswith("MTBLS"))
        self.assertTrue(payload["studies"][0]["title"])


if __name__ == "__main__":
    unittest.main()
