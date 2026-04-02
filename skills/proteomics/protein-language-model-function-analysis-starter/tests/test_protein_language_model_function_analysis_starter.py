from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "proteomics" / "protein-language-model-function-analysis-starter"
SCRIPT = SKILL / "scripts" / "run_protein_language_model_function_analysis.py"
WRAPPER = SKILL / "scripts" / "run_frontier_starter.py"


class ProteinLanguageModelFunctionAnalysisTests(unittest.TestCase):
    def run_script(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *args],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )

    def test_mock_backend_produces_expected_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            embeddings_out = tmp / "embeddings.tsv"
            summary_out = tmp / "summary.json"
            result = self.run_script(
                "--input", str(SKILL / "examples" / "toy_sequences.fasta"),
                "--labels", str(SKILL / "examples" / "toy_labels.tsv"),
                "--config", str(SKILL / "examples" / "analysis_config.json"),
                "--embeddings-out", str(embeddings_out),
                "--summary-out", str(summary_out),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["skill_slug"], "protein-language-model-function-analysis-starter")
            self.assertEqual(payload["backend"], "mock")
            self.assertEqual(payload["sequence_count"], 6)
            self.assertEqual(payload["embedding_dim"], 32)
            self.assertEqual(payload["label_summary"]["label_count"], 3)
            self.assertEqual(payload["holdout_accuracy"], 1.0)
            self.assertEqual(payload["predictions"][0]["predicted_label"], payload["predictions"][0]["observed_label"])

            with embeddings_out.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 6)
            self.assertIn("dim_031", rows[0])

    def test_wrapper_runs_bundled_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_out = Path(tmpdir) / "wrapper.json"
            completed = subprocess.run(
                [sys.executable, str(WRAPPER), "--out", str(summary_out)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=60,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["model_family"], "esm2")
            self.assertEqual(payload["nearest_neighbors"][0]["sequence_id"], "enzyme_like_a")

    def test_invalid_sequence_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = tmp / "broken.fasta"
            fasta.write_text(">bad\nMSTX*\n", encoding="utf-8")
            embeddings_out = tmp / "embeddings.tsv"
            summary_out = tmp / "summary.json"
            result = self.run_script(
                "--input", str(fasta),
                "--embeddings-out", str(embeddings_out),
                "--summary-out", str(summary_out),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unsupported residues", result.stderr)


if __name__ == "__main__":
    unittest.main()
