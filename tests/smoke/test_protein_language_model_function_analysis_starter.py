from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class ProteinLanguageModelFunctionAnalysisSmokeTests(unittest.TestCase):
    def test_skill_smoke_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            embeddings_out = tmp / "embeddings.tsv"
            summary_out = tmp / "summary.json"
            subprocess.run(
                [
                    "python3",
                    "skills/proteomics/protein-language-model-function-analysis-starter/scripts/run_protein_language_model_function_analysis.py",
                    "--input",
                    "skills/proteomics/protein-language-model-function-analysis-starter/examples/toy_sequences.fasta",
                    "--labels",
                    "skills/proteomics/protein-language-model-function-analysis-starter/examples/toy_labels.tsv",
                    "--config",
                    "skills/proteomics/protein-language-model-function-analysis-starter/examples/analysis_config.json",
                    "--embeddings-out",
                    str(embeddings_out),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["label_summary"]["labels_seen"], ["dna_binding_like", "enzyme_like", "membrane_like"])
            self.assertEqual(payload["holdout_accuracy"], 1.0)


if __name__ == "__main__":
    unittest.main()
