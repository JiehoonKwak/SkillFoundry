from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "transcriptomics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "pydeseq2-differential-expression-starter"
    / "scripts"
    / "run_pydeseq2_differential_expression.py"
)
COUNTS = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "pydeseq2-differential-expression-starter"
    / "examples"
    / "toy_counts.tsv"
)
METADATA = (
    ROOT
    / "skills"
    / "transcriptomics"
    / "pydeseq2-differential-expression-starter"
    / "examples"
    / "toy_metadata.tsv"
)


class PyDESeq2DifferentialExpressionTests(unittest.TestCase):
    def test_runtime_summary_contains_expected_ranked_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    str(COUNTS),
                    "--metadata",
                    str(METADATA),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=300,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_count"], 6)
            self.assertEqual(payload["gene_count"], 4)
            self.assertEqual(payload["top_gene"], "GeneD")
            self.assertEqual(payload["significant_gene_count"], 4)
            self.assertLess(payload["top_adjusted_p_value"], 1e-20)

    def test_rejects_missing_metadata_sample(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_metadata = Path(tmp_dir) / "bad_metadata.tsv"
            bad_metadata.write_text("sample\tcondition\nS1\tA\nS2\tA\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--counts",
                    str(COUNTS),
                    "--metadata",
                    str(bad_metadata),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=300,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Metadata is missing samples", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
