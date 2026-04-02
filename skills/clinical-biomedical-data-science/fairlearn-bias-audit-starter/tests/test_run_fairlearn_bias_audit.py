from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
SKILL_DIR = ROOT / "skills" / "clinical-biomedical-data-science" / "fairlearn-bias-audit-starter"
SCRIPT = SKILL_DIR / "scripts" / "run_fairlearn_bias_audit.py"
INPUT = SKILL_DIR / "examples" / "toy_fairness_cohort.tsv"
ASSET = SKILL_DIR / "assets" / "toy_fairlearn_bias_audit_summary.json"


class FairlearnBiasAuditTests(unittest.TestCase):
    def test_runtime_summary_matches_verified_asset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            expected = json.loads(ASSET.read_text(encoding="utf-8"))
            self.assertEqual(payload, expected)

    def test_relaxed_thresholds_clear_flags(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_out = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--demographic-parity-threshold",
                    "0.6",
                    "--equalized-odds-threshold",
                    "0.6",
                    "--minimum-ratio-threshold",
                    "0.0",
                    "--summary-out",
                    str(summary_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["flags"], [])
            self.assertIn("No configured fairness thresholds", payload["interpretation"])

    def test_rejects_missing_required_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_input = Path(tmp_dir) / "bad.tsv"
            bad_input.write_text(
                "patient_id\tgroup\tprediction\nP1\tA\t1\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(bad_input),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Missing required columns", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
