from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "agents" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-agents-and-automation"
    / "inspect-evaluation-harness-starter"
    / "scripts"
    / "run_inspect_evaluation_harness.py"
)
EXAMPLES = (
    ROOT
    / "skills"
    / "scientific-agents-and-automation"
    / "inspect-evaluation-harness-starter"
    / "examples"
)


class InspectEvaluationHarnessTests(unittest.TestCase):
    def test_candidate_beats_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            out_path = tmp_root / "summary.json"
            log_dir = tmp_root / "logs"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--cases",
                    str(EXAMPLES / "toy_eval_cases.json"),
                    "--summary-out",
                    str(out_path),
                    "--log-dir",
                    str(log_dir),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["case_count"], 4)
            self.assertAlmostEqual(payload["candidate_accuracy"], 1.0, places=6)
            self.assertAlmostEqual(payload["baseline_accuracy"], 0.25, places=6)
            self.assertGreater(payload["accuracy_gain"], 0.0)
            self.assertEqual(payload["candidate"]["log_file_count"], 1)
            self.assertEqual(payload["baseline"]["log_file_count"], 1)
            self.assertTrue(Path(payload["candidate"]["log_dir"]).exists())
            self.assertTrue(Path(payload["baseline"]["log_dir"]).exists())

    def test_missing_cases_file_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            completed = subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--cases",
                    "missing.json",
                    "--log-dir",
                    str(Path(tmp_dir) / "logs"),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("Cases file not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
