from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "statistical-and-machine-learning-foundations-for-science"
    / "arviz-posterior-diagnostics-starter"
    / "scripts"
    / "run_arviz_posterior_diagnostics.py"
)


class ArvizPosteriorDiagnosticsTests(unittest.TestCase):
    def test_runtime_summary_reports_good_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["parameter_count"], 3)
            self.assertLess(payload["max_rhat"], 1.02)
            self.assertGreater(payload["min_ess_bulk"], 500)

    def test_rejects_single_chain(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--chains", "1"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("chains must be at least 2", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
