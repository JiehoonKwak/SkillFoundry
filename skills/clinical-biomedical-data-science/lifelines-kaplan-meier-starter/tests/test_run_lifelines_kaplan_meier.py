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
    / "clinical-biomedical-data-science"
    / "lifelines-kaplan-meier-starter"
    / "scripts"
    / "run_lifelines_kaplan_meier.py"
)
INPUT = (
    ROOT
    / "skills"
    / "clinical-biomedical-data-science"
    / "lifelines-kaplan-meier-starter"
    / "examples"
    / "toy_survival_cohort.tsv"
)


class LifelinesKaplanMeierTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_medians(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_out = Path(tmp_dir) / "summary.json"
            png_out = Path(tmp_dir) / "plot.png"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--input",
                    str(INPUT),
                    "--summary-out",
                    str(summary_out),
                    "--png-out",
                    str(png_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["median_survival_by_group"]["A"], 5.0)
            self.assertEqual(payload["median_survival_by_group"]["B"], 4.0)
            self.assertEqual(payload["survival_at_times"]["A"]["6"], 0.25)
            self.assertEqual(payload["survival_at_times"]["B"]["4"], 0.5)
            self.assertTrue(png_out.exists())

    def test_rejects_missing_input(self) -> None:
        completed = subprocess.run(
            [
                str(PYTHON),
                str(SCRIPT),
                "--input",
                "missing.tsv",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Input cohort not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
