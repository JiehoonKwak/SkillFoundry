from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "scipy-statistical-testing-starter" / "scripts" / "run_scipy_statistical_testing.py"
EXAMPLE = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "scipy-statistical-testing-starter" / "examples" / "toy_groups.tsv"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"

SPEC = importlib.util.spec_from_file_location("scipy_testing", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class SciPyStatisticalTestingTests(unittest.TestCase):
    def test_load_groups(self) -> None:
        groups = MODULE.load_groups(EXAMPLE)
        self.assertEqual(sorted(groups), ["control", "treated"])
        self.assertEqual(len(groups["control"]), 4)

    def test_run_tests(self) -> None:
        payload = MODULE.run_tests(EXAMPLE)
        self.assertEqual(payload["higher_mean_group"], "treated")
        self.assertLess(payload["welch_ttest_pvalue"], 0.05)
        self.assertGreater(payload["cohen_d"], 0)

    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(STATS_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["higher_mean_group"], "treated")


if __name__ == "__main__":
    unittest.main()
