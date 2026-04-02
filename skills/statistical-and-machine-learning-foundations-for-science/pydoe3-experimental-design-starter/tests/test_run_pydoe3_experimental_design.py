from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "pydoe3-experimental-design-starter" / "scripts" / "run_pydoe3_experimental_design.py"
EXAMPLE = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "pydoe3-experimental-design-starter" / "examples" / "toy_factors.json"
STATS_PYTHON = ROOT / "slurm" / "envs" / "statistics" / "bin" / "python"


class PyDoe3ExperimentalDesignTests(unittest.TestCase):
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
            self.assertEqual(payload["row_count"], 4)
            self.assertEqual(payload["factor_names"], ["temperature_c", "ph"])
            self.assertEqual(payload["rows"][0]["temperature_c"], 20.0)


if __name__ == "__main__":
    unittest.main()
