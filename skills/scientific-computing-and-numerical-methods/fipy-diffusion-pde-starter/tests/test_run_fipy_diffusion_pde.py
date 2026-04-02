from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "numerics" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "scientific-computing-and-numerical-methods"
    / "fipy-diffusion-pde-starter"
    / "scripts"
    / "run_fipy_diffusion_pde.py"
)


class FipyDiffusionPdeTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_profile(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "fipy.json"
            subprocess.run(
                [str(PYTHON), str(SCRIPT), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["cell_count"], 20)
            self.assertEqual(payload["initial_mass"], 5.0)
            self.assertEqual(payload["final_mass"], 5.0)
            self.assertAlmostEqual(payload["center_value"], 0.000438, places=6)
            self.assertAlmostEqual(payload["leading_profile"][0], 0.997403, places=6)


if __name__ == "__main__":
    unittest.main()
