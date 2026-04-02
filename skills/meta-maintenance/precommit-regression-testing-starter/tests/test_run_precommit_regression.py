from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "meta-maintenance"
    / "precommit-regression-testing-starter"
    / "scripts"
    / "run_precommit_regression.py"
)


class PrecommitRegressionTestingTests(unittest.TestCase):
    def test_runtime_summary_reports_successful_hooks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workspace = Path(tmp_dir) / "workspace"
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--workspace",
                    str(workspace),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["returncode"], 0)
            self.assertEqual(payload["tracked_file_count"], 3)
            self.assertEqual(payload["hook_ids"], ["python-compile", "json-validate"])


if __name__ == "__main__":
    unittest.main()
