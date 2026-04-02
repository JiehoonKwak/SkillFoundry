from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "validate_sc_skill_experiments.py"


class ExperimentScSkillsToySmokeTests(unittest.TestCase):
    def test_validator_script_runs_on_subset(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            json_out = tmpdir_path / "report.json"
            markdown_out = tmpdir_path / "report.md"
            completed = subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--all",
                    "--limit",
                    "3",
                    "--run-local-tests",
                    "--scratch-root",
                    str(tmpdir_path / "artifacts"),
                    "--json-out",
                    str(json_out),
                    "--markdown-out",
                    str(markdown_out),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(json_out.read_text(encoding="utf-8"))
            self.assertEqual(payload["total_tasks"], 3)
            self.assertEqual(payload["passed_tasks"], 3)
            self.assertEqual(payload["failed_tasks"], [])
            self.assertTrue(markdown_out.exists())


if __name__ == "__main__":
    unittest.main()
