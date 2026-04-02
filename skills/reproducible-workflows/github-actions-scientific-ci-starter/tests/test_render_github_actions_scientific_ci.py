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
    / "reproducible-workflows"
    / "github-actions-scientific-ci-starter"
    / "scripts"
    / "render_github_actions_scientific_ci.py"
)


class GithubActionsScientificCiTests(unittest.TestCase):
    def test_renders_expected_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_root = Path(tmp_dir)
            workflow_out = tmp_root / "sciskill_ci.yml"
            summary_out = tmp_root / "summary.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--workflow-out",
                    str(workflow_out),
                    "--summary-out",
                    str(summary_out),
                    "--smoke-target",
                    "smoke-zarr",
                    "--smoke-target",
                    "smoke-openmm-md",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(summary_out.read_text(encoding="utf-8"))
            workflow_text = workflow_out.read_text(encoding="utf-8")
            self.assertTrue(payload["contains_validate"])
            self.assertTrue(payload["contains_make_test"])
            self.assertEqual(payload["smoke_target_count"], 2)
            self.assertIn("actions/checkout@v4", workflow_text)
            self.assertIn("make smoke-zarr", workflow_text)
            self.assertIn("make smoke-openmm-md", workflow_text)

    def test_uses_default_smokes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workflow_out = Path(tmp_dir) / "default.yml"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--workflow-out",
                    str(workflow_out),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            workflow_text = workflow_out.read_text(encoding="utf-8")
            self.assertIn("make smoke-zarr", workflow_text)
            self.assertIn("make smoke-openmm-md", workflow_text)
            self.assertIn("make smoke-optuna", workflow_text)


if __name__ == "__main__":
    unittest.main()
