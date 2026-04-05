from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "run_ref_skill_framework_batches.py"


def load_batch_module():
    spec = importlib.util.spec_from_file_location("sc_skill_batch_runner", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class ScSkillFrameworkBatchesIntegrationTests(unittest.TestCase):
    def test_print_plan_shows_all_tasks_in_batches(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "run",
                    "--label",
                    "sc-skill-test",
                    "--state-root",
                    tmpdir,
                    "--batch-size",
                    "4",
                    "--print-plan",
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(completed.stdout)
            self.assertEqual(payload["batch_size"], 4)
            self.assertEqual(len(payload["remaining_tasks"]), 17)
            self.assertEqual(len(payload["batch_plan"]), 5)
            self.assertEqual(payload["batch_plan"][0][0], "annotation")

    def test_status_returns_missing_for_unknown_label(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "status",
                    "--label",
                    "missing-label",
                    "--state-root",
                    tmpdir,
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            payload = json.loads(completed.stdout)
            self.assertEqual(payload["status"], "missing")

    def test_framework_command_includes_explicit_model_and_effort_overrides(self) -> None:
        module = load_batch_module()
        command = module.framework_command(
            "annotation",
            "none",
            framework_model="gpt-5.4",
            framework_reasoning_effort="xhigh",
            framework_profile=None,
            framework_codex_bin="codex",
            framework_full_auto=False,
            framework_extra_context="upgrade starter",
        )
        self.assertIn("--framework-model", command)
        self.assertIn("gpt-5.4", command)
        self.assertIn("--framework-reasoning-effort", command)
        self.assertIn("xhigh", command)
        self.assertIn("--framework-extra-context", command)


if __name__ == "__main__":
    unittest.main()
