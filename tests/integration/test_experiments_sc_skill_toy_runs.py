from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_ROOT = ROOT / "experiments" / "sc_skills"
MANIFEST_PATH = EXPERIMENT_ROOT / "batch_design_manifest.json"


class ScSkillToyRunIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_all_sc_skill_wrappers_run_on_toy_inputs(self) -> None:
        for task in self.manifest["tasks"]:
            skill_dir = ROOT / task["target_dir"]
            script_name = f"run_{task['task_slug'].replace('-', '_')}.py"
            script_path = skill_dir / "scripts" / script_name
            input_path = skill_dir / "examples" / "toy_input.json"
            with self.subTest(task=task["task_slug"]), tempfile.TemporaryDirectory() as tmpdir:
                completed = subprocess.run(
                    ["python3", str(script_path), "--input", str(input_path), "--outdir", tmpdir],
                    cwd=ROOT,
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                for deliverable in task["deliverables"]:
                    self.assertTrue((Path(tmpdir) / deliverable).exists(), f"Missing output {deliverable} for {task['task_slug']}")


if __name__ == "__main__":
    unittest.main()
