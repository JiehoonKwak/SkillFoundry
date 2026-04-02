from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REF_ROOT = ROOT / "ref" / "skill"
EXPERIMENT_ROOT = ROOT / "experiments" / "sc_skills"
MANIFEST_PATH = EXPERIMENT_ROOT / "batch_design_manifest.json"


class ScSkillExperimentContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        cls.tasks = cls.manifest["tasks"]

    def test_manifest_covers_all_reference_skill_files(self) -> None:
        manifest_sources = {task["source_reference"] for task in self.tasks}
        repo_sources = {str(path.relative_to(ROOT)) for path in sorted(REF_ROOT.glob("*.md"))}
        self.assertEqual(manifest_sources, repo_sources)

    def test_manifest_targets_sc_skill_root(self) -> None:
        for task in self.tasks:
            self.assertTrue(task["target_dir"].startswith("experiments/sc_skills/"), task["target_dir"])

    def test_materialized_packages_have_standard_experiment_scaffold(self) -> None:
        for task in self.tasks:
            target_dir = ROOT / task["target_dir"]
            self.assertTrue(target_dir.is_dir(), f"Missing experiment dir: {target_dir}")
            required_files = {
                target_dir / "SKILL.md",
                target_dir / "refs.md",
                target_dir / "design_prompt.md",
                target_dir / "source_reference.md",
                target_dir / "metadata.yaml",
                target_dir / "examples" / "toy_input.json",
                target_dir / "assets" / "README.md",
                target_dir / "scripts" / f"run_{task['task_slug'].replace('-', '_')}.py",
                target_dir / "scripts" / "run_exercise.py",
                target_dir / "scripts" / "validate_outputs.py",
                target_dir / "tests" / f"test_{task['task_slug'].replace('-', '_')}.py",
                target_dir / "tests" / "test_contract.py",
            }
            for path in required_files:
                self.assertTrue(path.exists(), f"Missing scaffold path: {path}")

    def test_source_reference_copy_matches_original(self) -> None:
        for task in self.tasks:
            source_text = (ROOT / task["source_reference"]).read_text(encoding="utf-8")
            copied_text = (ROOT / task["target_dir"] / "source_reference.md").read_text(encoding="utf-8")
            self.assertEqual(copied_text, source_text)

    def test_skill_mentions_package_scaffold_and_toy_validation(self) -> None:
        required_phrases = [
            "Compatible with Codex, Claude Code",
            "## Package scaffold",
            "## Toy validation",
            "## Starter scope",
            "scripts/",
            "examples/toy_input.json",
        ]
        for task in self.tasks:
            skill_text = (ROOT / task["target_dir"] / "SKILL.md").read_text(encoding="utf-8")
            for phrase in required_phrases:
                self.assertIn(phrase, skill_text, f"{phrase!r} missing in {task['task_slug']}")

    def test_design_prompt_scopes_framework_and_scaffold_expectations(self) -> None:
        for task in self.tasks:
            prompt_text = (ROOT / task["target_dir"] / "design_prompt.md").read_text(encoding="utf-8")
            self.assertIn(task["target_dir"], prompt_text)
            self.assertIn("Do not edit `registry/`, `skills/`, `site/`, `README.md`, `experiments.md`, or planning files", prompt_text)
            self.assertIn("metadata.yaml", prompt_text)
            self.assertIn("examples/", prompt_text)
            self.assertIn("tests/", prompt_text)
            self.assertIn("Do not use `examples/toy_input.json` to store final deliverables", prompt_text)
            self.assertIn("Replace the shared contract-only runner with a task-specific runner", prompt_text)
            self.assertIn("Implement 2 to 4 method-shaped computations", prompt_text)
            self.assertIn("Starter scope", prompt_text)

    def test_batch_helper_print_only_mode(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                "scripts/batch_design_experiment_skills.py",
                "--task",
                "annotation",
                "--print-only",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        selected = json.loads(completed.stdout)
        self.assertEqual(selected, ["annotation"])


if __name__ == "__main__":
    unittest.main()
