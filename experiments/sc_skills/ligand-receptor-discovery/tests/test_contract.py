from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT_PATH = SKILL_DIR / "scripts" / "run_ligand_receptor_discovery.py"
SPEC = importlib.util.spec_from_file_location("ligand_receptor_runner", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
validate_input = MODULE.validate_input


class SkillContractMetadataTests(unittest.TestCase):
    def test_toy_input_stays_raw_and_matches_expected_invariants(self) -> None:
        payload = json.loads((SKILL_DIR / "examples" / "toy_input.json").read_text(encoding="utf-8"))
        self.assertNotIn("deliverables", payload)
        self.assertNotIn("markdown", payload)
        validate_input(payload)
        self.assertEqual(len(payload["cells"]), payload["expected_invariants"]["cell_count"])
        self.assertEqual(len({cell["group"] for cell in payload["cells"]}), payload["expected_invariants"]["group_count"])
        self.assertEqual(len(payload["ligand_receptor_catalog"]), payload["expected_invariants"]["pair_count"])

    def test_metadata_declares_deliverables_and_qc_artifacts(self) -> None:
        metadata = json.loads((SKILL_DIR / "metadata.yaml").read_text(encoding="utf-8"))
        artifact_paths = [item["path"] for item in metadata["deliverables"] + metadata.get("qc_artifacts", [])]
        self.assertEqual(len(artifact_paths), len(set(artifact_paths)))
        self.assertIn("candidate_pairs.tsv", artifact_paths)
        self.assertIn("qc_pair_scores.tsv", artifact_paths)
        self.assertIn("run_summary.json", artifact_paths)
        self.assertIn("computed", metadata["starter_scope"])


if __name__ == "__main__":
    unittest.main()
