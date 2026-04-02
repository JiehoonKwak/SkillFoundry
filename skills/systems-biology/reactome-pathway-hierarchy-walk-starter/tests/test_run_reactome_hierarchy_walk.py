from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = ROOT / "skills" / "systems-biology" / "reactome-pathway-hierarchy-walk-starter" / "scripts" / "run_reactome_hierarchy_walk.py"
SPEC = importlib.util.spec_from_file_location("reactome_hierarchy_walk", SCRIPT_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ReactomeHierarchyWalkTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_canonical_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "hierarchy.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT_PATH),
                    "--species",
                    "9606",
                    "--stable-id",
                    "R-HSA-141409",
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
            self.assertEqual(payload["stable_id"], "R-HSA-141409")
            self.assertEqual(payload["top_level_pathway"], "Cell Cycle")
            self.assertEqual(payload["ancestor_count"], 5)
            self.assertEqual(payload["direct_children_count"], 0)

    def test_asset_fallback_is_available_for_canonical_target(self) -> None:
        payload = MODULE.load_asset_fallback("9606", "R-HSA-141409", "forced test fallback")
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["result_origin"], "asset_fallback")
        self.assertEqual(payload["fallback_reason"], "forced test fallback")
        self.assertEqual(payload["stable_id"], "R-HSA-141409")

    def test_rejects_unknown_stable_id(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT_PATH),
                "--species",
                "9606",
                "--stable-id",
                "R-HSA-NOTREAL",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        message = completed.stderr + completed.stdout
        self.assertTrue(
            "was not found" in message or "Reactome hierarchy request failed" in message,
            msg=message,
        )


if __name__ == "__main__":
    unittest.main()
