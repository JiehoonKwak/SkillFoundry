from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SKILL = ROOT / "skills" / "statistical-and-machine-learning-foundations-for-science" / "foundation-models-for-biology-chemistry-science-starter"


class FrontierStarterTests(unittest.TestCase):
    def test_starter_summary_contains_leaf_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "summary.json"
            subprocess.run(
                ["python3", str(SKILL / "scripts" / "run_frontier_starter.py"), "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["skill_slug"], "foundation-models-for-biology-chemistry-science-starter")
            self.assertEqual(payload["leaf_slug"], "foundation-models-for-biology-chemistry-science")
            self.assertEqual(payload["source_resource_ids"], ["scgpt-github"])


if __name__ == "__main__":
    unittest.main()
