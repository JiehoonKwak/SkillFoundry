from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "meta-maintenance" / "rapidfuzz-skill-deduplication-starter" / "scripts" / "run_rapidfuzz_skill_deduplication.py"
EXAMPLE = ROOT / "skills" / "meta-maintenance" / "rapidfuzz-skill-deduplication-starter" / "examples" / "toy_skills.tsv"
MAINT_PYTHON = ROOT / "slurm" / "envs" / "maintenance" / "bin" / "python"


class RapidfuzzSkillDeduplicationTests(unittest.TestCase):
    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(MAINT_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--threshold", "85", "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["skill_count"], 4)
            self.assertGreaterEqual(len(payload["candidate_pairs"]), 2)
            self.assertGreaterEqual(payload["candidate_pairs"][0]["score"], 85)


if __name__ == "__main__":
    unittest.main()
