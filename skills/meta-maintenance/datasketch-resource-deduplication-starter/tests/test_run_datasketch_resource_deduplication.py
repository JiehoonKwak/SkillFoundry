from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "meta-maintenance" / "datasketch-resource-deduplication-starter" / "scripts" / "run_datasketch_resource_deduplication.py"
EXAMPLE = ROOT / "skills" / "meta-maintenance" / "datasketch-resource-deduplication-starter" / "examples" / "toy_resources.jsonl"
MAINT_PYTHON = ROOT / "slurm" / "envs" / "maintenance" / "bin" / "python"


class DatasketchResourceDeduplicationTests(unittest.TestCase):
    def test_script_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [str(MAINT_PYTHON), str(SCRIPT), "--input", str(EXAMPLE), "--threshold", "0.5", "--out", str(out_path)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["resource_count"], 3)
            self.assertEqual(len(payload["candidate_pairs"]), 1)


if __name__ == "__main__":
    unittest.main()
