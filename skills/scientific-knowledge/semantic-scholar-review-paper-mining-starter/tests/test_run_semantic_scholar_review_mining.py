from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_semantic_scholar_review_mining.py"
EXAMPLE = ROOT / "examples" / "paper_metadata.json"


class SemanticScholarReviewMiningTests(unittest.TestCase):
    def test_review_titles_are_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "reviews.json"
            subprocess.run(["python3", str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)], check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["review_paper_count"], 3)
            self.assertEqual(payload["review_papers"][0]["paper_id"], "R1")

    def test_invalid_limit_fails(self) -> None:
        result = subprocess.run(
            ["python3", str(SCRIPT), "--input", str(EXAMPLE), "--limit", "0"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--limit must be >= 1", result.stderr or result.stdout)


if __name__ == "__main__":
    unittest.main()
