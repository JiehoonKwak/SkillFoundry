from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_semantic_scholar_paper_triage.py"
EXAMPLE = ROOT / "examples" / "candidate_papers.json"


class SemanticScholarPaperTriageTests(unittest.TestCase):
    def test_triage_ranks_atlas_papers_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "triage.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    str(EXAMPLE),
                    "--query",
                    "single-cell RNA-seq atlas integration",
                    "--out",
                    str(out_path),
                ],
                check=True,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["candidate_count"], 5)
            self.assertEqual(payload["top_candidates"][0]["paper_id"], "P1")

    def test_invalid_limit_fails(self) -> None:
        result = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                str(EXAMPLE),
                "--query",
                "single-cell",
                "--limit",
                "0",
            ],
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--limit must be >= 1", result.stderr or result.stdout)


if __name__ == "__main__":
    unittest.main()
