from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]


def run_json_command(*args: str) -> dict:
    completed = subprocess.run(
        ["python3", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=180,
    )
    return json.loads(completed.stdout)


class MultiSourceLiteratureBriefTests(unittest.TestCase):
    def test_query_returns_all_source_blocks(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/multi-source-literature-brief/scripts/build_literature_brief.py",
            "--query",
            "single-cell RNA-seq",
            "--limit",
            "1",
        )
        self.assertEqual(payload["limit_per_source"], 1)
        self.assertEqual(sorted(payload["sources"].keys()), ["crossref", "europepmc", "openalex", "pubmed"])
        self.assertEqual(len(payload["sources"]["openalex"]), 1)
        self.assertEqual(len(payload["sources"]["europepmc"]), 1)
        self.assertEqual(len(payload["sources"]["crossref"]), 1)
        self.assertEqual(len(payload["sources"]["pubmed"]), 1)

    def test_cli_writes_json_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "brief.json"
            subprocess.run(
                [
                    "python3",
                    "skills/scientific-knowledge/multi-source-literature-brief/scripts/build_literature_brief.py",
                    "--query",
                    "single-cell RNA-seq",
                    "--limit",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertIn("overlap_summary", payload)
            self.assertEqual(payload["query"], "single-cell RNA-seq")


if __name__ == "__main__":
    unittest.main()
