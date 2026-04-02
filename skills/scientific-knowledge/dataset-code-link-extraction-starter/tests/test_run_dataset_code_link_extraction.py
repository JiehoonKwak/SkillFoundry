from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_dataset_code_link_extraction.py"
EXAMPLE = ROOT / "examples" / "paper_text.md"


class DatasetCodeLinkExtractionTests(unittest.TestCase):
    def test_links_are_classified(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "links.json"
            subprocess.run(["python3", str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)], check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["url_count"], 4)
            self.assertEqual(len(payload["dataset_links"]), 2)
            self.assertEqual(len(payload["code_links"]), 1)
            self.assertEqual(len(payload["package_links"]), 1)


if __name__ == "__main__":
    unittest.main()
