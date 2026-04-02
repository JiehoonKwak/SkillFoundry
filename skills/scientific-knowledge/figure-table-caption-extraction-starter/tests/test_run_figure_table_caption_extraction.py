from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_figure_table_caption_extraction.py"
EXAMPLE = ROOT / "examples" / "paper_excerpt.txt"


class FigureTableCaptionExtractionTests(unittest.TestCase):
    def test_captions_are_extracted(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "captions.json"
            subprocess.run(["python3", str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)], check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["figure_caption_count"], 2)
            self.assertEqual(payload["table_caption_count"], 2)
            self.assertIn("UMAP projection", payload["figure_captions"][0]["caption"])


if __name__ == "__main__":
    unittest.main()
