from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_benchmark_table_mining.py"
EXAMPLE = ROOT / "examples" / "benchmark_note.md"
SAMPLE_NOTE = ROOT / "examples" / "sample_benchmark_notes.md"
SAMPLE_ROWS = ROOT / "assets" / "sample_benchmark_rows.json"


class BenchmarkTableMiningTests(unittest.TestCase):
    def test_best_method_is_identified(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "benchmark.json"
            subprocess.run(["python3", str(SCRIPT), "--input", str(EXAMPLE), "--out", str(out_path)], check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["row_count"], 3)
            self.assertEqual(payload["best_method"]["method"], "scVI")
            self.assertEqual(payload["metric_name"], "ARI")

    def test_multi_table_sample_matches_canonical_asset(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sample_benchmark.json"
            subprocess.run(["python3", str(SCRIPT), "--input", str(SAMPLE_NOTE), "--out", str(out_path)], check=True)
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            expected = json.loads(SAMPLE_ROWS.read_text(encoding="utf-8"))
            self.assertEqual(payload["table_count"], 2)
            self.assertEqual(payload["selected_table_index"], 1)
            self.assertEqual(payload["metric_name"], "Score")
            self.assertEqual(payload["columns"], expected["columns"])
            self.assertEqual(payload["rows"], expected["rows"])
            self.assertEqual(payload["best_method"]["model"], "scBERT")


if __name__ == "__main__":
    unittest.main()
