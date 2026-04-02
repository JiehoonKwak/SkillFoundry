from __future__ import annotations

import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "systems-biology" / "quickgo-term-search" / "scripts" / "search_quickgo_terms.py"


def load_script_module():
    spec = importlib.util.spec_from_file_location("quickgo_term_search", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = load_script_module()


class QuickGOTermSearchTests(unittest.TestCase):
    def test_normalize_positive_int_rejects_non_positive_values(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.normalize_positive_int(0, "limit")

    def test_compact_term_record_extracts_expected_fields(self) -> None:
        compact = MODULE.compact_term_record(
            {
                "id": "GO:0006915",
                "name": "apoptotic process",
                "aspect": "biological_process",
                "isObsolete": False,
                "definition": {"text": "Programmed cell death."},
            }
        )
        self.assertEqual(compact["id"], "GO:0006915")
        self.assertEqual(compact["name"], "apoptotic process")
        self.assertEqual(compact["aspect"], "biological_process")
        self.assertFalse(compact["is_obsolete"])
        self.assertEqual(compact["definition"], "Programmed cell death.")

    def test_live_quickgo_query_returns_a_go_term(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "quickgo.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--query",
                    "apoptosis",
                    "--limit",
                    "1",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            payload = MODULE.json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["query"], "apoptosis")
            self.assertEqual(payload["term_count"], 1)
            self.assertGreaterEqual(payload["number_of_hits"], 1)
            self.assertTrue(payload["terms"][0]["id"].startswith("GO:"))


if __name__ == "__main__":
    unittest.main()
