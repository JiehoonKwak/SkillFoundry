from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "proteomics" / "interpro-entry-summary" / "scripts" / "fetch_interpro_entry.py"


class InterProEntrySummaryTests(unittest.TestCase):
    def test_runtime_summary_for_known_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "interpro.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--accession",
                    "IPR000023",
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["accession"], "IPR000023")
            self.assertEqual(payload["type"], "domain")
            self.assertGreaterEqual(payload["go_term_count"], 2)
            self.assertIn("Phosphofructokinase", payload["name"])

    def test_rejects_invalid_accession(self) -> None:
        completed = subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--accession",
                "IPR999999999",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("InterPro accession not found", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
