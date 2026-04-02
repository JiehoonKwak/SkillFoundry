from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = (
    ROOT
    / "skills"
    / "systems-biology"
    / "string-interaction-partners-starter"
    / "scripts"
    / "run_string_interaction_partners.py"
)
EXAMPLE = (
    ROOT
    / "skills"
    / "systems-biology"
    / "string-interaction-partners-starter"
    / "examples"
    / "tp53_query.txt"
)
ASSET = (
    ROOT
    / "skills"
    / "systems-biology"
    / "string-interaction-partners-starter"
    / "assets"
    / "tp53_interaction_partners.json"
)


class StringInteractionPartnersTests(unittest.TestCase):
    def test_runtime_tp53_query_returns_partner_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "tp53_partners.json"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--identifier-file",
                    str(EXAMPLE),
                    "--species",
                    "9606",
                    "--limit",
                    "5",
                    "--required-score",
                    "700",
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
            self.assertEqual(payload["query_preferred_name"], "TP53")
            self.assertEqual(payload["partner_count"], 5)
            self.assertIn(payload["result_origin"], {"live_api", "asset_fallback"})
            partner_names = [item["partner_preferred_name"] for item in payload["top_partners"]]
            self.assertIn("SFN", partner_names)
            self.assertGreaterEqual(payload["score_max"], payload["score_min"])

    def test_asset_matches_shape(self) -> None:
        payload = json.loads(ASSET.read_text(encoding="utf-8"))
        self.assertEqual(payload["query_preferred_name"], "TP53")
        self.assertEqual(payload["partner_count"], 5)
        self.assertEqual(len(payload["top_partners"]), 5)

    def test_empty_identifier_file_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_query = Path(tmp_dir) / "empty.txt"
            bad_query.write_text("\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--identifier-file",
                    str(bad_query),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("must contain at least one non-empty line", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
