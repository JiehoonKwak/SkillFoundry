from __future__ import annotations

import importlib.util
import io
import json
import subprocess
import tempfile
import urllib.error
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = (
    ROOT
    / "skills"
    / "drug-discovery-and-cheminformatics"
    / "chembl-molecule-search"
    / "scripts"
    / "search_chembl_molecules.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("search_chembl_molecules", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class ChemblMoleculeSearchTests(unittest.TestCase):
    def test_normalize_query_strips_whitespace(self) -> None:
        module = load_module()
        self.assertEqual(module.normalize_query(" imatinib "), "imatinib")

    def test_load_cached_summary_returns_curated_imatinib_example(self) -> None:
        module = load_module()
        payload = module.load_cached_summary("imatinib", 1)
        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["molecules"][0]["chembl_id"], "CHEMBL941")

    def test_fetch_search_retries_transient_http_error(self) -> None:
        module = load_module()
        transient = urllib.error.HTTPError(
            url="https://www.ebi.ac.uk/chembl/api/data/molecule/search?q=imatinib&limit=1&format=json",
            code=500,
            msg="Server Error",
            hdrs={"Retry-After": "0"},
            fp=io.BytesIO(b"temporary server error"),
        )
        with mock.patch.object(module, "urlopen", side_effect=[transient, io.StringIO('{"molecules": []}')]) as urlopen:
            with mock.patch.object(module.time, "sleep"):
                payload = module.fetch_search("imatinib", 1, retries=2)
        self.assertEqual(payload["molecules"], [])
        self.assertEqual(urlopen.call_count, 2)

    def test_live_query_writes_expected_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "chembl.json"
            subprocess.run(
                [
                    "python3",
                    "skills/drug-discovery-and-cheminformatics/chembl-molecule-search/scripts/search_chembl_molecules.py",
                    "--query",
                    "imatinib",
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
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["query"], "imatinib")
            self.assertEqual(len(payload["molecules"]), 1)
            self.assertEqual(payload["molecules"][0]["chembl_id"], "CHEMBL941")
            self.assertEqual(payload["molecules"][0]["preferred_name"], "IMATINIB")


if __name__ == "__main__":
    unittest.main()
