from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run_json_command(*args: str) -> dict:
    completed = subprocess.run(
        ["python3", *args],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return json.loads(completed.stdout)


class SecondWaveApiSkillSmokeTests(unittest.TestCase):
    def test_crossref_search_returns_items(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/crossref-metadata-search/scripts/search_crossref.py",
            "--query-title",
            "single-cell RNA-seq",
            "--rows",
            "1",
        )
        self.assertEqual(payload["message-type"], "work-list")
        self.assertGreaterEqual(len(payload["message"]["items"]), 1)

    def test_rcsb_search_returns_identifiers(self) -> None:
        payload = run_json_command(
            "skills/structural-biology/rcsb-pdb-search/scripts/search_rcsb.py",
            "--query",
            "hemoglobin",
            "--rows",
            "1",
        )
        self.assertIn("result_set", payload)
        self.assertGreaterEqual(len(payload["result_set"]), 1)

    def test_rcsb_entry_summary_returns_metadata(self) -> None:
        payload = run_json_command(
            "skills/structural-biology/rcsb-pdb-entry-summary/scripts/fetch_pdb_entry.py",
            "--entry-id",
            "4HHB",
        )
        self.assertIn("struct", payload)
        self.assertIn("audit_author", payload)

    def test_ensembl_lookup_returns_gene(self) -> None:
        payload = run_json_command(
            "skills/genomics/ensembl-gene-lookup/scripts/lookup_gene.py",
            "--symbol",
            "BRCA2",
            "--species",
            "homo_sapiens",
        )
        self.assertEqual(payload["lookup"]["display_name"], "BRCA2")
        self.assertTrue(payload["lookup"]["id"].startswith("ENSG"))
        self.assertGreaterEqual(len(payload["xrefs"]), 1)

    def test_ncbi_gene_search_returns_human_gene(self) -> None:
        payload = run_json_command(
            "skills/genomics/ncbi-gene-search/scripts/search_ncbi_gene.py",
            "--symbol",
            "BRCA1",
            "--species",
            "homo sapiens",
            "--retmax",
            "1",
        )
        self.assertEqual(payload["search"]["ids"], ["672"])
        self.assertEqual(payload["genes"][0]["symbol"], "BRCA1")
        self.assertEqual(payload["genes"][0]["organism"], "Homo sapiens")

    def test_pubmed_search_returns_ids(self) -> None:
        payload = run_json_command(
            "skills/scientific-knowledge/ncbi-pubmed-search/scripts/search_pubmed.py",
            "--term",
            "single-cell RNA-seq",
            "--retmax",
            "1",
        )
        self.assertIn("search", payload)
        self.assertIn("summary", payload)
        self.assertGreaterEqual(len(payload["search"]["esearchresult"]["idlist"]), 1)

    def test_clinicaltrials_search_returns_studies(self) -> None:
        payload = run_json_command(
            "skills/clinical-biomedical-data-science/clinicaltrials-study-search/scripts/search_clinicaltrials.py",
            "--condition",
            "melanoma",
            "--page-size",
            "1",
        )
        self.assertIn("studies", payload)
        self.assertGreaterEqual(len(payload["studies"]), 1)


if __name__ == "__main__":
    unittest.main()
