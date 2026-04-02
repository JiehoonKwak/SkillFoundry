from __future__ import annotations

import importlib.util
import json
import tempfile
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MODULE_PATH = ROOT / "skills" / "genomics" / "ncbi-gene-search" / "scripts" / "search_ncbi_gene.py"
BRCA1_ASSET = ROOT / "skills" / "genomics" / "ncbi-gene-search" / "assets" / "brca1_gene_summary.json"


def load_module():
    spec = importlib.util.spec_from_file_location("ncbi_gene_search_tests", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def raw_gene_record(compact_gene: dict) -> dict:
    return {
        "uid": compact_gene.get("gene_id"),
        "name": compact_gene.get("symbol"),
        "description": compact_gene.get("description"),
        "nomenclaturename": compact_gene.get("official_name"),
        "nomenclaturesymbol": compact_gene.get("official_symbol"),
        "organism": {
            "scientificname": compact_gene.get("organism"),
            "commonname": compact_gene.get("common_name"),
            "taxid": compact_gene.get("taxid"),
        },
        "chromosome": compact_gene.get("chromosome"),
        "maplocation": compact_gene.get("map_location"),
        "otheraliases": ", ".join(compact_gene.get("aliases", [])),
        "mim": compact_gene.get("mim", []),
        "genomicinfo": [compact_gene.get("genomic_info", {})],
    }


def fixture_payload(symbol: str) -> tuple[dict, dict]:
    if symbol == "BRCA1":
        compact_gene = json.loads(BRCA1_ASSET.read_text(encoding="utf-8"))["genes"][0]
    elif symbol == "TP53":
        compact_gene = {
            "gene_id": "7157",
            "symbol": "TP53",
            "description": "tumor protein p53",
            "official_name": "tumor protein p53",
            "official_symbol": "TP53",
            "organism": "Homo sapiens",
            "common_name": "human",
            "taxid": 9606,
            "chromosome": "17",
            "map_location": "17p13.1",
            "aliases": ["P53", "BCC7", "LFS1"],
            "mim": [],
            "summary": "Synthetic local fixture used to exercise CLI output writing.",
            "genomic_info": {
                "chraccver": "NC_000017.11",
                "chrstart": 7661779,
                "chrstop": 7687550,
                "exoncount": 11,
            },
        }
    else:
        raise ValueError(symbol)

    search_ids = [str(compact_gene["gene_id"])]
    search_payload = {
        "esearchresult": {
            "count": "1",
            "idlist": search_ids,
            "querytranslation": f"{symbol}[sym] AND \"Homo sapiens\"[Organism]",
        }
    }
    summary_payload = {
        "result": {
            "uids": search_ids,
            str(compact_gene["gene_id"]): raw_gene_record(compact_gene),
        }
    }
    return search_payload, summary_payload


class NcbiGeneSearchTests(unittest.TestCase):
    def test_brca1_query_returns_human_gene(self) -> None:
        module = load_module()
        search_payload, summary_payload = fixture_payload("BRCA1")

        def _search_gene(query_symbol: str, species: str, retmax: int, email: str | None):
            query = module.build_query(query_symbol, species)
            return query, "mock://search", search_payload

        def _summarize_gene_ids(ids: list[str], email: str | None):
            return "mock://summary", summary_payload

        module.search_gene = _search_gene
        module.summarize_gene_ids = _summarize_gene_ids
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "brca1.json"
            argv = ["search_ncbi_gene.py", "--symbol", "BRCA1", "--species", "homo sapiens", "--retmax", "1", "--out", str(out_path)]
            old_argv = sys.argv
            sys.argv = argv
            try:
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = old_argv
            payload = json.loads(out_path.read_text(encoding="utf-8"))
        self.assertEqual(payload["search"]["ids"], ["672"])
        self.assertEqual(payload["genes"][0]["symbol"], "BRCA1")
        self.assertEqual(payload["genes"][0]["organism"], "Homo sapiens")

    def test_cli_writes_json_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "gene.json"
            module = load_module()
            search_payload, summary_payload = fixture_payload("TP53")

            def _search_gene(query_symbol: str, species: str, retmax: int, email: str | None):
                query = module.build_query(query_symbol, species)
                return query, "mock://search", search_payload

            def _summarize_gene_ids(ids: list[str], email: str | None):
                return "mock://summary", summary_payload

            module.search_gene = _search_gene
            module.summarize_gene_ids = _summarize_gene_ids
            argv = ["search_ncbi_gene.py", "--symbol", "TP53", "--species", "homo sapiens", "--retmax", "1", "--out", str(out_path)]
            old_argv = sys.argv
            sys.argv = argv
            try:
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = old_argv
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["genes"][0]["symbol"], "TP53")
            self.assertEqual(payload["genes"][0]["organism"], "Homo sapiens")


if __name__ == "__main__":
    unittest.main()
