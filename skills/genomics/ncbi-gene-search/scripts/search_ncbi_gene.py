#!/usr/bin/env python3
"""Search NCBI Gene via E-utilities and return compact gene summaries."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
USER_AGENT = "SciSkillUniverse/0.2"


def fetch_json(url: str) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 3:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = float(retry_after) if retry_after else float(1 + attempt)
            time.sleep(delay)
        except URLError as exc:
            last_error = exc
            if attempt == 3:
                raise
            time.sleep(float(1 + attempt))
    raise RuntimeError(f"NCBI Gene request failed after retries: {last_error}")


def build_common_params(email: str | None) -> dict[str, str]:
    params = {"tool": "SciSkillUniverse", "retmode": "json"}
    if email:
        params["email"] = email
    return params


def build_query(symbol: str, species: str) -> str:
    return f"{symbol}[sym] AND {species}[orgn]"


def search_gene(symbol: str, species: str, retmax: int, email: str | None) -> tuple[str, str, dict[str, Any]]:
    query = build_query(symbol, species)
    params = build_common_params(email)
    params.update({"db": "gene", "term": query, "retmax": str(retmax)})
    url = f"{ESEARCH_URL}?{urlencode(params)}"
    return query, url, fetch_json(url)


def summarize_gene_ids(ids: list[str], email: str | None) -> tuple[str | None, dict[str, Any]]:
    if not ids:
        return None, {"result": {"uids": []}}
    params = build_common_params(email)
    params.update({"db": "gene", "id": ",".join(ids)})
    url = f"{ESUMMARY_URL}?{urlencode(params)}"
    return url, fetch_json(url)


def compact_gene_record(record: dict[str, Any]) -> dict[str, Any]:
    genomic_info = record.get("genomicinfo") or []
    primary_genomic = genomic_info[0] if genomic_info else {}
    organism = record.get("organism") or {}
    aliases = [
        alias.strip()
        for alias in str(record.get("otheraliases") or "").split(",")
        if alias.strip()
    ]
    return {
        "gene_id": record.get("uid"),
        "symbol": record.get("name"),
        "description": record.get("description"),
        "official_name": record.get("nomenclaturename"),
        "official_symbol": record.get("nomenclaturesymbol"),
        "organism": organism.get("scientificname"),
        "common_name": organism.get("commonname"),
        "taxid": organism.get("taxid"),
        "chromosome": record.get("chromosome"),
        "map_location": record.get("maplocation"),
        "aliases": aliases,
        "mim": record.get("mim") or [],
        "summary": record.get("summary"),
        "genomic_info": {
            "chraccver": primary_genomic.get("chraccver"),
            "chrstart": primary_genomic.get("chrstart"),
            "chrstop": primary_genomic.get("chrstop"),
            "exoncount": primary_genomic.get("exoncount"),
        },
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbol", required=True, help="Gene symbol, for example BRCA1")
    parser.add_argument("--species", default="homo sapiens", help="Organism name for the query")
    parser.add_argument("--retmax", type=int, default=5, help="Maximum number of Gene IDs to request")
    parser.add_argument("--email", default=None, help="Optional contact email for NCBI etiquette")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.retmax < 1 or args.retmax > 20:
        raise SystemExit("--retmax must be between 1 and 20")

    query, search_url, search_payload = search_gene(args.symbol, args.species, args.retmax, args.email)
    ids = search_payload.get("esearchresult", {}).get("idlist", [])
    summary_url, summary_payload = summarize_gene_ids(ids, args.email)
    result_block = summary_payload.get("result", {})
    payload = {
        "symbol": args.symbol,
        "species": args.species,
        "query": query,
        "search": {
            "count": search_payload.get("esearchresult", {}).get("count"),
            "ids": ids,
            "query_translation": search_payload.get("esearchresult", {}).get("querytranslation"),
        },
        "genes": [
            compact_gene_record(result_block[uid])
            for uid in result_block.get("uids", [])
            if uid in result_block
        ],
        "source_urls": {
            "search": search_url,
            "summary": summary_url,
        },
    }
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
