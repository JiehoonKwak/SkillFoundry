#!/usr/bin/env python3
"""Search PubMed via NCBI E-utilities and fetch summary metadata."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_json(url: str) -> dict:
    request = Request(url, headers={"User-Agent": "SciSkillUniverse/0.2"})
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
    raise RuntimeError(f"PubMed request failed after retries: {last_error}")


def build_common_params(email: str | None) -> dict[str, str]:
    params = {"tool": "SciSkillUniverse", "retmode": "json"}
    if email:
        params["email"] = email
    return params


def search_pubmed(term: str, retmax: int, email: str | None) -> dict:
    params = build_common_params(email)
    params.update({"db": "pubmed", "term": term, "retmax": str(retmax)})
    return fetch_json(f"{ESEARCH_URL}?{urlencode(params)}")


def summarize_pubmed(ids: list[str], email: str | None) -> dict:
    if not ids:
        return {}
    params = build_common_params(email)
    params.update({"db": "pubmed", "id": ",".join(ids)})
    return fetch_json(f"{ESUMMARY_URL}?{urlencode(params)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--term", required=True, help="PubMed search term")
    parser.add_argument("--retmax", type=int, default=5, help="Maximum number of PMIDs to request")
    parser.add_argument("--email", default=None, help="Optional contact email for NCBI etiquette")
    parser.add_argument("--out", type=Path, default=None, help="Optional output path")
    args = parser.parse_args()

    if args.retmax < 1 or args.retmax > 20:
        raise SystemExit("--retmax must be between 1 and 20")

    search_payload = search_pubmed(args.term, args.retmax, args.email)
    ids = search_payload.get("esearchresult", {}).get("idlist", [])
    summary_payload = summarize_pubmed(ids, args.email)
    payload = {"term": args.term, "search": search_payload, "summary": summary_payload}
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
