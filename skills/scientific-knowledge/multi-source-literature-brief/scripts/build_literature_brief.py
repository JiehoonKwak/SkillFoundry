#!/usr/bin/env python3
"""Build a normalized literature brief across multiple official scholarly APIs."""

from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


OPENALEX_URL = "https://api.openalex.org/works"
EUROPEPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF_URL = "https://api.crossref.org/works"
NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
USER_AGENT = "SciSkillUniverse/0.3"


def fetch_json(url: str, retries: int = 3) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt + 1 == retries:
                raise
            retry_after = exc.headers.get("Retry-After")
            time.sleep(float(retry_after) if retry_after else float(1 + attempt))
        except URLError as exc:
            last_error = exc
            if attempt + 1 == retries:
                raise
            time.sleep(float(1 + attempt))
    raise RuntimeError(f"Request failed after retries: {last_error}")


def first_list_value(value: Any) -> Any:
    if isinstance(value, list) and value:
        return value[0]
    return None


def normalize_doi(value: Any) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip()
    if text.lower().startswith("https://doi.org/"):
        text = text[16:]
    return text.lower()


def normalize_title(value: Any) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return None
    return " ".join(value.lower().split())


def short_authors(names: list[str], limit: int = 5) -> list[str]:
    return [name for name in names if name][:limit]


def query_openalex(query: str, limit: int, mailto: str | None) -> tuple[str, dict[str, Any], list[dict[str, Any]]]:
    params = {"search": query, "per-page": str(limit)}
    if mailto:
        params["mailto"] = mailto
    url = f"{OPENALEX_URL}?{urlencode(params)}"
    payload = fetch_json(url)
    results = []
    for record in payload.get("results", [])[:limit]:
        results.append(
            {
                "source": "openalex",
                "source_id": record.get("id"),
                "title": record.get("display_name"),
                "year": record.get("publication_year"),
                "doi": normalize_doi((record.get("ids") or {}).get("doi")),
                "authors": short_authors(
                    [
                        ((item.get("author") or {}).get("display_name"))
                        for item in record.get("authorships") or []
                    ]
                ),
                "venue": (((record.get("primary_location") or {}).get("source") or {}).get("display_name")),
                "identifier": record.get("id"),
                "url": record.get("id"),
            }
        )
    return url, payload, results


def query_europepmc(query: str, limit: int) -> tuple[str, dict[str, Any], list[dict[str, Any]]]:
    params = {"query": query, "pageSize": str(limit), "format": "json"}
    url = f"{EUROPEPMC_URL}?{urlencode(params)}"
    payload = fetch_json(url)
    results = []
    for record in ((payload.get("resultList") or {}).get("result") or [])[:limit]:
        results.append(
            {
                "source": "europepmc",
                "source_id": record.get("id"),
                "title": record.get("title"),
                "year": record.get("pubYear"),
                "doi": normalize_doi(record.get("doi")),
                "authors": short_authors([record.get("authorString")]),
                "venue": record.get("journalTitle"),
                "identifier": f"{record.get('source')}:{record.get('id')}",
                "url": f"https://europepmc.org/article/{record.get('source')}/{record.get('id')}",
            }
        )
    return url, payload, results


def query_crossref(query: str, limit: int) -> tuple[str, dict[str, Any], list[dict[str, Any]]]:
    url = f"{CROSSREF_URL}?{urlencode({'query.title': query, 'rows': limit})}"
    payload = fetch_json(url)
    results = []
    for record in ((payload.get("message") or {}).get("items") or [])[:limit]:
        published = first_list_value((record.get("published-print") or {}).get("date-parts"))
        if not published:
            published = first_list_value((record.get("published-online") or {}).get("date-parts"))
        if not published:
            published = first_list_value((record.get("created") or {}).get("date-parts"))
        results.append(
            {
                "source": "crossref",
                "source_id": record.get("DOI"),
                "title": first_list_value(record.get("title")),
                "year": first_list_value(published) if isinstance(published, list) else None,
                "doi": normalize_doi(record.get("DOI")),
                "authors": short_authors(
                    [
                        " ".join(part for part in [author.get("given"), author.get("family")] if part)
                        for author in record.get("author") or []
                    ]
                ),
                "venue": first_list_value(record.get("container-title")),
                "identifier": record.get("DOI"),
                "url": record.get("URL"),
            }
        )
    return url, payload, results


def build_ncbi_params(email: str | None) -> dict[str, str]:
    params = {"tool": "SciSkillUniverse", "retmode": "json"}
    if email:
        params["email"] = email
    return params


def query_pubmed(query: str, limit: int, email: str | None) -> tuple[str, str | None, dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    search_params = build_ncbi_params(email)
    search_params.update({"db": "pubmed", "term": query, "retmax": str(limit)})
    search_url = f"{NCBI_ESEARCH_URL}?{urlencode(search_params)}"
    search_payload = fetch_json(search_url, retries=4)
    ids = (search_payload.get("esearchresult") or {}).get("idlist") or []

    summary_payload: dict[str, Any] = {}
    summary_url = None
    if ids:
        summary_params = build_ncbi_params(email)
        summary_params.update({"db": "pubmed", "id": ",".join(ids)})
        summary_url = f"{NCBI_ESUMMARY_URL}?{urlencode(summary_params)}"
        summary_payload = fetch_json(summary_url, retries=4)

    results = []
    summary_result = summary_payload.get("result") or {}
    for uid in summary_result.get("uids", [])[:limit]:
        record = summary_result.get(uid) or {}
        results.append(
            {
                "source": "pubmed",
                "source_id": uid,
                "title": record.get("title"),
                "year": record.get("pubdate"),
                "doi": None,
                "authors": short_authors([author.get("name") for author in record.get("authors") or []]),
                "venue": record.get("source"),
                "identifier": uid,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/",
            }
        )
    return search_url, summary_url, search_payload, summary_payload, results


def build_overlap_summary(items: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    doi_groups: dict[str, list[str]] = defaultdict(list)
    title_groups: dict[str, list[str]] = defaultdict(list)
    for item in items:
        label = f"{item['source']}:{item['identifier']}"
        doi = normalize_doi(item.get("doi"))
        title = normalize_title(item.get("title"))
        if doi:
            doi_groups[doi].append(label)
        if title:
            title_groups[title].append(label)
    return {
        "doi_matches": [
            {"doi": doi, "matches": matches}
            for doi, matches in sorted(doi_groups.items())
            if len(matches) > 1
        ],
        "title_matches": [
            {"title": title, "matches": matches}
            for title, matches in sorted(title_groups.items())
            if len(matches) > 1
        ],
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
    parser.add_argument("--query", required=True, help="Shared literature query")
    parser.add_argument("--limit", type=int, default=3, help="Maximum normalized results to keep per source")
    parser.add_argument("--email", default=None, help="Optional contact email for NCBI etiquette")
    parser.add_argument("--mailto", default=None, help="Optional OpenAlex polite-pool email")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.limit < 1 or args.limit > 10:
        raise SystemExit("--limit must be between 1 and 10")

    openalex_url, openalex_payload, openalex_results = query_openalex(args.query, args.limit, args.mailto)
    europepmc_url, europepmc_payload, europepmc_results = query_europepmc(args.query, args.limit)
    crossref_url, crossref_payload, crossref_results = query_crossref(args.query, args.limit)
    pubmed_search_url, pubmed_summary_url, pubmed_search, pubmed_summary, pubmed_results = query_pubmed(args.query, args.limit, args.email)

    merged_results = openalex_results + europepmc_results + crossref_results + pubmed_results
    payload = {
        "query": args.query,
        "limit_per_source": args.limit,
        "counts": {
            "openalex": {"raw_hits": (openalex_payload.get("meta") or {}).get("count"), "normalized_results": len(openalex_results)},
            "europepmc": {"raw_hits": ((europepmc_payload.get("hitCount")) or (europepmc_payload.get("request") or {}).get("pageSize")), "normalized_results": len(europepmc_results)},
            "crossref": {"raw_hits": ((crossref_payload.get("message") or {}).get("total-results")), "normalized_results": len(crossref_results)},
            "pubmed": {"raw_hits": ((pubmed_search.get("esearchresult") or {}).get("count")), "normalized_results": len(pubmed_results)},
        },
        "sources": {
            "openalex": openalex_results,
            "europepmc": europepmc_results,
            "crossref": crossref_results,
            "pubmed": pubmed_results,
        },
        "overlap_summary": build_overlap_summary(merged_results),
        "merged_results": merged_results,
        "source_urls": {
            "openalex": openalex_url,
            "europepmc": europepmc_url,
            "crossref": crossref_url,
            "pubmed_search": pubmed_search_url,
            "pubmed_summary": pubmed_summary_url,
        },
    }
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
