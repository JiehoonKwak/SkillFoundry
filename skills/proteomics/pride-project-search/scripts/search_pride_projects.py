#!/usr/bin/env python3
"""Search PRIDE Archive v2 projects and return compact summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects"
DEFAULT_KEYWORD = "phosphoproteomics"
USER_AGENT = "SciSkillUniverse/0.2"


def normalize_keyword(keyword: str) -> str:
    normalized = keyword.strip()
    if not normalized:
        raise ValueError("Keyword must not be empty.")
    return normalized


def normalize_page_size(page_size: int) -> int:
    if page_size <= 0:
        raise ValueError("Page size must be a positive integer.")
    return page_size


def normalize_page(page: int) -> int:
    if page < 0:
        raise ValueError("Page must be zero or greater.")
    return page


def build_query_url(keyword: str, page_size: int, page: int) -> str:
    query = urlencode({"keyword": keyword, "page": page, "pageSize": page_size})
    return f"{API_ROOT}?{query}"


def fetch_projects(keyword: str, page_size: int, page: int) -> list[dict[str, Any]]:
    normalized_keyword = normalize_keyword(keyword)
    normalized_page_size = normalize_page_size(page_size)
    normalized_page = normalize_page(page)
    request = Request(
        build_query_url(normalized_keyword, normalized_page_size, normalized_page),
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except HTTPError as exc:
        message = exc.read().decode("utf-8", errors="replace").strip()
        raise SystemExit(f"PRIDE request failed for keyword {keyword!r}: HTTP {exc.code} {message}") from exc
    except URLError as exc:
        raise SystemExit(f"PRIDE request failed for keyword {keyword!r}: {exc.reason}") from exc

    if not isinstance(payload, list):
        raise SystemExit(f"Unexpected PRIDE response shape: expected list, got {type(payload).__name__}")
    return payload


def _clean_text(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = " ".join(value.split())
    if not normalized:
        return None
    return normalized


def _string_list(items: Any, *, limit: int | None = None) -> list[str]:
    if not isinstance(items, list):
        return []
    values: list[str] = []
    for item in items:
        cleaned = _clean_text(item)
        if cleaned and cleaned not in values:
            values.append(cleaned)
        if limit is not None and len(values) >= limit:
            break
    return values


def summarize_project(project: dict[str, Any]) -> dict[str, Any]:
    accession = _clean_text(project.get("accession"))
    return {
        "accession": accession,
        "title": _clean_text(project.get("title")),
        "project_description": _clean_text(project.get("projectDescription")),
        "submission_type": _clean_text(project.get("submissionType")),
        "publication_date": _clean_text(project.get("publicationDate")),
        "submission_date": _clean_text(project.get("submissionDate")),
        "updated_date": _clean_text(project.get("updatedDate")),
        "organisms": _string_list(project.get("organisms"), limit=10),
        "keywords": _string_list(project.get("keywords"), limit=10),
        "references": _string_list(project.get("references"), limit=10),
        "experiment_types": _string_list(project.get("experimentTypes"), limit=10),
        "project_url": f"https://www.ebi.ac.uk/pride/archive/projects/{accession}" if accession else None,
    }


def build_summary(projects: list[dict[str, Any]], keyword: str, page_size: int, page: int) -> dict[str, Any]:
    normalized_keyword = normalize_keyword(keyword)
    normalized_page_size = normalize_page_size(page_size)
    normalized_page = normalize_page(page)
    summaries = [summarize_project(project) for project in projects]

    return {
        "keyword": normalized_keyword,
        "page": normalized_page,
        "page_size": normalized_page_size,
        "project_count": len(summaries),
        "projects": summaries,
        "source_url": API_ROOT,
        "query_url": build_query_url(normalized_keyword, normalized_page_size, normalized_page),
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
    parser.add_argument("--keyword", default=DEFAULT_KEYWORD, help="PRIDE project keyword query.")
    parser.add_argument("--page-size", type=int, default=3, help="Maximum number of projects to return.")
    parser.add_argument("--page", type=int, default=0, help="Zero-based PRIDE results page.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        projects = fetch_projects(args.keyword, args.page_size, args.page)
        summary = build_summary(projects, args.keyword, args.page_size, args.page)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
