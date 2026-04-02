#!/usr/bin/env python3
"""Search public ENCODE experiments and return compact metadata summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://www.encodeproject.org/search/"
DEFAULT_SEARCH_TERM = "ATAC-seq"
USER_AGENT = "SciSkillUniverse/0.2"


def normalize_search_term(search_term: str) -> str:
    normalized = " ".join(search_term.split())
    if not normalized:
        raise ValueError("Search term must not be empty.")
    return normalized


def normalize_limit(limit: int) -> int:
    if limit <= 0:
        raise ValueError("Limit must be a positive integer.")
    return limit


def build_query_url(search_term: str, limit: int) -> str:
    query = urlencode(
        {
            "searchTerm": search_term,
            "type": "Experiment",
            "format": "json",
            "limit": limit,
            "frame": "object",
        }
    )
    return f"{API_ROOT}?{query}"


def fetch_experiments(search_term: str, limit: int) -> dict[str, Any]:
    normalized_search_term = normalize_search_term(search_term)
    normalized_limit = normalize_limit(limit)
    request = Request(
        build_query_url(normalized_search_term, normalized_limit),
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
        raise SystemExit(f"ENCODE request failed for term {search_term!r}: HTTP {exc.code} {message}") from exc
    except URLError as exc:
        raise SystemExit(f"ENCODE request failed for term {search_term!r}: {exc.reason}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"Unexpected ENCODE response shape: expected dict, got {type(payload).__name__}")
    return payload


def _clean_text(value: Any) -> str | None:
    if isinstance(value, str):
        normalized = " ".join(value.split())
        return normalized or None
    return None


def _nested_text(item: dict[str, Any], *keys: str) -> str | None:
    current: Any = item
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return _clean_text(current)


def summarize_experiment(experiment: dict[str, Any]) -> dict[str, Any]:
    accession = _clean_text(experiment.get("accession"))
    experiment_path = _clean_text(experiment.get("@id"))
    files = experiment.get("files")
    file_count = len(files) if isinstance(files, list) else None
    return {
        "accession": accession,
        "assay_title": _clean_text(experiment.get("assay_title")),
        "biosample_summary": _clean_text(experiment.get("biosample_summary"))
        or _nested_text(experiment, "biosample_ontology", "term_name"),
        "target_label": _nested_text(experiment, "target", "label")
        or _nested_text(experiment, "target", "title"),
        "assembly": _clean_text(experiment.get("assembly")),
        "status": _clean_text(experiment.get("status")),
        "lab_title": _nested_text(experiment, "lab", "title"),
        "award_project": _nested_text(experiment, "award", "project"),
        "file_count": file_count,
        "experiment_url": f"https://www.encodeproject.org{experiment_path}" if experiment_path else None,
    }


def build_summary(payload: dict[str, Any], search_term: str, limit: int) -> dict[str, Any]:
    experiments = payload.get("@graph")
    if not isinstance(experiments, list):
        raise SystemExit("Unexpected ENCODE response shape: missing @graph list.")
    summaries = [summarize_experiment(item) for item in experiments]
    return {
        "search_term": normalize_search_term(search_term),
        "limit": normalize_limit(limit),
        "total": payload.get("total"),
        "experiment_count": len(summaries),
        "experiments": summaries,
        "source_url": API_ROOT,
        "query_url": build_query_url(search_term, limit),
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
    parser.add_argument("--search-term", default=DEFAULT_SEARCH_TERM, help="ENCODE free-text search term.")
    parser.add_argument("--limit", type=int, default=3, help="Maximum number of experiments to summarize.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        payload = fetch_experiments(args.search_term, args.limit)
        summary = build_summary(payload, args.search_term, args.limit)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
