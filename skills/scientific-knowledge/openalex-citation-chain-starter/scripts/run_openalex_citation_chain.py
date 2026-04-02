#!/usr/bin/env python3
"""Fetch a small citation chain summary for a seed OpenAlex work."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qsl, quote, urlencode, urlparse, urlunparse
from urllib.request import Request, urlopen


API_ROOT = "https://api.openalex.org/works"
DEFAULT_WORK_ID = "10.1038/nature12373"
RETRYABLE_STATUS_CODES = {408, 425, 429, 500, 502, 503, 504}
USER_AGENT = "SciSkillUniverse/0.2"
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "hallmarks_citation_chain.json"


class OpenAlexCitationError(RuntimeError):
    """Raised when the OpenAlex citation chain query fails."""


def normalize_work_id(raw: str) -> str:
    normalized = raw.strip()
    if not normalized:
        raise ValueError("Work identifier must not be empty.")
    return normalized


def _is_openalex_work_id(value: str) -> bool:
    if value.startswith("https://openalex.org/W"):
        return True
    if value.startswith("W") and value[1:].isdigit():
        return True
    return False


def build_lookup_url(work_id: str, mailto: str | None = None) -> str:
    normalized = normalize_work_id(work_id)
    if _is_openalex_work_id(normalized):
        target = normalized if normalized.startswith("https://") else f"https://openalex.org/{normalized}"
    elif normalized.startswith("https://doi.org/") or normalized.startswith("http://doi.org/"):
        target = normalized.replace("http://", "https://", 1)
    elif normalized.startswith("10."):
        target = f"https://doi.org/{normalized}"
    else:
        target = normalized
    url = f"{API_ROOT}/{quote(target, safe='')}"
    if mailto:
        url = f"{url}?{urlencode({'mailto': mailto})}"
    return url


def append_params(url: str, **params: Any) -> str:
    parsed = urlparse(url)
    query = dict(parse_qsl(parsed.query, keep_blank_values=True))
    for key, value in params.items():
        if value is None:
            continue
        query[key] = str(value)
    return urlunparse(parsed._replace(query=urlencode(query)))


def fetch_json(url: str) -> dict[str, Any]:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise OpenAlexCitationError(f"OpenAlex request failed: HTTP {exc.code} {message}") from exc
        except URLError as exc:
            if attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise OpenAlexCitationError(f"OpenAlex request failed: {exc.reason}") from exc
    raise OpenAlexCitationError(f"OpenAlex request failed after retries: {url}")


def work_summary(payload: dict[str, Any]) -> dict[str, Any]:
    authorships = payload.get("authorships") or []
    authors = []
    for authorship in authorships[:5]:
        author = authorship.get("author") or {}
        display_name = author.get("display_name")
        if display_name:
            authors.append(display_name)
    concepts = []
    for concept in payload.get("concepts") or []:
        display_name = concept.get("display_name")
        score = concept.get("score")
        if display_name is None:
            continue
        concepts.append({"display_name": display_name, "score": score})
    return {
        "id": payload.get("id"),
        "title": payload.get("display_name"),
        "doi": payload.get("doi"),
        "publication_year": payload.get("publication_year"),
        "type": payload.get("type"),
        "cited_by_count": payload.get("cited_by_count"),
        "referenced_works_count": len(payload.get("referenced_works") or []),
        "authors": authors,
        "concepts": concepts[:5],
    }


def fetch_work_by_id(work_id: str, mailto: str | None = None) -> dict[str, Any]:
    return fetch_json(build_lookup_url(work_id, mailto=mailto))


def fetch_work_list(filter_key: str, openalex_id: str, limit: int, mailto: str | None = None) -> list[dict[str, Any]]:
    payload = fetch_json(
        append_params(
            API_ROOT,
            filter=f"{filter_key}:{openalex_id.rsplit('/', 1)[-1]}",
            **{"per-page": limit, "mailto": mailto},
        )
    )
    results = payload.get("results") or []
    return [work_summary(result) for result in results if isinstance(result, dict)][:limit]


def build_summary(seed_payload: dict[str, Any], limit: int, mailto: str | None = None) -> dict[str, Any]:
    openalex_id = str(seed_payload.get("id") or "")
    if not openalex_id:
        raise OpenAlexCitationError("Seed work did not expose an OpenAlex ID.")
    return {
        "seed_work": work_summary(seed_payload),
        "seed_work_resolved_id": openalex_id,
        "requested_limit": limit,
        "top_referenced_works": fetch_work_list("cited_by", openalex_id, limit=limit, mailto=mailto),
        "top_citing_works": fetch_work_list("cites", openalex_id, limit=limit, mailto=mailto),
        "source_url": build_lookup_url(openalex_id, mailto=mailto),
        "cited_by_filter_url": append_params(
            API_ROOT,
            filter=f"cited_by:{openalex_id.rsplit('/', 1)[-1]}",
            **{"per-page": limit},
        ),
        "cites_filter_url": append_params(
            API_ROOT,
            filter=f"cites:{openalex_id.rsplit('/', 1)[-1]}",
            **{"per-page": limit},
        ),
        "result_origin": "live_api",
    }


def load_asset_fallback(work_id: str, reason: str, *, limit: int = 3) -> dict[str, Any] | None:
    normalized = normalize_work_id(work_id).removeprefix("https://doi.org/")
    if normalized != DEFAULT_WORK_ID or not ASSET_FALLBACK.exists():
        return None
    payload = json.loads(ASSET_FALLBACK.read_text(encoding="utf-8"))
    payload["result_origin"] = "asset_fallback"
    payload["fallback_reason"] = reason
    payload["requested_limit"] = limit
    payload["top_referenced_works"] = (payload.get("top_referenced_works") or [])[:limit]
    payload["top_citing_works"] = (payload.get("top_citing_works") or [])[:limit]
    payload["cited_by_filter_url"] = append_params(
        API_ROOT,
        filter=f"cited_by:{payload['seed_work_resolved_id'].rsplit('/', 1)[-1]}",
        **{"per-page": limit},
    )
    payload["cites_filter_url"] = append_params(
        API_ROOT,
        filter=f"cites:{payload['seed_work_resolved_id'].rsplit('/', 1)[-1]}",
        **{"per-page": limit},
    )
    return payload


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-id", default=DEFAULT_WORK_ID, help="DOI or OpenAlex work ID to inspect.")
    parser.add_argument("--limit", type=int, default=3, help="Number of citing and referenced works to summarize.")
    parser.add_argument("--mailto", default=None, help="Optional polite-pool email.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.limit < 1 or args.limit > 10:
        raise SystemExit("--limit must be between 1 and 10.")

    try:
        seed_payload = fetch_work_by_id(args.work_id, mailto=args.mailto)
        summary = build_summary(seed_payload, limit=args.limit, mailto=args.mailto)
    except (ValueError, OpenAlexCitationError) as exc:
        fallback = load_asset_fallback(args.work_id, str(exc), limit=args.limit)
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        summary = fallback
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
