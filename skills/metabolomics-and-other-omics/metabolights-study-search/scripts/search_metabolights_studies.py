#!/usr/bin/env python3
"""Search MetaboLights studies and return compact metadata summaries."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


LIST_API_ROOT = "https://www.ebi.ac.uk/metabolights/ws/studies"
DETAIL_API_ROOT = "https://www.ebi.ac.uk/metabolights/ws/studies"
DEFAULT_QUERY = "diabetes"
USER_AGENT = "SciSkillUniverse/0.2"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "diabetes_study.json"


class MetaboLightsRequestError(RuntimeError):
    """Raised when the live MetaboLights API cannot be queried reliably."""


def normalize_query(query: str) -> str:
    normalized = " ".join(query.split())
    if not normalized:
        raise ValueError("Query must not be empty.")
    return normalized


def normalize_page(page: int) -> int:
    if page <= 0:
        raise ValueError("Page must be a positive integer.")
    return page


def normalize_rows(rows: int) -> int:
    if rows <= 0:
        raise ValueError("Rows must be a positive integer.")
    return rows


def build_query_url(query: str, page: int, rows: int) -> str:
    params = urlencode({"query": query, "page": page, "rows": rows})
    return f"{LIST_API_ROOT}?{params}"


def _request_json(url: str) -> dict[str, Any]:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                payload = json.load(response)
            break
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            last_error = exc
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise MetaboLightsRequestError(f"MetaboLights request failed for {url}: HTTP {exc.code} {message}") from exc
        except URLError as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise MetaboLightsRequestError(f"MetaboLights request failed for {url}: {exc.reason}") from exc
    else:  # pragma: no cover - guarded by explicit raises above
        raise MetaboLightsRequestError(f"MetaboLights request failed for {url}: {last_error}")
    if not isinstance(payload, dict):
        raise SystemExit(f"Unexpected MetaboLights response shape from {url}: {type(payload).__name__}")
    return payload


def fetch_accessions(query: str, page: int, rows: int) -> dict[str, Any]:
    normalized_query = normalize_query(query)
    normalized_page = normalize_page(page)
    normalized_rows = normalize_rows(rows)
    return _request_json(build_query_url(normalized_query, normalized_page, normalized_rows))


def fetch_study_detail(accession: str) -> dict[str, Any]:
    return _request_json(f"{DETAIL_API_ROOT}/{accession}")


def _clean_text(value: Any) -> str | None:
    if isinstance(value, str):
        normalized = " ".join(value.split())
        return normalized or None
    return None


def _factor_names(study_record: dict[str, Any]) -> list[str]:
    factors = study_record.get("factors")
    if not isinstance(factors, list):
        return []
    values: list[str] = []
    for factor in factors:
        if not isinstance(factor, dict):
            continue
        name = _clean_text(factor.get("factorName"))
        if name and name not in values:
            values.append(name)
    return values[:8]


def summarize_study(accession: str, payload: dict[str, Any]) -> dict[str, Any]:
    mtbls_study = payload.get("mtblsStudy")
    isa = payload.get("isaInvestigation")
    if not isinstance(mtbls_study, dict) or not isinstance(isa, dict):
        raise SystemExit(f"Unexpected MetaboLights detail shape for {accession}")
    studies = isa.get("studies")
    study_record = studies[0] if isinstance(studies, list) and studies else {}
    if not isinstance(study_record, dict):
        study_record = {}
    return {
        "accession": accession,
        "title": _clean_text(study_record.get("title")) or _clean_text(isa.get("title")),
        "description": _clean_text(study_record.get("description")) or _clean_text(isa.get("description")),
        "study_status": _clean_text(mtbls_study.get("studyStatus")),
        "study_category": _clean_text(mtbls_study.get("studyCategory")),
        "submission_date": _clean_text(study_record.get("submissionDate")) or _clean_text(isa.get("submissionDate")),
        "public_release_date": _clean_text(study_record.get("publicReleaseDate")) or _clean_text(isa.get("publicReleaseDate")),
        "factor_names": _factor_names(study_record),
        "assay_count": len(study_record.get("assays", [])) if isinstance(study_record.get("assays"), list) else 0,
        "publication_count": len(study_record.get("publications", [])) if isinstance(study_record.get("publications"), list) else 0,
        "person_count": len(study_record.get("people", [])) if isinstance(study_record.get("people"), list) else 0,
        "study_url": _clean_text(mtbls_study.get("studyHttpUrl")),
        "ftp_url": _clean_text(mtbls_study.get("studyFtpUrl")),
    }


def build_summary(search_payload: dict[str, Any], query: str, page: int, rows: int) -> dict[str, Any]:
    accessions = search_payload.get("content")
    total_studies = search_payload.get("studies")
    if not isinstance(accessions, list):
        raise SystemExit("Unexpected MetaboLights response shape: missing content accession list.")
    limited_accessions = accessions[: normalize_rows(rows)]
    summaries = [summarize_study(accession, fetch_study_detail(accession)) for accession in limited_accessions]
    return {
        "query": normalize_query(query),
        "page": normalize_page(page),
        "rows": normalize_rows(rows),
        "total_studies": total_studies,
        "study_count": len(summaries),
        "studies": summaries,
        "source_url": LIST_API_ROOT,
        "query_url": build_query_url(query, page, rows),
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def load_asset_fallback(query: str, page: int, rows: int, reason: str) -> dict[str, Any] | None:
    if normalize_query(query).lower() != DEFAULT_QUERY:
        return None
    if normalize_page(page) != 1 or normalize_rows(rows) != 1 or not ASSET_FALLBACK.exists():
        return None
    payload = json.loads(ASSET_FALLBACK.read_text(encoding="utf-8"))
    payload["fallback_reason"] = reason
    payload["result_origin"] = "asset_fallback"
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", default=DEFAULT_QUERY, help="MetaboLights study query text.")
    parser.add_argument("--page", type=int, default=1, help="One-based results page.")
    parser.add_argument("--rows", type=int, default=3, help="Maximum number of accessions to resolve.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        search_payload = fetch_accessions(args.query, args.page, args.rows)
        summary = build_summary(search_payload, args.query, args.page, args.rows)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    except MetaboLightsRequestError as exc:
        fallback = load_asset_fallback(args.query, args.page, args.rows, str(exc))
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        summary = fallback
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
