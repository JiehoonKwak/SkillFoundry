#!/usr/bin/env python3
"""Search ChEMBL molecules and return a compact compound summary."""

from __future__ import annotations

import argparse
import json
import socket
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API_ROOT = "https://www.ebi.ac.uk/chembl/api/data/molecule/search"
DEFAULT_QUERY = "imatinib"
USER_AGENT = "SciSkillUniverse/0.2"
ASSETS_DIR = Path(__file__).resolve().parents[1] / "assets"


def normalize_query(query: str) -> str:
    normalized = query.strip()
    if not normalized:
        raise ValueError("Query must not be empty.")
    return normalized


def build_request_url(query: str, limit: int) -> str:
    params = urlencode({"q": query, "limit": limit, "format": "json"})
    return f"{API_ROOT}?{params}"


def fetch_search(query: str, limit: int, retries: int = 4, timeout: float = 10.0) -> dict[str, Any]:
    normalized = normalize_query(query)
    request = Request(
        build_request_url(query=normalized, limit=limit),
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except HTTPError as exc:
            last_error = exc
            if exc.code in {408, 429, 500, 502, 503, 504} and attempt + 1 < retries:
                retry_after = exc.headers.get("Retry-After")
                delay = float(retry_after) if retry_after else float(1 + attempt)
                time.sleep(delay)
                continue
            message = exc.read().decode("utf-8", errors="replace").strip()
            raise SystemExit(f"ChEMBL request failed for {normalized}: HTTP {exc.code} {message}") from exc
        except URLError as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(float(1 + attempt))
                continue
            raise SystemExit(f"ChEMBL request failed for {normalized}: {exc.reason}") from exc
        except (TimeoutError, socket.timeout) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(float(1 + attempt))
                continue
            raise SystemExit(f"ChEMBL request failed for {normalized}: {exc}") from exc
    raise RuntimeError(f"ChEMBL request failed after retries: {last_error}")


def _dedupe_synonyms(records: Any, limit: int = 5) -> list[str]:
    values: list[str] = []
    for record in records or []:
        if not isinstance(record, dict):
            continue
        value = record.get("molecule_synonym")
        if isinstance(value, str) and value and value not in values:
            values.append(value)
        if len(values) >= limit:
            break
    return values


def _molecule_summary(record: dict[str, Any]) -> dict[str, Any]:
    properties = record.get("molecule_properties") or {}
    structures = record.get("molecule_structures") or {}
    return {
        "chembl_id": record.get("molecule_chembl_id"),
        "preferred_name": record.get("pref_name"),
        "molecule_type": record.get("molecule_type"),
        "max_phase": record.get("max_phase"),
        "availability_type": record.get("availability_type"),
        "alogp": properties.get("alogp"),
        "full_mwt": properties.get("full_mwt"),
        "molecular_species": properties.get("molecular_species"),
        "canonical_smiles": structures.get("canonical_smiles"),
        "standard_inchi_key": structures.get("standard_inchi_key"),
        "synonyms": _dedupe_synonyms(record.get("molecule_synonyms")),
    }


def build_summary(payload: dict[str, Any], query: str, limit: int) -> dict[str, Any]:
    molecules = payload.get("molecules") or []
    page_meta = payload.get("page_meta") or {}
    summaries = [_molecule_summary(record) for record in molecules if isinstance(record, dict)]
    return {
        "query": query,
        "limit": limit,
        "page_meta": page_meta,
        "result_count": len(summaries),
        "molecules": summaries,
        "source_url": build_request_url(query=query, limit=limit),
    }


def load_cached_summary(query: str, limit: int) -> dict[str, Any] | None:
    normalized = normalize_query(query)
    for asset_path in sorted(ASSETS_DIR.glob("*.json")):
        payload = json.loads(asset_path.read_text(encoding="utf-8"))
        if payload.get("query") == normalized and payload.get("limit") == limit:
            return payload
    return None


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", default=DEFAULT_QUERY, help="Free-text molecule query")
    parser.add_argument("--limit", type=int, default=3, help="Maximum number of molecules to return")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    if args.limit < 1:
        raise SystemExit("--limit must be positive.")

    normalized_query = normalize_query(args.query)
    try:
        payload = fetch_search(normalized_query, args.limit)
        summary = build_summary(payload, query=normalized_query, limit=args.limit)
    except SystemExit as exc:
        summary = load_cached_summary(normalized_query, args.limit)
        if summary is None:
            raise
        summary = dict(summary)
        summary["fallback_reason"] = str(exc)
        summary["source_mode"] = "cached_fallback"
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
