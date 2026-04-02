#!/usr/bin/env python3
"""Fetch a concise Reactome event or pathway summary by stable ID."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


API_ROOT = "https://reactome.org/ContentService/data/query"
DEFAULT_STABLE_ID = "R-HSA-141409"
USER_AGENT = "SciSkillUniverse/0.2"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "r_hsa_141409_summary.json"


class ReactomeRequestError(RuntimeError):
    """Raised when the live Reactome API cannot be queried reliably."""


def normalize_stable_id(stable_id: str) -> str:
    normalized = stable_id.strip().upper()
    if not normalized:
        raise ValueError("Stable ID must not be empty.")
    return normalized


def fetch_event(stable_id: str) -> dict[str, Any]:
    normalized = normalize_stable_id(stable_id)
    request = Request(
        f"{API_ROOT}/{quote(normalized, safe='')}",
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomeRequestError(f"Reactome request failed for {normalized}: HTTP {exc.code} {message}") from exc
        except URLError as exc:
            if attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomeRequestError(f"Reactome request failed for {normalized}: {exc.reason}") from exc
    raise ReactomeRequestError(f"Reactome request failed for {normalized}: exhausted retries")


def _first_text(items: Any) -> str | None:
    if not isinstance(items, list):
        return None
    for item in items:
        if isinstance(item, dict):
            text = item.get("text")
            if isinstance(text, str) and text.strip():
                return text.strip()
    return None


def _display_name(record: Any) -> str | None:
    if not isinstance(record, dict):
        return None
    for key in ("displayName", "name"):
        value = record.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def _participants(records: Any) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    participants: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        participants.append(
            {
                "stable_id": record.get("stId"),
                "display_name": _display_name(record),
                "schema_class": record.get("schemaClass"),
                "db_id": record.get("dbId"),
            }
        )
    return participants


def _review_status(record: Any) -> str | None:
    if isinstance(record, str) and record.strip():
        return record.strip()
    if not isinstance(record, dict):
        return None
    return _display_name(record) or record.get("definition")


def _literature_references(records: Any) -> list[dict[str, Any]]:
    if not isinstance(records, list):
        return []
    references: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        references.append(
            {
                "display_name": _display_name(record),
                "pubmed_id": record.get("pubMedIdentifier"),
                "db_id": record.get("dbId"),
            }
        )
    return references


def build_summary(payload: dict[str, Any]) -> dict[str, Any]:
    stable_id = payload.get("stId") or payload.get("stIdVersion") or DEFAULT_STABLE_ID
    summary = {
        "stable_id": stable_id,
        "stable_id_version": payload.get("stIdVersion"),
        "db_id": payload.get("dbId"),
        "display_name": payload.get("displayName"),
        "schema_class": payload.get("schemaClass") or payload.get("className"),
        "species": payload.get("speciesName") or _display_name(payload.get("species")),
        "review_status": _review_status(payload.get("reviewStatus")),
        "is_inferred": bool(payload.get("isInferred")),
        "is_in_disease": bool(payload.get("isInDisease")),
        "release_date": payload.get("releaseDate"),
        "category": payload.get("category"),
        "compartments": [
            name
            for name in (_display_name(record) for record in payload.get("compartment") or [])
            if name
        ],
        "inputs": _participants(payload.get("input")),
        "outputs": _participants(payload.get("output")),
        "literature_references": _literature_references(payload.get("literatureReference")),
        "summary_text": _first_text(payload.get("summation")),
        "source_url": f"{API_ROOT}/{quote(str(stable_id), safe='')}",
    }
    return summary


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def load_asset_fallback(stable_id: str, reason: str) -> dict[str, Any] | None:
    if normalize_stable_id(stable_id) != DEFAULT_STABLE_ID or not ASSET_FALLBACK.exists():
        return None
    payload = json.loads(ASSET_FALLBACK.read_text(encoding="utf-8"))
    payload["fallback_reason"] = reason
    payload["result_origin"] = "asset_fallback"
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stable-id",
        default=DEFAULT_STABLE_ID,
        help=f"Reactome stable ID, for example {DEFAULT_STABLE_ID}",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    try:
        payload = fetch_event(args.stable_id)
        summary = build_summary(payload)
    except ReactomeRequestError as exc:
        fallback = load_asset_fallback(args.stable_id, str(exc))
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        summary = fallback
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
