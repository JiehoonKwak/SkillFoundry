#!/usr/bin/env python3
"""Resolve a GBIF species match and fetch a small occurrence sample."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


USER_AGENT = "SciSkillUniverse/0.2"
SPECIES_MATCH_URL = "https://api.gbif.org/v1/species/match"
OCCURRENCE_SEARCH_URL = "https://api.gbif.org/v1/occurrence/search"
CANONICAL_QUERY = ("Puma concolor", "US", 3)
ASSET_PATH = Path(__file__).resolve().parents[1] / "assets" / "puma_concolor_us_occurrences.json"


def fetch_json(url: str, retries: int = 3) -> dict:
    for attempt in range(retries):
        request = Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except (HTTPError, URLError) as exc:
            if attempt == retries - 1:
                raise exc
            time.sleep(1.0 + attempt)
    raise AssertionError("unreachable")  # pragma: no cover


def build_summary(scientific_name: str, country: str | None, limit: int) -> dict[str, object]:
    species_params = urlencode({"name": scientific_name})
    species_match = fetch_json(f"{SPECIES_MATCH_URL}?{species_params}")
    occurrence_params = {"scientificName": scientific_name, "limit": limit}
    if country:
        occurrence_params["country"] = country
    occurrence_payload = fetch_json(f"{OCCURRENCE_SEARCH_URL}?{urlencode(occurrence_params)}")
    results = occurrence_payload.get("results", [])
    countries_seen = sorted({record.get("countryCode") for record in results if record.get("countryCode")})
    record_summaries = [
        {
            "gbif_id": record.get("gbifID"),
            "scientific_name": record.get("scientificName"),
            "country": record.get("countryCode"),
            "basis_of_record": record.get("basisOfRecord"),
        }
        for record in results
    ]
    return {
        "query_scientific_name": scientific_name,
        "country": country,
        "limit": int(limit),
        "matched_usage_key": species_match.get("usageKey"),
        "matched_scientific_name": species_match.get("scientificName"),
        "match_type": species_match.get("matchType"),
        "occurrence_count": len(results),
        "end_of_records": bool(occurrence_payload.get("endOfRecords", False)),
        "countries_seen": countries_seen,
        "record_summaries": record_summaries,
    }


def load_fallback() -> dict[str, object]:
    return json.loads(ASSET_PATH.read_text(encoding="utf-8"))


def maybe_run_with_fallback(scientific_name: str, country: str | None, limit: int) -> dict[str, object]:
    try:
        return build_summary(scientific_name=scientific_name, country=country, limit=limit)
    except (HTTPError, URLError):
        if (scientific_name, country, limit) == CANONICAL_QUERY and ASSET_PATH.exists():
            payload = load_fallback()
            payload["fallback_used"] = True
            return payload
        raise


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scientific-name", required=True, help="Scientific name to resolve and query.")
    parser.add_argument("--country", default=None, help="Optional ISO country code filter.")
    parser.add_argument("--limit", type=int, default=3, help="Number of occurrence records to request.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.limit < 1 or args.limit > 20:
        raise SystemExit("--limit must be between 1 and 20.")

    payload = maybe_run_with_fallback(
        scientific_name=args.scientific_name,
        country=args.country,
        limit=args.limit,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
