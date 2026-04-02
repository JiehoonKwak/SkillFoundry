#!/usr/bin/env python3
"""Resolve a gene symbol via Ensembl REST and fetch its lookup record."""

from __future__ import annotations

import argparse
import copy
import json
import socket
import time
from http.client import RemoteDisconnected
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


API_ROOT = "https://rest.ensembl.org"
USER_AGENT = "SciSkillUniverse/0.2"
ASSETS_DIR = Path(__file__).resolve().parents[1] / "assets"


def fetch_json(url: str, retries: int = 3, delay: float = 1.0, timeout: float = 8.0) -> dict | list:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
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
                time.sleep(delay * (attempt + 1))
                continue
            message = exc.read().decode("utf-8", errors="replace").strip()
            raise RuntimeError(f"Ensembl request failed: HTTP {exc.code} {message}") from exc
        except (URLError, TimeoutError, socket.timeout, ConnectionResetError, RemoteDisconnected) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(delay * (attempt + 1))
                continue
            reason = getattr(exc, "reason", str(exc))
            raise RuntimeError(f"Ensembl request failed: {reason}") from exc
    raise RuntimeError(f"Ensembl request failed after retries: {last_error}")


def xrefs_for_symbol(species: str, symbol: str) -> tuple[list[dict], str]:
    url = f"{API_ROOT}/xrefs/symbol/{quote(species)}/{quote(symbol)}?content-type=application/json"
    payload = fetch_json(url)
    return [item for item in payload if item.get("type") == "gene"], url


def lookup_gene_by_symbol(species: str, symbol: str) -> tuple[dict, str]:
    url = f"{API_ROOT}/lookup/symbol/{quote(species)}/{quote(symbol)}?content-type=application/json"
    return fetch_json(url), url


def fallback_xrefs_from_lookup(lookup: dict) -> list[dict]:
    gene_id = lookup.get("id")
    if not gene_id:
        return []
    return [{"id": gene_id, "type": "gene"}]


def load_cached_payload(species: str, symbol: str) -> tuple[dict | None, Path | None]:
    asset_path = ASSETS_DIR / f"{symbol.lower()}_lookup.json"
    if not asset_path.exists():
        return None, None
    payload = json.loads(asset_path.read_text(encoding="utf-8"))
    query = payload.get("query", {})
    if query.get("species") != species or str(query.get("symbol", "")).upper() != symbol.upper():
        return None, None
    return payload, asset_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbol", required=True, help="Gene symbol, for example BRCA1")
    parser.add_argument("--species", default="homo_sapiens", help="Ensembl species name")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    try:
        lookup, lookup_url = lookup_gene_by_symbol(args.species, args.symbol)
        xrefs, xrefs_url = xrefs_for_symbol(args.species, args.symbol)
        xrefs_source = "live"
        if not xrefs:
            xrefs = fallback_xrefs_from_lookup(lookup)
            xrefs_source = "lookup_fallback"
        payload = {
            "query": {"species": args.species, "symbol": args.symbol},
            "xrefs": xrefs,
            "lookup": lookup,
            "source_mode": "live",
            "source_urls": {"lookup": lookup_url, "xrefs": xrefs_url},
            "xrefs_source": xrefs_source,
        }
    except RuntimeError as exc:
        cached_payload, asset_path = load_cached_payload(args.species, args.symbol)
        if cached_payload is None or asset_path is None:
            raise SystemExit(str(exc)) from exc
        payload = copy.deepcopy(cached_payload)
        payload["source_mode"] = "asset_fallback"
        payload["fallback_reason"] = str(exc)
        payload["source_urls"] = {"asset": str(asset_path)}
        payload["xrefs_source"] = "asset"

    if not payload.get("xrefs"):
        raise SystemExit(f"No gene xrefs found for {args.symbol} in {args.species}")
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
