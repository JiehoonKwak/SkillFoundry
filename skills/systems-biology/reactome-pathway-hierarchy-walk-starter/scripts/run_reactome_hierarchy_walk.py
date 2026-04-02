#!/usr/bin/env python3
"""Walk the Reactome hierarchy tree for a target stable ID."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


API_ROOT = "https://reactome.org/ContentService/data/eventsHierarchy"
DEFAULT_SPECIES = "9606"
DEFAULT_STABLE_ID = "R-HSA-141409"
USER_AGENT = "SciSkillUniverse/0.2"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
ASSET_FALLBACK = Path(__file__).resolve().parents[1] / "assets" / "r_hsa_141409_hierarchy.json"


class ReactomeHierarchyError(RuntimeError):
    """Raised when the live Reactome hierarchy cannot be queried reliably."""


def normalize_species(species: str) -> str:
    normalized = species.strip()
    if not normalized:
        raise ValueError("Species identifier must not be empty.")
    return normalized


def normalize_stable_id(stable_id: str) -> str:
    normalized = stable_id.strip().upper()
    if not normalized:
        raise ValueError("Stable ID must not be empty.")
    return normalized


def fetch_hierarchy(species: str) -> list[dict[str, Any]]:
    normalized_species = normalize_species(species)
    request = Request(
        f"{API_ROOT}/{normalized_species}",
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                payload = json.load(response)
            break
        except HTTPError as exc:
            message = exc.read().decode("utf-8", errors="replace").strip()
            if exc.code in RETRYABLE_STATUS_CODES and attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomeHierarchyError(
                f"Reactome hierarchy request failed for species {normalized_species}: HTTP {exc.code} {message}"
            ) from exc
        except URLError as exc:
            if attempt < 2:
                time.sleep(attempt + 1)
                continue
            raise ReactomeHierarchyError(
                f"Reactome hierarchy request failed for species {normalized_species}: {exc.reason}"
            ) from exc
    else:  # pragma: no cover
        raise ReactomeHierarchyError(f"Reactome hierarchy request failed for species {normalized_species}")
    if not isinstance(payload, list):
        raise ReactomeHierarchyError(
            f"Unexpected Reactome hierarchy payload for species {normalized_species}: {type(payload).__name__}"
        )
    return payload


def find_path(nodes: list[dict[str, Any]], target_stable_id: str) -> list[dict[str, Any]] | None:
    normalized_target = normalize_stable_id(target_stable_id)

    def _visit(node: dict[str, Any], path: list[dict[str, Any]]) -> list[dict[str, Any]] | None:
        next_path = path + [node]
        if str(node.get("stId", "")).upper() == normalized_target:
            return next_path
        children = node.get("children") or []
        if not isinstance(children, list):
            return None
        for child in children:
            if not isinstance(child, dict):
                continue
            found = _visit(child, next_path)
            if found is not None:
                return found
        return None

    for node in nodes:
        found = _visit(node, [])
        if found is not None:
            return found
    return None


def count_descendants(node: dict[str, Any]) -> int:
    total = 0
    children = node.get("children") or []
    if not isinstance(children, list):
        return 0
    for child in children:
        if not isinstance(child, dict):
            continue
        total += 1 + count_descendants(child)
    return total


def summarize_path(species: str, path: list[dict[str, Any]]) -> dict[str, Any]:
    node = path[-1]
    children = node.get("children") or []
    direct_children = []
    if isinstance(children, list):
        for child in children[:10]:
            if not isinstance(child, dict):
                continue
            direct_children.append(
                {
                    "stable_id": child.get("stId"),
                    "display_name": child.get("name"),
                    "type": child.get("type"),
                }
            )
    return {
        "species_id": normalize_species(species),
        "species_name": path[0].get("species"),
        "stable_id": node.get("stId"),
        "display_name": node.get("name"),
        "type": node.get("type"),
        "top_level_pathway": path[0].get("name"),
        "ancestor_count": len(path) - 1,
        "ancestor_path": [
            {
                "stable_id": item.get("stId"),
                "display_name": item.get("name"),
                "type": item.get("type"),
            }
            for item in path
        ],
        "direct_children_count": len(children) if isinstance(children, list) else 0,
        "direct_children": direct_children,
        "descendant_count": count_descendants(node),
        "source_url": f"{API_ROOT}/{normalize_species(species)}",
        "result_origin": "live_api",
    }


def load_asset_fallback(species: str, stable_id: str, reason: str) -> dict[str, Any] | None:
    if normalize_species(species) != DEFAULT_SPECIES or normalize_stable_id(stable_id) != DEFAULT_STABLE_ID:
        return None
    if not ASSET_FALLBACK.exists():
        return None
    payload = json.loads(ASSET_FALLBACK.read_text(encoding="utf-8"))
    payload["fallback_reason"] = reason
    payload["result_origin"] = "asset_fallback"
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
    parser.add_argument("--species", default=DEFAULT_SPECIES, help="Reactome species identifier, default 9606.")
    parser.add_argument("--stable-id", default=DEFAULT_STABLE_ID, help="Reactome stable ID to locate in the hierarchy.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        hierarchy = fetch_hierarchy(args.species)
        path = find_path(hierarchy, args.stable_id)
        if path is None:
            raise SystemExit(f"Stable ID {normalize_stable_id(args.stable_id)} was not found in Reactome hierarchy {normalize_species(args.species)}.")
        summary = summarize_path(args.species, path)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    except ReactomeHierarchyError as exc:
        fallback = load_asset_fallback(args.species, args.stable_id, str(exc))
        if fallback is None:
            raise SystemExit(str(exc)) from exc
        summary = fallback
    write_json(summary, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
