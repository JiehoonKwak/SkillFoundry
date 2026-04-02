#!/usr/bin/env python3
"""Audit selected resource URLs from the registry and summarize reachability."""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
REGISTRY = ROOT / "registry"
USER_AGENT = "SciSkillUniverseLinkAudit/1.0"


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def selected_resources(resource_ids: list[str] | None) -> list[dict]:
    resources = load_jsonl(REGISTRY / "resources_dedup.jsonl")
    if not resource_ids:
        return resources
    wanted = set(resource_ids)
    selected = [resource for resource in resources if resource["resource_id"] in wanted]
    missing = sorted(wanted - {resource["resource_id"] for resource in selected})
    if missing:
        raise SystemExit(f"Unknown resource IDs: {', '.join(missing)}")
    return selected


def fetch_status(url: str, timeout: float = 12.0, attempts: int = 3) -> dict:
    last_result = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status_code = getattr(response, "status", 200)
                final_url = response.geturl()
                content_type = response.headers.get("Content-Type")
                return {
                    "ok": 200 <= status_code < 400,
                    "status_code": int(status_code),
                    "final_url": final_url,
                    "content_type": content_type,
                    "error": None,
                }
        except urllib.error.HTTPError as exc:
            last_result = {
                "ok": False,
                "status_code": int(exc.code),
                "final_url": url,
                "content_type": exc.headers.get("Content-Type"),
                "error": f"HTTPError: {exc.reason}",
            }
            if exc.code not in {429, 500, 502, 503, 504} or attempt == attempts:
                return last_result
        except urllib.error.URLError as exc:
            last_result = {
                "ok": False,
                "status_code": None,
                "final_url": url,
                "content_type": None,
                "error": f"URLError: {exc.reason}",
            }
            if attempt == attempts:
                return last_result
        time.sleep(0.5 * attempt)
    assert last_result is not None
    return last_result


def audit_resources(resource_ids: list[str] | None) -> dict:
    results = []
    for resource in selected_resources(resource_ids):
        status = fetch_status(resource["url"])
        results.append(
            {
                "resource_id": resource["resource_id"],
                "canonical_name": resource["canonical_name"],
                "url": resource["url"],
                **status,
            }
        )
    ok_count = sum(1 for item in results if item["ok"])
    failing = [item["resource_id"] for item in results if not item["ok"]]
    return {
        "checked_count": len(results),
        "ok_count": ok_count,
        "failing_count": len(failing),
        "failing_resource_ids": failing,
        "results": results,
    }


def write_payload(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resource-id", action="append", default=None, help="Registry resource ID to audit. Repeatable.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    payload = audit_resources(resource_ids=args.resource_id)
    write_payload(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
