#!/usr/bin/env python3
"""Run the repository smoke target for every registered skill and summarize the matrix."""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

from skill_suite_utils import ROOT, load_skills, map_skill_to_smoke_targets, parse_make_targets


def choose_target(skill: dict, smoke_map: dict[str, list[str]]) -> str | None:
    targets = smoke_map.get(skill["skill_id"], [])
    if not targets:
        return None
    return sorted(targets)[0]


def run_target(target: str, timeout: int) -> dict:
    started = time.monotonic()
    completed = subprocess.run(
        ["make", target],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    duration = round(time.monotonic() - started, 3)
    return {
        "target": target,
        "returncode": completed.returncode,
        "duration_seconds": duration,
        "stdout_tail": completed.stdout.strip().splitlines()[-20:],
        "stderr_tail": completed.stderr.strip().splitlines()[-20:],
    }


def write_json(payload: dict, path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if path is None:
        print(text)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output path.")
    parser.add_argument("--markdown-out", type=Path, default=None, help="Optional Markdown output path.")
    parser.add_argument("--timeout", type=int, default=1200, help="Per-target timeout in seconds.")
    parser.add_argument("--dry-run", action="store_true", help="Only resolve the skill-to-smoke mapping.")
    args = parser.parse_args()

    skills = load_skills()
    smoke_targets = parse_make_targets()
    smoke_map = map_skill_to_smoke_targets(skills, smoke_targets)

    records: list[dict] = []
    for skill in skills:
        target = choose_target(skill, smoke_map)
        record = {
            "skill_id": skill["skill_id"],
            "slug": skill["slug"],
            "status": skill["status"],
            "target": target,
        }
        if args.dry_run or target is None:
            record["result"] = None
        else:
            record["result"] = run_target(target, timeout=args.timeout)
        records.append(record)

    passed = 0
    failed: list[str] = []
    missing: list[str] = []
    for record in records:
        result = record["result"]
        if record["target"] is None:
            missing.append(record["skill_id"])
        elif result is None:
            passed += 1
        elif result["returncode"] == 0:
            passed += 1
        else:
            failed.append(record["skill_id"])

    payload = {
        "summary": {
            "skill_count": len(records),
            "targets_resolved": len(records) - len(missing),
            "passed": passed,
            "failed": len(failed),
            "missing_targets": missing,
            "failed_skills": failed,
            "dry_run": args.dry_run,
        },
        "skills": records,
    }

    write_json(payload, args.json_out)
    if args.markdown_out is not None:
        lines = [
            "# Skill Smoke Matrix",
            "",
            f"- Skills evaluated: `{payload['summary']['skill_count']}`",
            f"- Targets resolved: `{payload['summary']['targets_resolved']}`",
            f"- Passed: `{payload['summary']['passed']}`",
            f"- Failed: `{payload['summary']['failed']}`",
            "",
        ]
        if failed:
            lines.append("## Failed Skills")
            lines.append("")
            for skill_id in failed:
                lines.append(f"- `{skill_id}`")
        args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if failed or missing:
        return 1
    print(
        "Skill smoke matrix passed: "
        f"{payload['summary']['passed']} skills"
        + (" (dry run mapping only)." if args.dry_run else ".")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
