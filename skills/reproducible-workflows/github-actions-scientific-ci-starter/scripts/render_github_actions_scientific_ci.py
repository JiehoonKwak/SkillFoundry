#!/usr/bin/env python3
"""Render a conservative GitHub Actions workflow for this repository."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_SMOKE_TARGETS = ["smoke-zarr", "smoke-openmm-md", "smoke-optuna"]


def build_workflow(smoke_targets: list[str]) -> str:
    smoke_block = "\n".join(f"          make {target}" for target in smoke_targets)
    return (
        "name: SciSkillUniverse CI\n\n"
        "on:\n"
        "  push:\n"
        "    branches:\n"
        "      - main\n"
        "  pull_request:\n\n"
        "jobs:\n"
        "  validate-test:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v4\n"
        "      - uses: actions/setup-python@v5\n"
        "        with:\n"
        "          python-version: \"3.11\"\n"
        "      - name: Repository validation\n"
        "        run: make validate\n"
        "      - name: Build site\n"
        "        run: make build-site\n"
        "      - name: Full test suite\n"
        "        run: make test\n"
        "      - name: Targeted smoke commands\n"
        "        run: |\n"
        f"{smoke_block}\n"
    )


def summarize_workflow(workflow_text: str, workflow_out: Path, smoke_targets: list[str]) -> dict:
    return {
        "workflow_path": str(workflow_out),
        "job_names": ["validate-test"],
        "smoke_targets": smoke_targets,
        "smoke_target_count": len(smoke_targets),
        "contains_validate": "make validate" in workflow_text,
        "contains_build_site": "make build-site" in workflow_text,
        "contains_make_test": "make test" in workflow_text,
        "contains_checkout": "actions/checkout@v4" in workflow_text,
        "contains_setup_python": "actions/setup-python@v5" in workflow_text,
        "line_count": len(workflow_text.splitlines()),
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-out", type=Path, required=True, help="Rendered workflow YAML output path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument(
        "--smoke-target",
        action="append",
        default=None,
        help="Optional make smoke target to include. Repeat for multiple targets.",
    )
    args = parser.parse_args()

    smoke_targets = args.smoke_target or DEFAULT_SMOKE_TARGETS
    workflow_text = build_workflow(smoke_targets)
    args.workflow_out.parent.mkdir(parents=True, exist_ok=True)
    args.workflow_out.write_text(workflow_text, encoding="utf-8")
    write_json(summarize_workflow(workflow_text, args.workflow_out, smoke_targets), args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
