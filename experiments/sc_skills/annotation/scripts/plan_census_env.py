#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


CENSUS_REQUIRES_PYTHON = ">=3.10,<3.13"


def parse_version(text: str) -> tuple[int, int, int]:
    parts = text.split(".")
    if len(parts) < 2:
        raise ValueError(f"Invalid version string: {text}")
    normalized = [int(parts[0]), int(parts[1]), int(parts[2]) if len(parts) > 2 else 0]
    return tuple(normalized[:3])  # type: ignore[return-value]


def build_plan(*, python_version: str, conda_available: bool, python312_available: bool) -> dict[str, object]:
    parsed = parse_version(python_version)
    host_supported = (3, 10) <= parsed < (3, 13)

    if host_supported:
        recommended_route = "host_python_census_api"
        rationale = "The host interpreter is already inside the published support window for cellxgene-census."
    elif conda_available:
        recommended_route = "conda_py312_census_env"
        rationale = "The host interpreter is blocked, but conda is available so a supported Python 3.12 environment can be created locally."
    elif python312_available:
        recommended_route = "python312_venv_census_env"
        rationale = "The host interpreter is blocked, but python3.12 is present so a dedicated venv can be created locally."
    else:
        recommended_route = "discover_export_fallback"
        rationale = "No supported local Python runtime is available, so the agent should use the CELLxGENE Discover metadata-export route and local ranking helpers."

    return {
        "host_python_version": python_version,
        "cellxgene_census_requires_python": CENSUS_REQUIRES_PYTHON,
        "host_python_supported_for_cellxgene_census": host_supported,
        "conda_available": conda_available,
        "python312_available": python312_available,
        "recommended_route": recommended_route,
        "rationale": rationale,
        "routes": {
            "host_python_census_api": {
                "when": "Use this only when the current Python is already between 3.10 and 3.12.",
                "commands": [
                    "python3 -m pip install --upgrade pip",
                    "python3 -m pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
                ],
            },
            "conda_py312_census_env": {
                "when": "Preferred fallback on Python 3.13 hosts when conda is available.",
                "commands": [
                    "conda create -y -n sc-annotation-py312 python=3.12",
                    "conda activate sc-annotation-py312",
                    "pip install --upgrade pip",
                    "pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
                ],
            },
            "python312_venv_census_env": {
                "when": "Use when python3.12 exists on the host but conda is unavailable.",
                "commands": [
                    "python3.12 -m venv .venv-census312",
                    "source .venv-census312/bin/activate",
                    "pip install --upgrade pip",
                    "pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
                ],
            },
            "discover_export_fallback": {
                "when": "Use when a supported Census Python runtime cannot be created locally.",
                "commands": [
                    "Open the CELLxGENE Discover datasets browser in a real browser.",
                    "Search for '<species> <tissue> normal' and apply filters for organism, tissue, disease/control state, and suspension_type=cell.",
                    "Export metadata columns matching experiments/sc_skills/annotation/examples/reference_metadata_export.tsv.",
                    "Run: python3 experiments/sc_skills/annotation/scripts/rank_reference_metadata.py --metadata <export.tsv> --species <species> --tissue <tissue> --condition normal --out 04_reference_ranked.tsv --summary-out 05_reference_selection.json",
                    "Download the selected dataset source H5AD from the dataset page or record its source URI for later materialization in a supported environment.",
                ],
            },
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Plan a supported CELLxGENE/Census runtime route for the annotation workflow.")
    parser.add_argument("--python-version", default=sys.version.split()[0], help="Host Python version string to evaluate.")
    parser.add_argument("--conda", choices=["auto", "yes", "no"], default="auto")
    parser.add_argument("--python312", choices=["auto", "yes", "no"], default="auto")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.conda == "auto":
        conda_available = shutil.which("conda") is not None
    else:
        conda_available = args.conda == "yes"
    if args.python312 == "auto":
        python312_available = shutil.which("python3.12") is not None
    else:
        python312_available = args.python312 == "yes"

    payload = build_plan(
        python_version=args.python_version,
        conda_available=conda_available,
        python312_available=python312_available,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(str(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
