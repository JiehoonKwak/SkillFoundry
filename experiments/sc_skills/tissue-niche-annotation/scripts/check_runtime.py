#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import importlib.util
import json
import shutil
import sys
from pathlib import Path


DEFAULT_MODULES = [
    "anndata",
    "scanpy",
    "squidpy",
    "cellcharter",
    "utag",
]

MODULE_PURPOSES = {
    "anndata": "AnnData I/O and metadata inspection",
    "scanpy": "shared AnnData preprocessing and basic graph utilities",
    "squidpy": "spatial graph construction and neighborhood enrichment",
    "cellcharter": "neighborhood aggregation and niche clustering",
    "utag": "graph-aware microenvironment grouping",
}

KNOWN_INSTALL_HINTS = {
    "anndata": "pip install anndata",
    "scanpy": "pip install 'scanpy[leiden]'",
    "squidpy": "pip install squidpy",
    "cellcharter": "pip install cellcharter",
    "utag": "pip install 'scanpy>=1.9.3,<1.12' 'squidpy>=1.4' && pip install git+https://github.com/ElementoLab/utag.git@main",
}


def _base_distribution_name(module_name: str) -> str:
    return module_name.split(".", 1)[0]


def _module_available(module_name: str) -> bool:
    try:
        return importlib.util.find_spec(module_name) is not None
    except ModuleNotFoundError:
        return False


def _module_version(module_name: str) -> str | None:
    base_name = _base_distribution_name(module_name)
    try:
        return importlib.metadata.version(base_name)
    except importlib.metadata.PackageNotFoundError:
        try:
            module = importlib.import_module(base_name)
        except Exception:
            return None
        return getattr(module, "__version__", None)


def _module_report(module_name: str) -> dict[str, str | bool | None]:
    available = _module_available(module_name)
    return {
        "module": module_name,
        "available": available,
        "version": _module_version(module_name) if available else None,
        "purpose": MODULE_PURPOSES.get(module_name, "user-requested runtime dependency"),
        "install_hint": KNOWN_INSTALL_HINTS.get(module_name, f"pip install {module_name}"),
    }


def build_runtime_report(modules: list[str]) -> dict[str, object]:
    checked = [_module_report(module_name) for module_name in modules]
    missing_modules = [str(item["module"]) for item in checked if not item["available"]]
    available_modules = [str(item["module"]) for item in checked if item["available"]]
    conda_available = shutil.which("conda") is not None

    bootstrap_commands = [
        {
            "name": "base-spatial-stack",
            "command": "python3 -m venv .venv && source .venv/bin/activate && pip install --upgrade pip && pip install anndata scanpy pandas numpy scipy scikit-learn squidpy",
        },
        {
            "name": "conda-spatial-stack",
            "command": "conda create -y -n sc-tissue-niche python=3.12 && conda activate sc-tissue-niche && pip install --upgrade pip && pip install anndata scanpy pandas numpy scipy scikit-learn squidpy",
        },
        {
            "name": "cellcharter-stack",
            "command": "source .venv/bin/activate && pip install cellcharter",
        },
        {
            "name": "utag-stack",
            "command": "python3 -m venv .venv-utag && source .venv-utag/bin/activate && pip install --upgrade pip && pip install 'scanpy>=1.9.3,<1.12' 'squidpy>=1.4' && pip install git+https://github.com/ElementoLab/utag.git@main",
        },
    ]

    fallback_guidance = [
        {
            "missing_if_any": ["squidpy"],
            "next_step": "Use the bundled graph and neighborhood surrogate in scripts/run_tissue_niche_annotation.py only for provisional workflow rehearsal; do not claim Squidpy enrichment results.",
        },
        {
            "missing_if_any": ["cellcharter"],
            "next_step": "Continue with neighborhood profiles plus prototype-style niche scoring and label the report as no-CellCharter fallback.",
        },
        {
            "missing_if_any": ["utag"],
            "next_step": "Continue with graph-aware fallback grouping and record that the run did not execute UTAG proper.",
        },
    ]

    return {
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "conda_available": conda_available,
        "checked_modules": checked,
        "available_modules": available_modules,
        "missing_modules": missing_modules,
        "bootstrap_commands": bootstrap_commands,
        "fallback_guidance": fallback_guidance,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check whether the public runtime stack for tissue niche annotation is present.")
    parser.add_argument("--modules", nargs="*", default=DEFAULT_MODULES, help="Modules to inspect.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args(argv)

    report = build_runtime_report(args.modules)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.out is None:
        print(text)
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text + "\n", encoding="utf-8")
    print(str(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
