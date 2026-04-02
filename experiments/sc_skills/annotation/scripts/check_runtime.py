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
    "scanpy.external",
    "cellxgene_census",
    "harmonypy",
    "celltypist",
]

MODULE_PURPOSES = {
    "anndata": "AnnData I/O and metadata inspection",
    "scanpy": "preprocessing, PCA, gene scoring, and neighbor graphs",
    "scanpy.external": "Harmony integration entry point via scanpy.external.pp.harmony_integrate",
    "cellxgene_census": "live CELLxGENE Census dataset discovery and source H5AD download",
    "harmonypy": "Harmony backend if scanpy.external.harmony_integrate is available",
    "celltypist": "optional probability and majority-voting review anchor",
}

KNOWN_INSTALL_HINTS = {
    "anndata": "pip install anndata",
    "scanpy": "pip install 'scanpy[leiden]'",
    "scanpy.external": "pip install 'scanpy[leiden]' harmonypy",
    "cellxgene_census": "pip install cellxgene-census",
    "harmonypy": "pip install harmonypy",
    "celltypist": "pip install celltypist",
}

CENSUS_REQUIRES_PYTHON = ">=3.10,<3.13"


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
    host_version = sys.version_info[:3]
    host_version_str = ".".join(str(part) for part in host_version)
    python_supported_for_census = (3, 10) <= host_version < (3, 13)
    conda_available = shutil.which("conda") is not None
    python312_available = shutil.which("python3.12") is not None

    bootstrap_commands = [
        {
            "name": "conda-census-py312",
            "when": "Preferred when the host interpreter is Python 3.13 or when cellxgene-census is missing.",
            "command": "conda create -y -n sc-annotation-py312 python=3.12 && conda activate sc-annotation-py312 && pip install --upgrade pip && pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
        },
        {
            "name": "python312-venv-census",
            "when": "Use this when python3.12 exists on the host but conda is unavailable.",
            "command": "python3.12 -m venv .venv-census312 && source .venv-census312/bin/activate && pip install --upgrade pip && pip install 'scanpy[leiden]' anndata pandas numpy scipy scikit-learn cellxgene-census harmonypy celltypist",
        },
        {
            "name": "discover-export-fallback",
            "when": "Use this when no supported Python 3.10-3.12 runtime can be created locally.",
            "command": "Use the CELLxGENE Discover browser to search/filter datasets, export metadata to TSV with the columns in examples/reference_metadata_export.tsv, then rank candidates locally with scripts/rank_reference_metadata.py.",
        },
    ]

    fallback_guidance = [
        {
            "missing_if_any": ["host_python_support_for_cellxgene_census"],
            "next_step": "The local interpreter is outside the current cellxgene-census support window. Create a Python 3.12 conda/venv environment first, or switch to the Discover metadata-export fallback route.",
        },
        {
            "missing_if_any": ["cellxgene_census"],
            "next_step": "Either install cellxgene-census in a supported Python 3.10-3.12 environment or export CELLxGENE Discover metadata and rank it locally with scripts/rank_reference_metadata.py.",
        },
        {
            "missing_if_any": ["scanpy.external", "harmonypy"],
            "next_step": "Use the deterministic Harmony-style surrogate embedded in scripts/run_annotation.py only to validate workflow shape; do not treat it as production batch correction.",
        },
    ]

    return {
        "python_executable": sys.executable,
        "python_version": host_version_str,
        "cellxgene_census_requires_python": CENSUS_REQUIRES_PYTHON,
        "host_python_supported_for_cellxgene_census": python_supported_for_census,
        "conda_available": conda_available,
        "python312_available": python312_available,
        "checked_modules": checked,
        "available_modules": available_modules,
        "missing_modules": missing_modules,
        "bootstrap_commands": bootstrap_commands,
        "fallback_guidance": fallback_guidance,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check whether the public runtime stack for the annotation workflow is present.")
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
