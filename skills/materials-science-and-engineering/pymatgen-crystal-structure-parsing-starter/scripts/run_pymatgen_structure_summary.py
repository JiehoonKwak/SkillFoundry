#!/usr/bin/env python3
"""Parse one crystal structure file with pymatgen and summarize it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_stack():
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:  # pragma: no cover - exercised by runtime usage
        raise SystemExit(
            "This script requires the dedicated materials prefix at slurm/envs/materials. "
            "Run it with slurm/envs/materials/bin/python."
        ) from exc
    return Structure, SpacegroupAnalyzer


def summarize_structure(input_path: Path) -> dict:
    if not input_path.exists():
        raise ValueError(f"Input structure file does not exist: {input_path}")
    Structure, SpacegroupAnalyzer = load_stack()
    try:
        structure = Structure.from_file(str(input_path))
    except Exception as exc:  # pragma: no cover - parser errors are input-dependent
        raise ValueError(f"Failed to parse structure file {input_path}: {exc}") from exc

    composition = structure.composition
    lattice = structure.lattice
    species_counts = {str(element): int(count) for element, count in composition.element_composition.items()}

    summary = {
        "input_file": str(input_path.resolve()),
        "formula": structure.formula,
        "reduced_formula": composition.reduced_formula,
        "anonymous_formula": composition.anonymized_formula,
        "site_count": len(structure),
        "species_counts": dict(sorted(species_counts.items())),
        "lattice": {
            "a": round(float(lattice.a), 6),
            "b": round(float(lattice.b), 6),
            "c": round(float(lattice.c), 6),
            "alpha": round(float(lattice.alpha), 6),
            "beta": round(float(lattice.beta), 6),
            "gamma": round(float(lattice.gamma), 6),
        },
        "volume": round(float(structure.volume), 6),
        "density": round(float(structure.density), 6),
    }

    try:
        analyzer = SpacegroupAnalyzer(structure)
        summary["symmetry"] = {
            "space_group_symbol": analyzer.get_space_group_symbol(),
            "space_group_number": int(analyzer.get_space_group_number()),
            "crystal_system": analyzer.get_crystal_system(),
        }
    except Exception as exc:  # pragma: no cover - deterministic example succeeds
        summary["symmetry"] = {
            "space_group_symbol": None,
            "space_group_number": None,
            "crystal_system": None,
            "warning": str(exc),
        }
    return summary


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Path to a structure file such as CIF.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        payload = summarize_structure(args.input)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
