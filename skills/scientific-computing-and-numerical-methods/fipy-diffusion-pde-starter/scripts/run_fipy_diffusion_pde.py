#!/usr/bin/env python3
"""Solve a small 1D diffusion problem with FiPy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from fipy import CellVariable, DiffusionTerm, Grid1D, TransientTerm


def build_summary() -> dict[str, object]:
    mesh = Grid1D(nx=20, dx=1.0)
    phi = CellVariable(name="phi", mesh=mesh, value=0.0)
    phi.setValue(1.0, where=mesh.x < 5.0)
    equation = TransientTerm() == DiffusionTerm(coeff=0.5)
    initial_mass = float(phi.value.sum())
    for _ in range(10):
        equation.solve(var=phi, dt=0.2)
    return {
        "cell_count": int(len(phi.value)),
        "initial_mass": round(initial_mass, 6),
        "final_mass": round(float(phi.value.sum()), 6),
        "center_value": round(float(phi.value[len(phi.value) // 2]), 6),
        "min_value": round(float(phi.value.min()), 6),
        "max_value": round(float(phi.value.max()), 6),
        "leading_profile": [round(float(value), 6) for value in phi.value[:5]],
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = build_summary()
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
