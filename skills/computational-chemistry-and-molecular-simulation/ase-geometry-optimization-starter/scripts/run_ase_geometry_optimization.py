#!/usr/bin/env python3
"""Optimize a tiny argon dimer with ASE and summarize the result."""

from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.optimize import BFGS


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def run_optimization(input_path: Path) -> dict:
    config = load_config(input_path)
    atoms = Atoms("Ar2", positions=[[0.0, 0.0, 0.0], [config["initial_distance"], 0.0, 0.0]])
    atoms.calc = LennardJones(epsilon=config["epsilon"], sigma=config["sigma"], rc=10.0)
    initial_distance = atoms.get_distance(0, 1)
    initial_energy = atoms.get_potential_energy()
    logfile = io.StringIO()
    optimizer = BFGS(atoms, logfile=logfile)
    optimizer.run(fmax=1e-6, steps=80)
    final_distance = atoms.get_distance(0, 1)
    final_energy = atoms.get_potential_energy()
    return {
        "input_path": str(input_path.resolve()),
        "initial_distance": round(float(initial_distance), 6),
        "final_distance": round(float(final_distance), 6),
        "initial_energy": round(float(initial_energy), 8),
        "final_energy": round(float(final_energy), 8),
        "energy_drop": round(float(initial_energy - final_energy), 8),
        "step_count": len([line for line in logfile.getvalue().splitlines() if line.strip()]) - 1,
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
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    payload = run_optimization(args.input)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
