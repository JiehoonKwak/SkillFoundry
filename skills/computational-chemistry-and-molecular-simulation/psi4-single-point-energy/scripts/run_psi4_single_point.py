#!/usr/bin/env python3
"""Run a deterministic Psi4 single-point energy on a tiny water molecule."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


WATER_GEOMETRY = """\
0 1
O 0.000000 0.000000 0.000000
H 0.000000 0.757160 0.586260
H 0.000000 -0.757160 0.586260
"""

WATER_ATOMS = [
    {"element": "O", "x": 0.0, "y": 0.0, "z": 0.0},
    {"element": "H", "x": 0.0, "y": 0.75716, "z": 0.58626},
    {"element": "H", "x": 0.0, "y": -0.75716, "z": 0.58626},
]


def import_psi4():
    try:
        import psi4  # type: ignore
    except ImportError as exc:  # pragma: no cover - exercised via runtime test
        raise SystemExit(
            "Psi4 is not importable. Run this script with slurm/envs/psi4/bin/python "
            "or install the dedicated prefix first."
        ) from exc
    return psi4


def build_molecule_metadata() -> dict[str, Any]:
    return {
        "name": "water",
        "formula": "H2O",
        "charge": 0,
        "multiplicity": 1,
        "atom_count": len(WATER_ATOMS),
        "geometry_angstrom": WATER_ATOMS,
    }


def run_single_point(method: str = "hf", basis: str = "sto-3g") -> dict[str, Any]:
    psi4 = import_psi4()
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.be_quiet()
    psi4.set_num_threads(1)
    psi4.set_memory("512 MB")
    psi4.set_options(
        {
            "reference": "rhf",
            "scf_type": "pk",
            "e_convergence": 1e-10,
            "d_convergence": 1e-10,
        }
    )

    molecule = psi4.geometry(WATER_GEOMETRY)
    energy = float(psi4.energy(f"{method}/{basis}", molecule=molecule))

    return {
        "psi4_version": getattr(psi4, "__version__", "unknown"),
        "job_label": f"water-{method}-{basis}".replace("/", "-"),
        "mode": "single_point",
        "method": method,
        "basis": basis,
        "molecule": build_molecule_metadata(),
        "energy_hartree": round(energy, 12),
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", default="hf", help="Electronic-structure method, default: hf")
    parser.add_argument("--basis", default="sto-3g", help="Basis set, default: sto-3g")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    payload = run_single_point(method=args.method, basis=args.basis)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
