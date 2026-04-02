#!/usr/bin/env python3
"""Run a deterministic OpenMM minimization on a tiny toy system."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from openmm import Context, HarmonicBondForce, LocalEnergyMinimizer, NonbondedForce, System, VerletIntegrator, unit
from openmm import Platform


def scalar_energy(value) -> float:
    return float(value.value_in_unit(unit.kilojoule_per_mole))


def vector_nm(position) -> list[float]:
    return [round(float(component.value_in_unit(unit.nanometer)), 6) for component in position]


def build_system() -> tuple[System, list]:
    system = System()
    for _ in range(3):
        system.addParticle(12.0)

    bonds = HarmonicBondForce()
    bonds.addBond(0, 1, 0.15, 300.0)
    bonds.addBond(1, 2, 0.15, 300.0)
    system.addForce(bonds)

    nonbonded = NonbondedForce()
    for charge in (0.0, 0.0, 0.0):
        nonbonded.addParticle(charge, 0.30, 0.20)
    system.addForce(nonbonded)

    positions = [
        unit.Quantity((0.00, 0.00, 0.00), unit.nanometer),
        unit.Quantity((0.45, 0.00, 0.00), unit.nanometer),
        unit.Quantity((0.95, 0.10, 0.00), unit.nanometer),
    ]
    return system, positions


def summarize_minimization() -> dict:
    system, positions = build_system()
    integrator = VerletIntegrator(0.001)
    platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, platform)
    context.setPositions(positions)

    initial_state = context.getState(getEnergy=True, getPositions=True)
    LocalEnergyMinimizer.minimize(context, tolerance=1.0)
    final_state = context.getState(getEnergy=True, getPositions=True)

    initial_positions = [vector_nm(position) for position in initial_state.getPositions()]
    final_positions = [vector_nm(position) for position in final_state.getPositions()]

    return {
        "openmm_version": Platform.getOpenMMVersion(),
        "platform": platform.getName(),
        "particle_count": 3,
        "initial_potential_energy_kj_mol": round(scalar_energy(initial_state.getPotentialEnergy()), 6),
        "final_potential_energy_kj_mol": round(scalar_energy(final_state.getPotentialEnergy()), 6),
        "initial_positions_nm": initial_positions,
        "final_positions_nm": final_positions,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    payload = summarize_minimization()
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

