#!/usr/bin/env python3
"""Run a short deterministic OpenMM Langevin dynamics trajectory on a toy system."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from openmm import Context, HarmonicBondForce, LangevinIntegrator, Platform, System, unit


def build_system() -> tuple[System, list]:
    system = System()
    for _ in range(2):
        system.addParticle(39.9)
    bond_force = HarmonicBondForce()
    bond_force.addBond(0, 1, 0.30, 200.0)
    system.addForce(bond_force)
    positions = [
        unit.Quantity((0.0, 0.0, 0.0), unit.nanometer),
        unit.Quantity((0.6, 0.0, 0.0), unit.nanometer),
    ]
    return system, positions


def vector_nm(position) -> list[float]:
    return [float(component.value_in_unit(unit.nanometer)) for component in position]


def bond_length_nm(positions) -> float:
    first = vector_nm(positions[0])
    second = vector_nm(positions[1])
    return round(math.sqrt(sum((left - right) ** 2 for left, right in zip(first, second))), 6)


def potential_energy_kj_mol(state) -> float:
    return round(float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)), 6)


def run_dynamics(total_steps: int, report_interval: int, temperature_k: float, seed: int) -> dict:
    if total_steps < 1:
        raise ValueError("total_steps must be positive")
    if report_interval < 1 or report_interval > total_steps:
        raise ValueError("report_interval must be positive and not exceed total_steps")

    system, positions = build_system()
    integrator = LangevinIntegrator(
        temperature_k * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, platform)
    context.setPositions(positions)
    context.setVelocitiesToTemperature(temperature_k * unit.kelvin, seed)

    frames = []
    for step in range(0, total_steps + 1, report_interval):
        if step > 0:
            integrator.step(report_interval)
        state = context.getState(getEnergy=True, getPositions=True)
        frame_positions = state.getPositions()
        frames.append(
            {
                "step": step,
                "bond_length_nm": bond_length_nm(frame_positions),
                "potential_energy_kj_mol": potential_energy_kj_mol(state),
            }
        )

    return {
        "openmm_version": Platform.getOpenMMVersion(),
        "platform": platform.getName(),
        "temperature_k": float(temperature_k),
        "seed": int(seed),
        "total_steps": int(total_steps),
        "report_interval": int(report_interval),
        "frame_count": len(frames),
        "initial_bond_length_nm": frames[0]["bond_length_nm"],
        "final_bond_length_nm": frames[-1]["bond_length_nm"],
        "min_potential_energy_kj_mol": min(frame["potential_energy_kj_mol"] for frame in frames),
        "max_potential_energy_kj_mol": max(frame["potential_energy_kj_mol"] for frame in frames),
        "frames": frames,
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--steps", type=int, default=50, help="Total integration steps.")
    parser.add_argument("--report-interval", type=int, default=10, help="Frame reporting interval.")
    parser.add_argument("--temperature-k", type=float, default=300.0, help="Target temperature in Kelvin.")
    parser.add_argument("--seed", type=int, default=7, help="Velocity seed.")
    args = parser.parse_args()

    try:
        payload = run_dynamics(args.steps, args.report_interval, args.temperature_k, args.seed)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
