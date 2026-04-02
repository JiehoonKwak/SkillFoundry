#!/usr/bin/env python3
"""Assign a simple OpenMM force field to a tiny PDB topology."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from openmm import NonbondedForce, Platform, unit
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "two_waters.pdb"


def summarize_assignment(input_path: Path) -> dict:
    pdb = PDBFile(str(input_path))
    topology = pdb.topology
    force_field = ForceField("tip3p.xml")
    system = force_field.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    forces = [force.__class__.__name__ for force in system.getForces()]
    total_mass = sum(
        system.getParticleMass(index).value_in_unit(unit.dalton)
        for index in range(system.getNumParticles())
    )
    nonbonded_particles = 0
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded_particles = force.getNumParticles()
            break
    residues = [residue.name for residue in topology.residues()]
    return {
        "input_file": str(input_path),
        "openmm_version": Platform.getOpenMMVersion(),
        "residue_count": len(residues),
        "residue_names": residues,
        "atom_count": topology.getNumAtoms(),
        "bond_count": topology.getNumBonds(),
        "particle_count": system.getNumParticles(),
        "constraint_count": system.getNumConstraints(),
        "force_classes": forces,
        "nonbonded_particle_count": nonbonded_particles,
        "total_mass_dalton": round(float(total_mass), 6),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="PDB topology to parameterize.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input PDB not found: {args.input}")

    payload = summarize_assignment(args.input.resolve())
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
