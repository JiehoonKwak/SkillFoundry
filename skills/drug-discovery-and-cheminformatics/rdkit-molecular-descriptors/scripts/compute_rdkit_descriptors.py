#!/usr/bin/env python3
"""Compute a compact RDKit descriptor summary for one SMILES string."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


def descriptor_payload(smiles: str, name: str | None) -> dict:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise SystemExit(f"Invalid SMILES: {smiles}")
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    return {
        "name": name,
        "input_smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": round(Descriptors.MolWt(molecule), 4),
        "exact_molecular_weight": round(rdMolDescriptors.CalcExactMolWt(molecule), 6),
        "logp": round(Descriptors.MolLogP(molecule), 4),
        "tpsa": round(rdMolDescriptors.CalcTPSA(molecule), 4),
        "hba": int(Lipinski.NumHAcceptors(molecule)),
        "hbd": int(Lipinski.NumHDonors(molecule)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(molecule)),
        "ring_count": int(rdMolDescriptors.CalcNumRings(molecule)),
        "heavy_atom_count": int(molecule.GetNumHeavyAtoms()),
        "atom_count": int(molecule.GetNumAtoms()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="Input SMILES string")
    parser.add_argument("--name", default=None, help="Optional molecule name")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    payload = descriptor_payload(args.smiles, args.name)
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

