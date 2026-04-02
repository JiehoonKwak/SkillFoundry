#!/usr/bin/env python3
"""Standardize a molecule with RDKit MolStandardize."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_toolkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires the chem-tools prefix at slurm/envs/chem-tools.") from exc
    return Chem, rdMolDescriptors, rdMolStandardize


def formal_charge(molecule) -> int:
    return int(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()))


def standardize_payload(smiles: str, name: str | None) -> dict:
    Chem, rdMolDescriptors, rdMolStandardize = load_toolkit()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise SystemExit(f"Invalid SMILES: {smiles}")

    cleaned = rdMolStandardize.Cleanup(molecule)
    parent = rdMolStandardize.FragmentParent(cleaned)
    uncharger = rdMolStandardize.Uncharger()
    uncharged = uncharger.uncharge(parent)
    tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
    canonical_tautomer = tautomer_enumerator.Canonicalize(uncharged)

    return {
        "name": name,
        "input_smiles": smiles,
        "cleaned_smiles": Chem.MolToSmiles(cleaned, canonical=True),
        "fragment_parent_smiles": Chem.MolToSmiles(parent, canonical=True),
        "uncharged_smiles": Chem.MolToSmiles(uncharged, canonical=True),
        "canonical_tautomer_smiles": Chem.MolToSmiles(canonical_tautomer, canonical=True),
        "charge_before": formal_charge(molecule),
        "charge_after": formal_charge(canonical_tautomer),
        "removed_fragments": Chem.MolToSmiles(cleaned, canonical=True) != Chem.MolToSmiles(parent, canonical=True),
        "formula": rdMolDescriptors.CalcMolFormula(canonical_tautomer),
        "heavy_atom_count": int(canonical_tautomer.GetNumHeavyAtoms()),
    }


def resolve_smiles(smiles: str | None, smiles_file: Path | None) -> str:
    if smiles is not None:
        return smiles
    if smiles_file is None:
        raise SystemExit("Provide --smiles or --smiles-file.")
    for line in smiles_file.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    raise SystemExit(f"No SMILES found in file: {smiles_file}")


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--smiles", help="Input SMILES string to standardize.")
    input_group.add_argument("--smiles-file", type=Path, help="Path to a file containing one SMILES string.")
    parser.add_argument("--name", default=None, help="Optional molecule name.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = standardize_payload(resolve_smiles(args.smiles, args.smiles_file), args.name)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
