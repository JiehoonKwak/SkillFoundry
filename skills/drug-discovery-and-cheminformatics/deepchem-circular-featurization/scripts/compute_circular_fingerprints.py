#!/usr/bin/env python3
"""Compute compact DeepChem circular fingerprints for a small list of SMILES strings."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
from pathlib import Path
from typing import Any


DEFAULT_SMILES = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
]


def normalize_smiles_list(smiles_values: list[str] | None) -> list[str]:
    values = DEFAULT_SMILES if not smiles_values else smiles_values
    normalized: list[str] = []
    for item in values:
        cleaned = item.strip()
        if not cleaned:
            raise ValueError("SMILES values must not be empty.")
        normalized.append(cleaned)
    return normalized


def normalize_positive_int(value: int, label: str) -> int:
    if value <= 0:
        raise ValueError(f"{label} must be a positive integer.")
    return value


def load_toolkit():
    capture_out = io.StringIO()
    capture_err = io.StringIO()
    try:
        with contextlib.redirect_stdout(capture_out), contextlib.redirect_stderr(capture_err):
            import deepchem as dc
            from rdkit import Chem, RDLogger, rdBase
    except ImportError as exc:
        raise SystemExit(
            "DeepChem circular featurization requires the dedicated prefix at "
            "slurm/envs/deepchem with both deepchem and rdkit installed."
        ) from exc
    RDLogger.DisableLog("rdApp.*")
    return dc, Chem, rdBase


def featurize_smiles(smiles_values: list[str], size: int, radius: int) -> dict[str, Any]:
    dc, Chem, rdBase = load_toolkit()
    featurizer = dc.feat.CircularFingerprint(size=size, radius=radius)
    molecules = []
    for index, smiles in enumerate(smiles_values, start=1):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise SystemExit(f"Invalid SMILES at position {index}: {smiles}")
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        fingerprint = featurizer.featurize([molecule])[0]
        bit_vector = [int(value) for value in fingerprint.tolist()]
        on_bits = [bit_index for bit_index, value in enumerate(bit_vector) if value]
        molecules.append(
            {
                "molecule_index": index,
                "input_smiles": smiles,
                "canonical_smiles": canonical_smiles,
                "bit_vector": bit_vector,
                "on_bits": on_bits,
                "on_bit_count": len(on_bits),
            }
        )
    return {
        "featurizer": "CircularFingerprint",
        "deepchem_version": dc.__version__,
        "rdkit_version": rdBase.rdkitVersion,
        "size": size,
        "radius": radius,
        "molecule_count": len(molecules),
        "molecules": molecules,
    }


def write_json(payload: dict[str, Any], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smiles",
        action="append",
        default=None,
        help="Input SMILES string. Repeat for multiple molecules. Defaults to aspirin and caffeine if omitted.",
    )
    parser.add_argument("--size", type=int, default=32, help="Fingerprint size in bits.")
    parser.add_argument("--radius", type=int, default=2, help="Circular fingerprint radius.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    try:
        smiles_values = normalize_smiles_list(args.smiles)
        size = normalize_positive_int(args.size, "size")
        radius = normalize_positive_int(args.radius, "radius")
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    payload = featurize_smiles(smiles_values, size=size, radius=radius)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
