#!/usr/bin/env python
"""Generate a small deterministic conformer ensemble with RDKit."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "molecules.tsv"
DEFAULT_RANDOM_SEED = 20260314


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="TSV with columns name and smiles.")
    parser.add_argument("--num-confs", type=int, default=4, help="How many conformers to embed per molecule.")
    parser.add_argument("--out", help="Optional JSON output path.")
    return parser.parse_args()


def error(message: str) -> int:
    print(message, file=sys.stderr)
    return 1


def read_rows(input_path: Path) -> list[dict[str, str]]:
    if not input_path.exists():
        raise ValueError(f"Input file does not exist: {input_path}")
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input file is empty or missing a header row.")
        required = {"name", "smiles"}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Input TSV must contain columns {sorted(required)}; missing {sorted(missing)}."
            )
        rows: list[dict[str, str]] = []
        for row_number, row in enumerate(reader, start=2):
            name = (row.get("name") or "").strip()
            smiles = (row.get("smiles") or "").strip()
            if not name or not smiles:
                raise ValueError(f"Missing name or smiles at row {row_number}.")
            rows.append({"row_number": row_number, "name": name, "smiles": smiles})
    if not rows:
        raise ValueError("Input TSV does not contain any molecule rows.")
    return rows


def summarize_molecule(row: dict[str, str], num_confs: int) -> dict[str, object]:
    if num_confs < 1:
        raise ValueError("--num-confs must be positive.")
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        raise ValueError(f"Invalid SMILES at row {row['row_number']} ({row['name']}): {row['smiles']}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = DEFAULT_RANDOM_SEED
    params.pruneRmsThresh = -1.0
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params))
    if not conf_ids:
        raise ValueError(f"Failed to embed any conformers for {row['name']}.")
    optimization = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=500)
    ranked = []
    for conf_id, (_, energy) in zip(conf_ids, optimization, strict=True):
        ranked.append({"conformer_id": int(conf_id), "uff_energy": round(float(energy), 6)})
    ranked.sort(key=lambda item: (item["uff_energy"], item["conformer_id"]))
    lowest = ranked[0]
    return {
        "name": row["name"],
        "input_smiles": row["smiles"],
        "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True),
        "heavy_atom_count": int(Chem.RemoveHs(mol).GetNumHeavyAtoms()),
        "conformer_count": len(conf_ids),
        "lowest_energy_conformer_id": lowest["conformer_id"],
        "lowest_uff_energy": lowest["uff_energy"],
        "ranked_conformers": ranked,
    }


def build_summary(input_path: Path, num_confs: int) -> dict[str, object]:
    rows = read_rows(input_path)
    molecules = [summarize_molecule(row, num_confs=num_confs) for row in rows]
    return {
        "input_file": str(input_path),
        "num_confs_requested": num_confs,
        "rdkit_version": rdBase.rdkitVersion,
        "molecule_count": len(molecules),
        "molecules": molecules,
    }


def write_payload(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=False)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).resolve()
    out_path = Path(args.out).resolve() if args.out else None
    try:
        payload = build_summary(input_path, num_confs=args.num_confs)
    except ValueError as exc:
        return error(str(exc))
    write_payload(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
