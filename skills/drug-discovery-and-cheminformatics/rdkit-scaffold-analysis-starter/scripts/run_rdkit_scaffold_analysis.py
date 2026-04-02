#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold


SCRIPT_PATH = Path(__file__).resolve()
SKILL_DIR = SCRIPT_PATH.parents[1]
DEFAULT_INPUT = SKILL_DIR / "examples" / "molecules.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze Murcko scaffolds for a small TSV of molecules."
    )
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
        help="TSV with columns name and smiles",
    )
    parser.add_argument(
        "--out",
        help="Optional JSON output path",
    )
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
            if not name:
                raise ValueError(f"Missing name at row {row_number}.")
            if not smiles:
                raise ValueError(f"Missing smiles at row {row_number}.")
            rows.append({"row_number": row_number, "name": name, "smiles": smiles})
    if not rows:
        raise ValueError("Input TSV does not contain any molecule rows.")
    return rows


def scaffold_smiles(mol: Chem.Mol) -> str | None:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(scaffold, canonical=True)


def generic_scaffold_smiles(murcko_smiles: str | None) -> str | None:
    if not murcko_smiles:
        return None
    scaffold = Chem.MolFromSmiles(murcko_smiles)
    if scaffold is None:
        return None
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(generic, canonical=True)


def analyze_rows(rows: list[dict[str, str]]) -> dict[str, object]:
    molecules: list[dict[str, object]] = []
    murcko_groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    generic_groups: dict[str, list[dict[str, object]]] = defaultdict(list)

    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            raise ValueError(
                f"Invalid SMILES at row {row['row_number']} ({row['name']}): {row['smiles']}"
            )
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        murcko = scaffold_smiles(mol)
        generic = generic_scaffold_smiles(murcko)
        molecule_summary = {
            "row_number": row["row_number"],
            "name": row["name"],
            "input_smiles": row["smiles"],
            "canonical_smiles": canonical_smiles,
            "murcko_scaffold": murcko,
            "generic_scaffold": generic,
        }
        molecules.append(molecule_summary)
        murcko_groups[murcko or "(none)"].append(molecule_summary)
        generic_groups[generic or "(none)"].append(molecule_summary)

    molecules.sort(key=lambda item: (item["row_number"], item["name"]))

    scaffold_groups = []
    for murcko_key in sorted(
        murcko_groups,
        key=lambda key: (-len(murcko_groups[key]), key),
    ):
        members = sorted(entry["name"] for entry in murcko_groups[murcko_key])
        generic_values = {
            entry["generic_scaffold"]
            for entry in murcko_groups[murcko_key]
            if entry["generic_scaffold"] is not None
        }
        scaffold_groups.append(
            {
                "murcko_scaffold": None if murcko_key == "(none)" else murcko_key,
                "generic_scaffold": sorted(generic_values)[0] if generic_values else None,
                "count": len(members),
                "members": members,
            }
        )

    generic_scaffold_groups = []
    for generic_key in sorted(
        generic_groups,
        key=lambda key: (-len(generic_groups[key]), key),
    ):
        members = sorted(entry["name"] for entry in generic_groups[generic_key])
        murcko_values = {
            entry["murcko_scaffold"]
            for entry in generic_groups[generic_key]
            if entry["murcko_scaffold"] is not None
        }
        generic_scaffold_groups.append(
            {
                "generic_scaffold": None if generic_key == "(none)" else generic_key,
                "count": len(members),
                "members": members,
                "murcko_scaffolds": sorted(murcko_values),
            }
        )

    summary = {
        "molecule_count": len(molecules),
        "valid_molecule_count": len(molecules),
        "unique_murcko_scaffolds": sum(
            1 for key in murcko_groups if key != "(none)"
        ),
        "unique_generic_scaffolds": sum(
            1 for key in generic_groups if key != "(none)"
        ),
        "largest_scaffold_group_size": max(
            (group["count"] for group in scaffold_groups),
            default=0,
        ),
    }

    return {
        "rdkit_version": rdBase.rdkitVersion,
        "molecules": molecules,
        "scaffold_groups": scaffold_groups,
        "generic_scaffold_groups": generic_scaffold_groups,
        "summary": summary,
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
        rows = read_rows(input_path)
        payload = analyze_rows(rows)
    except ValueError as exc:
        return error(str(exc))

    payload = {
        "input_file": str(input_path),
        **payload,
    }
    write_payload(payload, out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
