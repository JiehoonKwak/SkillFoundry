#!/usr/bin/env python3
"""Featurize SMILES strings into DeepChem molecular graph metadata."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import deepchem as dc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="TSV with molecule_id and smiles columns")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()

    with args.input.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    featurizer = dc.feat.MolGraphConvFeaturizer()
    features = featurizer.featurize([row["smiles"] for row in rows])
    payload = {
        "deepchem_version": dc.__version__,
        "input_file": str(args.input),
        "molecule_count": len(rows),
        "graphs": [
            {
                "molecule_id": row["molecule_id"],
                "smiles": row["smiles"],
                "node_count": int(graph.num_nodes),
                "edge_count": int(graph.num_edges),
                "node_feature_count": int(graph.num_node_features),
            }
            for row, graph in zip(rows, features, strict=True)
        ],
    }
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
