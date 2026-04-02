#!/usr/bin/env python3
"""Compare two toy phylogenetic trees with scikit-bio."""

from __future__ import annotations

import argparse
import json
from io import StringIO
from pathlib import Path

from skbio import TreeNode


LEFT_NEWICK = "((A:1,B:1):1,(C:1,D:1):1);"
RIGHT_NEWICK = "((A:1,C:1):1,(B:1,D:1):1);"


def total_branch_length(tree: TreeNode) -> float:
    return float(sum(node.length or 0.0 for node in tree.postorder()))


def build_summary() -> dict[str, object]:
    left = TreeNode.read(StringIO(LEFT_NEWICK))
    right = TreeNode.read(StringIO(RIGHT_NEWICK))
    left_tips = [tip.name for tip in left.tips()]
    right_tips = [tip.name for tip in right.tips()]
    return {
        "left_tip_names": left_tips,
        "right_tip_names": right_tips,
        "shared_tip_count": len(sorted(set(left_tips) & set(right_tips))),
        "left_total_branch_length": round(total_branch_length(left), 6),
        "right_total_branch_length": round(total_branch_length(right), 6),
        "robinson_foulds_distance": round(float(left.compare_rfd(right)), 6),
        "weighted_robinson_foulds_distance": round(float(left.compare_wrfd(right)), 6),
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    payload = build_summary()
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
