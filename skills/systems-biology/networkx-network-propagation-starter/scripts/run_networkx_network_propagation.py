#!/usr/bin/env python3
"""Run deterministic personalized PageRank on a small weighted interaction graph."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_network.tsv"
DEFAULT_SEEDS = Path(__file__).resolve().parents[1] / "examples" / "toy_seeds.txt"


def load_networkx():
    try:
        import networkx as nx
    except ImportError as exc:  # pragma: no cover - runtime surface
        raise SystemExit("This script requires networkx in the active Python environment.") from exc
    return nx


def load_edges(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"No edges found in {path}")
    return rows


def load_seeds(path: Path) -> list[str]:
    seeds = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not seeds:
        raise SystemExit(f"No seed nodes found in {path}")
    return seeds


def build_graph(edge_rows: list[dict[str, str]]):
    nx = load_networkx()
    graph = nx.DiGraph()
    for row in edge_rows:
        graph.add_edge(
            row["source"],
            row["target"],
            weight=float(row.get("weight") or 1.0),
        )
    return graph


def summarize_propagation(input_path: Path, seed_path: Path, alpha: float, top_k: int) -> dict:
    nx = load_networkx()
    graph = build_graph(load_edges(input_path))
    seeds = load_seeds(seed_path)
    missing = [seed for seed in seeds if seed not in graph]
    if missing:
        raise SystemExit(f"Seed nodes not present in the graph: {', '.join(missing)}")

    personalization = {node: 0.0 for node in graph.nodes()}
    for seed in seeds:
        personalization[seed] = 1.0 / len(seeds)

    scores = nx.pagerank(graph, alpha=alpha, personalization=personalization, weight="weight")
    ranked_nodes = [
        {"node": node, "score": round(float(score), 6), "is_seed": node in seeds}
        for node, score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    ]
    top_non_seed_nodes = [item for item in ranked_nodes if not item["is_seed"]][:top_k]
    return {
        "input_path": str(input_path),
        "seed_path": str(seed_path),
        "alpha": alpha,
        "seed_nodes": seeds,
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "top_k": top_k,
        "top_non_seed_nodes": top_non_seed_nodes,
        "top_nodes": ranked_nodes[:top_k],
        "score_sum": round(sum(scores.values()), 6),
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Weighted TSV edge list.")
    parser.add_argument("--seeds", type=Path, default=DEFAULT_SEEDS, help="Text file with one seed node per line.")
    parser.add_argument("--alpha", type=float, default=0.85, help="Personalized PageRank damping factor.")
    parser.add_argument("--top-k", type=int, default=5, help="How many nodes to keep in the ranked summaries.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input edge list not found: {args.input}")
    if not args.seeds.exists():
        raise SystemExit(f"Seed file not found: {args.seeds}")
    if not 0 < args.alpha < 1:
        raise SystemExit("--alpha must be between 0 and 1.")
    if args.top_k < 1:
        raise SystemExit("--top-k must be positive.")

    payload = summarize_propagation(args.input.resolve(), args.seeds.resolve(), alpha=args.alpha, top_k=args.top_k)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
