#!/usr/bin/env python3
"""Construct a small NetworkX graph and summarize key topology metrics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "examples" / "toy_pathway_edges.tsv"


def load_networkx():
    try:
        import networkx as nx
    except ImportError as exc:  # pragma: no cover - exercised through runtime tests
        raise SystemExit("This script requires networkx in the active Python environment.") from exc
    return nx


def load_edges(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"No edges found in {path}")
    return rows


def build_graph(edge_rows: list[dict]):
    nx = load_networkx()
    graph = nx.Graph()
    for row in edge_rows:
        graph.add_edge(
            row["source"],
            row["target"],
            interaction=row.get("interaction") or "unknown",
            weight=float(row.get("weight") or 1.0),
        )
    return graph


def summarize_graph(input_path: Path, source_node: str, target_node: str) -> dict:
    nx = load_networkx()
    graph = build_graph(load_edges(input_path))
    degree_centrality = nx.degree_centrality(graph)
    top_degree = sorted(
        (
            {"node": node, "degree_centrality": round(float(score), 6), "degree": int(graph.degree(node))}
            for node, score in degree_centrality.items()
        ),
        key=lambda item: (-item["degree_centrality"], item["node"]),
    )[:3]
    shortest_path = nx.shortest_path(graph, source=source_node, target=target_node)
    return {
        "input_path": str(input_path),
        "graph_type": "undirected",
        "node_count": graph.number_of_nodes(),
        "edge_count": graph.number_of_edges(),
        "connected_component_count": nx.number_connected_components(graph),
        "nodes": sorted(graph.nodes()),
        "top_degree_centrality": top_degree,
        "shortest_path_query": {"source": source_node, "target": target_node},
        "shortest_path_nodes": shortest_path,
        "shortest_path_length": len(shortest_path) - 1,
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
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="TSV edge list path.")
    parser.add_argument("--source-node", default="EGFR", help="Source node for the shortest-path example.")
    parser.add_argument("--target-node", default="STAT3", help="Target node for the shortest-path example.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if not args.input.exists():
        raise SystemExit(f"Input edge list not found: {args.input}")

    payload = summarize_graph(
        input_path=args.input,
        source_node=args.source_node,
        target_node=args.target_node,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
