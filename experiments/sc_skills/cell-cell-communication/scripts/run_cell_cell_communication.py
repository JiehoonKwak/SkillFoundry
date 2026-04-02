#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path
from typing import Any


SKILL_DIR = Path(__file__).resolve().parents[1]
RESULT_METHOD = "starter_consensus"
MARKDOWN_TITLE = "Cell-Cell Communication Starter Report"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata(skill_dir: Path) -> dict[str, Any]:
    return load_json(skill_dir / "metadata.yaml")


def require_keys(mapping: dict[str, Any], required: list[str], *, name: str) -> None:
    missing = [key for key in required if key not in mapping]
    if missing:
        raise AssertionError(f"Missing keys in {name}: {missing}")


def normalize_float(value: Any) -> float:
    return float(value)


def format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row[field]) for field in fieldnames})


def write_markdown(path: Path, title: str, sections: list[dict[str, Any]]) -> None:
    lines = [f"# {title}", ""]
    for section in sections:
        lines.append(f"## {section['name']}")
        lines.append("")
        for bullet in section.get("bullets", []):
            lines.append(f"- {bullet}")
        for paragraph in section.get("paragraphs", []):
            lines.append(paragraph)
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def validate_input(payload: dict[str, Any]) -> None:
    require_keys(
        payload,
        ["run_label", "expression_threshold", "null_permutations", "random_seed", "priority_top_n", "cells", "ligand_receptor_pairs"],
        name="toy input",
    )
    cells = payload["cells"]
    if len(cells) < 4:
        raise AssertionError("The starter expects at least four cells.")
    cell_ids = [cell["cell_id"] for cell in cells]
    if len(cell_ids) != len(set(cell_ids)):
        raise AssertionError("Cell IDs must be unique.")
    groups = {cell["group"] for cell in cells}
    if len(groups) < 2:
        raise AssertionError("The starter expects at least two groups.")
    for cell in cells:
        require_keys(cell, ["cell_id", "group", "expression"], name=f"cell {cell.get('cell_id', '<unknown>')}")
        if not isinstance(cell["expression"], dict) or not cell["expression"]:
            raise AssertionError(f"Cell {cell['cell_id']} must contain a non-empty expression mapping.")
    if not payload["ligand_receptor_pairs"]:
        raise AssertionError("At least one ligand-receptor pair is required.")
    genes = collect_genes(cells)
    for pair in payload["ligand_receptor_pairs"]:
        require_keys(pair, ["pair_id", "ligand", "receptor", "evidence"], name="ligand_receptor_pairs[]")
        if pair["ligand"] not in genes:
            raise AssertionError(f"Missing ligand gene in toy input: {pair['ligand']}")
        if pair["receptor"] not in genes:
            raise AssertionError(f"Missing receptor gene in toy input: {pair['receptor']}")
    adjacency = payload.get("adjacency", [])
    known_cells = set(cell_ids)
    for edge in adjacency:
        if len(edge) != 2:
            raise AssertionError(f"Adjacency edges must contain exactly two cell IDs: {edge}")
        if edge[0] not in known_cells or edge[1] not in known_cells:
            raise AssertionError(f"Adjacency references unknown cells: {edge}")


def collect_genes(cells: list[dict[str, Any]]) -> list[str]:
    genes: set[str] = set()
    for cell in cells:
        genes.update(cell["expression"].keys())
    return sorted(genes)


def cell_expression(cell: dict[str, Any], gene: str) -> float:
    return normalize_float(cell["expression"].get(gene, 0.0))


def assignment_from_cells(cells: list[dict[str, Any]]) -> dict[str, str]:
    return {cell["cell_id"]: str(cell["group"]) for cell in cells}


def compute_group_stats(cells: list[dict[str, Any]], group_assignment: dict[str, str], genes: list[str]) -> dict[str, dict[str, dict[str, float]]]:
    grouped_cells: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        group = group_assignment[cell["cell_id"]]
        grouped_cells.setdefault(group, []).append(cell)

    stats: dict[str, dict[str, dict[str, float]]] = {}
    for group, group_cells in sorted(grouped_cells.items()):
        cell_count = len(group_cells)
        stats[group] = {}
        for gene in genes:
            values = [cell_expression(cell, gene) for cell in group_cells]
            mean_expression = sum(values) / float(cell_count)
            expression_fraction = sum(value > 0 for value in values) / float(cell_count)
            stats[group][gene] = {
                "cell_count": float(cell_count),
                "mean_expression": mean_expression,
                "expression_fraction": expression_fraction,
            }
    return stats


def build_group_expression_rows(group_stats: dict[str, dict[str, dict[str, float]]], genes: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for group in sorted(group_stats):
        for gene in genes:
            stats = group_stats[group][gene]
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "cell_count": int(stats["cell_count"]),
                    "mean_expression": stats["mean_expression"],
                    "expression_fraction": stats["expression_fraction"],
                }
            )
    return rows


def row_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row["source_group"]),
        str(row["target_group"]),
        str(row["ligand"]),
        str(row["receptor"]),
    )


def compute_pair_rows(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    group_assignment: dict[str, str],
    genes: list[str],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, dict[str, float]]]]:
    threshold = normalize_float(payload["expression_threshold"])
    group_stats = compute_group_stats(cells, group_assignment, genes)
    groups = sorted(group_stats)
    rows: list[dict[str, Any]] = []
    for pair in payload["ligand_receptor_pairs"]:
        ligand = str(pair["ligand"])
        receptor = str(pair["receptor"])
        for source_group in groups:
            ligand_stats = group_stats[source_group][ligand]
            for target_group in groups:
                receptor_stats = group_stats[target_group][receptor]
                expression_pass = (
                    ligand_stats["expression_fraction"] >= threshold
                    and receptor_stats["expression_fraction"] >= threshold
                )
                min_score = min(ligand_stats["mean_expression"], receptor_stats["mean_expression"]) if expression_pass else 0.0
                geometric_mean = math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"]) if expression_pass else 0.0
                cpdb_mean_stat = ((ligand_stats["mean_expression"] + receptor_stats["mean_expression"]) / 2.0) if expression_pass else 0.0
                rows.append(
                    {
                        "pair_id": pair["pair_id"],
                        "source_group": source_group,
                        "target_group": target_group,
                        "ligand": ligand,
                        "receptor": receptor,
                        "evidence": pair["evidence"],
                        "expression_pass": expression_pass,
                        "ligand_mean": ligand_stats["mean_expression"],
                        "receptor_mean": receptor_stats["mean_expression"],
                        "ligand_fraction": ligand_stats["expression_fraction"],
                        "receptor_fraction": receptor_stats["expression_fraction"],
                        "liana_min_score": min_score,
                        "liana_geometric_mean": geometric_mean,
                        "cpdb_mean_stat": cpdb_mean_stat,
                    }
                )
    return rows, group_stats


def assign_rank_scores(rows: list[dict[str, Any]], field: str, output_field: str) -> None:
    ranked = [row for row in rows if row["expression_pass"]]
    ranked.sort(key=lambda row: (-normalize_float(row[field]), row_key(row)))
    total = len(ranked)
    if total == 0:
        raise AssertionError("No interactions passed the expression threshold.")
    for row in rows:
        row[output_field] = 0.0
    if total == 1:
        ranked[0][output_field] = 1.0
        return
    for position, row in enumerate(ranked, start=1):
        row[output_field] = 1.0 - ((position - 1) / float(total - 1))


def generate_permuted_assignments(cells: list[dict[str, Any]], permutations: int, seed: int) -> list[dict[str, str]]:
    labels = [str(cell["group"]) for cell in cells]
    baseline = tuple(labels)
    rng = random.Random(seed)
    seen: set[tuple[str, ...]] = set()
    assignments: list[dict[str, str]] = []
    max_attempts = max(50, permutations * 20)
    attempts = 0
    while len(assignments) < permutations and attempts < max_attempts:
        attempts += 1
        shuffled = labels[:]
        rng.shuffle(shuffled)
        candidate = tuple(shuffled)
        if candidate == baseline or candidate in seen:
            continue
        seen.add(candidate)
        assignments.append({cells[index]["cell_id"]: shuffled[index] for index in range(len(cells))})
    if len(assignments) != permutations:
        raise AssertionError("Could not generate the requested number of unique deterministic shuffles.")
    return assignments


def summarize_null_statistics(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    genes: list[str],
    observed_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    null_stats: dict[tuple[str, str, str, str], list[float]] = {row_key(row): [] for row in observed_rows}
    for assignment in generate_permuted_assignments(cells, int(payload["null_permutations"]), int(payload["random_seed"])):
        permuted_rows, _ = compute_pair_rows(payload, cells, assignment, genes)
        for row in permuted_rows:
            null_stats[row_key(row)].append(normalize_float(row["cpdb_mean_stat"]))

    summary_rows: list[dict[str, Any]] = []
    permutation_count = int(payload["null_permutations"])
    for row in observed_rows:
        stats = null_stats[row_key(row)]
        observed_stat = normalize_float(row["cpdb_mean_stat"])
        ge_count = sum(value >= observed_stat for value in stats)
        empirical_pvalue = (ge_count + 1) / float(permutation_count + 1)
        row["null_mean_stat"] = sum(stats) / float(len(stats))
        row["null_max_stat"] = max(stats)
        row["null_ge_count"] = ge_count
        row["cpdb_empirical_pvalue"] = empirical_pvalue
        summary_rows.append(
            {
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "observed_mean_stat": observed_stat,
                "null_mean_stat": row["null_mean_stat"],
                "null_max_stat": row["null_max_stat"],
                "null_ge_count": ge_count,
                "empirical_pvalue": empirical_pvalue,
            }
        )
    summary_rows.sort(key=lambda item: (normalize_float(item["empirical_pvalue"]), row_key(item)))
    return summary_rows


def compute_spatial_support(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    group_assignment: dict[str, str],
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    cells_by_id = {cell["cell_id"]: cell for cell in cells}
    adjacency = payload.get("adjacency", [])
    summary_rows: list[dict[str, Any]] = []
    if not adjacency:
        for row in rows:
            row["group_edges"] = 0
            row["supporting_edges"] = 0
            row["adjacency_fraction"] = 0.0
            row["spatial_support"] = "not_tested"
            summary_rows.append(
                {
                    "source_group": row["source_group"],
                    "target_group": row["target_group"],
                    "ligand": row["ligand"],
                    "receptor": row["receptor"],
                    "supporting_edges": 0,
                    "group_edges": 0,
                    "adjacency_fraction": 0.0,
                    "spatial_support": "not_tested",
                }
            )
        return summary_rows

    for row in rows:
        group_edges = 0
        supporting_edges = 0
        for first_cell_id, second_cell_id in adjacency:
            for source_id, target_id in ((first_cell_id, second_cell_id), (second_cell_id, first_cell_id)):
                if group_assignment[source_id] != row["source_group"] or group_assignment[target_id] != row["target_group"]:
                    continue
                group_edges += 1
                source_cell = cells_by_id[source_id]
                target_cell = cells_by_id[target_id]
                if cell_expression(source_cell, row["ligand"]) > 0 and cell_expression(target_cell, row["receptor"]) > 0:
                    supporting_edges += 1
        adjacency_fraction = supporting_edges / float(group_edges) if group_edges else 0.0
        if group_edges == 0:
            spatial_support = "not_observed"
        elif supporting_edges >= 2 and adjacency_fraction >= 0.5:
            spatial_support = "supported"
        elif supporting_edges > 0:
            spatial_support = "limited"
        else:
            spatial_support = "not_observed"
        row["group_edges"] = group_edges
        row["supporting_edges"] = supporting_edges
        row["adjacency_fraction"] = adjacency_fraction
        row["spatial_support"] = spatial_support
        summary_rows.append(
            {
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "supporting_edges": supporting_edges,
                "group_edges": group_edges,
                "adjacency_fraction": adjacency_fraction,
                "spatial_support": spatial_support,
            }
        )
    summary_rows.sort(key=lambda item: (-normalize_float(item["adjacency_fraction"]), row_key(item)))
    return summary_rows


def priority_reason(row: dict[str, Any]) -> str:
    reasons: list[str] = []
    if normalize_float(row["liana_rank_aggregate"]) >= 0.9:
        reasons.append("top LIANA-like aggregate")
    if normalize_float(row["cpdb_empirical_pvalue"]) <= 0.125:
        reasons.append("strong shuffle null separation")
    elif normalize_float(row["cpdb_empirical_pvalue"]) <= 0.25:
        reasons.append("moderate shuffle null separation")
    if row["spatial_support"] == "supported":
        reasons.append("adjacent expressing neighborhoods")
    elif row["spatial_support"] == "limited":
        reasons.append("some adjacency support")
    if row["evidence"] == "positive_control":
        reasons.append("positive-control catalog pair")
    if not reasons:
        reasons.append("passes expression filter")
    return "; ".join(reasons)


def finalize_scores(rows: list[dict[str, Any]]) -> None:
    assign_rank_scores(rows, "liana_min_score", "liana_min_rank")
    assign_rank_scores(rows, "liana_geometric_mean", "liana_geometric_rank")
    for row in rows:
        if not row["expression_pass"]:
            row["liana_rank_aggregate"] = 0.0
            row["score"] = 0.0
            row["method"] = RESULT_METHOD
            continue
        row["liana_rank_aggregate"] = (normalize_float(row["liana_min_rank"]) + normalize_float(row["liana_geometric_rank"])) / 2.0
        score = (
            0.55 * normalize_float(row["liana_rank_aggregate"])
            + 0.30 * (1.0 - normalize_float(row["cpdb_empirical_pvalue"]))
            + 0.15 * normalize_float(row["adjacency_fraction"])
        )
        row["score"] = score
        row["method"] = RESULT_METHOD


def sort_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            -normalize_float(row["score"]),
            normalize_float(row["cpdb_empirical_pvalue"]),
            row_key(row),
        ),
    )


def build_priority_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, Any]]:
    ranked_rows = [row for row in sort_rows(rows) if row["expression_pass"]][:top_n]
    priority_rows: list[dict[str, Any]] = []
    for rank, row in enumerate(ranked_rows, start=1):
        priority_rows.append(
            {
                "rank": rank,
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "priority_reason": priority_reason(row),
                "score": row["score"],
                "liana_rank_aggregate": row["liana_rank_aggregate"],
                "cpdb_empirical_pvalue": row["cpdb_empirical_pvalue"],
                "spatial_support": row["spatial_support"],
            }
        )
    return priority_rows


def build_communication_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ordered_rows = sort_rows(rows)
    deliverable_rows: list[dict[str, Any]] = []
    for row in ordered_rows:
        deliverable_rows.append(
            {
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "score": row["score"],
                "method": row["method"],
                "spatial_support": row["spatial_support"],
                "expression_pass": row["expression_pass"],
                "ligand_mean": row["ligand_mean"],
                "receptor_mean": row["receptor_mean"],
                "ligand_fraction": row["ligand_fraction"],
                "receptor_fraction": row["receptor_fraction"],
                "liana_min_score": row["liana_min_score"],
                "liana_geometric_mean": row["liana_geometric_mean"],
                "liana_rank_aggregate": row["liana_rank_aggregate"],
                "cpdb_mean_stat": row["cpdb_mean_stat"],
                "cpdb_empirical_pvalue": row["cpdb_empirical_pvalue"],
                "supporting_edges": row["supporting_edges"],
                "group_edges": row["group_edges"],
                "adjacency_fraction": row["adjacency_fraction"],
                "evidence": row["evidence"],
            }
        )
    return deliverable_rows


def build_report_sections(
    payload: dict[str, Any],
    communication_rows: list[dict[str, Any]],
    priority_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    supported_pairs = [row for row in communication_rows if row["spatial_support"] == "supported"]
    top_pair = priority_rows[0]
    return [
        {
            "name": "Run context",
            "bullets": [
                f"Starter label: {payload['run_label']}.",
                f"Toy cells: {len(payload['cells'])} across {len({cell['group'] for cell in payload['cells']})} groups.",
                f"Ligand-receptor catalog size: {len(payload['ligand_receptor_pairs'])}; deterministic shuffles: {payload['null_permutations']}.",
            ],
        },
        {
            "name": "Methods",
            "bullets": [
                "Computed group means and expression fractions from the raw toy cell matrix.",
                "Approximated LIANA aggregation by averaging rank scores from min-based and geometric-mean interaction magnitudes.",
                "Approximated CellPhoneDB statistics with a deterministic label-shuffle null over the toy groups.",
                "Computed a graph-based spatial cross-check from the provided adjacency list.",
            ],
        },
        {
            "name": "Prioritized interactions",
            "bullets": [
                f"Top pair: {top_pair['ligand']}->{top_pair['receptor']} from {top_pair['source_group']} to {top_pair['target_group']} (score {top_pair['score']:.3f}, empirical p-value {top_pair['cpdb_empirical_pvalue']:.3f}).",
                *[
                    f"Rank {row['rank']}: {row['ligand']}->{row['receptor']} {row['source_group']}->{row['target_group']} with {row['priority_reason']}."
                    for row in priority_rows[1:4]
                ],
            ],
        },
        {
            "name": "Spatial support",
            "bullets": [
                f"{len(supported_pairs)} interaction rows reached the starter's 'supported' adjacency threshold.",
                f"The top pair shows {communication_rows[0]['supporting_edges']} supporting directed neighbor edges out of {communication_rows[0]['group_edges']} group-specific neighbor opportunities.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter uses synthetic expression values, a tiny ligand-receptor catalog, and a toy adjacency graph.",
                "It does not reproduce the full LIANA, CellPhoneDB, or Squidpy method stacks, databases, or multiple-testing procedures.",
                "Use the outputs only as a deterministic shell-friendly starter path and upgrade to real public workflows for biological interpretation.",
            ],
        },
    ]


def write_outputs(
    payload: dict[str, Any],
    input_path: Path,
    outdir: Path,
    communication_rows: list[dict[str, Any]],
    priority_rows: list[dict[str, Any]],
    group_rows: list[dict[str, Any]],
    null_rows: list[dict[str, Any]],
    spatial_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    communication_path = outdir / "communication_results.tsv"
    priority_path = outdir / "priority_pairs.tsv"
    report_path = outdir / "interpretation_report.md"
    group_path = outdir / "qc_group_expression.tsv"
    null_path = outdir / "qc_pair_nulls.tsv"
    spatial_path = outdir / "qc_spatial_support.tsv"

    write_tsv(
        communication_path,
        communication_rows,
        [
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "score",
            "method",
            "spatial_support",
            "expression_pass",
            "ligand_mean",
            "receptor_mean",
            "ligand_fraction",
            "receptor_fraction",
            "liana_min_score",
            "liana_geometric_mean",
            "liana_rank_aggregate",
            "cpdb_mean_stat",
            "cpdb_empirical_pvalue",
            "supporting_edges",
            "group_edges",
            "adjacency_fraction",
            "evidence",
        ],
    )
    write_tsv(
        priority_path,
        priority_rows,
        [
            "rank",
            "ligand",
            "receptor",
            "source_group",
            "target_group",
            "priority_reason",
            "score",
            "liana_rank_aggregate",
            "cpdb_empirical_pvalue",
            "spatial_support",
        ],
    )
    write_tsv(
        group_path,
        group_rows,
        ["group", "gene", "cell_count", "mean_expression", "expression_fraction"],
    )
    write_tsv(
        null_path,
        null_rows,
        [
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "observed_mean_stat",
            "null_mean_stat",
            "null_max_stat",
            "null_ge_count",
            "empirical_pvalue",
        ],
    )
    write_tsv(
        spatial_path,
        spatial_rows,
        [
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "supporting_edges",
            "group_edges",
            "adjacency_fraction",
            "spatial_support",
        ],
    )
    write_markdown(report_path, MARKDOWN_TITLE, build_report_sections(payload, communication_rows, priority_rows))

    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path),
        "counts": {
            "cell_count": len(payload["cells"]),
            "group_count": len({cell["group"] for cell in payload["cells"]}),
            "catalog_pair_count": len(payload["ligand_receptor_pairs"]),
            "communication_rows": len(communication_rows),
            "priority_rows": len(priority_rows),
        },
        "top_pair": {
            "ligand": priority_rows[0]["ligand"],
            "receptor": priority_rows[0]["receptor"],
            "source_group": priority_rows[0]["source_group"],
            "target_group": priority_rows[0]["target_group"],
            "score": round(normalize_float(priority_rows[0]["score"]), 6),
            "empirical_pvalue": round(normalize_float(priority_rows[0]["cpdb_empirical_pvalue"]), 6),
        },
        "configuration": {
            "expression_threshold": payload["expression_threshold"],
            "null_permutations": payload["null_permutations"],
            "random_seed": payload["random_seed"],
            "priority_top_n": payload["priority_top_n"],
        },
        "written_files": sorted(path.name for path in outdir.iterdir() if path.is_file()),
    }
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["written_files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_outputs(payload: dict[str, Any], *, input_path: Path, outdir: Path) -> dict[str, Any]:
    validate_input(payload)
    cells = payload["cells"]
    genes = collect_genes(cells)
    observed_assignment = assignment_from_cells(cells)
    pair_rows, group_stats = compute_pair_rows(payload, cells, observed_assignment, genes)
    null_rows = summarize_null_statistics(payload, cells, genes, pair_rows)
    spatial_rows = compute_spatial_support(payload, cells, observed_assignment, pair_rows)
    finalize_scores(pair_rows)
    communication_rows = build_communication_rows(pair_rows)
    priority_rows = build_priority_rows(pair_rows, int(payload["priority_top_n"]))
    if not priority_rows:
        raise AssertionError("The starter could not prioritize any interaction rows.")
    group_rows = build_group_expression_rows(group_stats, genes)
    return write_outputs(payload, input_path, outdir, communication_rows, priority_rows, group_rows, null_rows, spatial_rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the portable cell-cell communication starter.")
    parser.add_argument("--input", type=Path, required=True, help="Path to the raw toy or project input JSON.")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory where starter outputs will be written.")
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    summary = build_outputs(payload, input_path=args.input, outdir=args.outdir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
