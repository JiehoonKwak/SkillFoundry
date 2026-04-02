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
REPORT_TITLE = "CellPhoneDB Portable Starter Report"
METHOD_LABEL = "portable_cpdb_starter"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def require_keys(mapping: dict[str, Any], required: list[str], *, name: str) -> None:
    missing = [key for key in required if key not in mapping]
    if missing:
        raise AssertionError(f"Missing keys in {name}: {missing}")


def normalize_float(value: Any) -> float:
    return float(value)


def format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
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


def collect_genes(cells: list[dict[str, Any]]) -> list[str]:
    genes: set[str] = set()
    for cell in cells:
        genes.update(cell["expression"].keys())
    return sorted(genes)


def validate_input(payload: dict[str, Any]) -> None:
    require_keys(
        payload,
        [
            "run_label",
            "expression_threshold",
            "null_permutations",
            "pvalue_cutoff",
            "random_seed",
            "ranking_top_n",
            "cells",
            "ligand_receptor_pairs",
            "degs",
        ],
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

    genes = collect_genes(cells)
    if not payload["ligand_receptor_pairs"]:
        raise AssertionError("At least one ligand-receptor pair is required.")
    for pair in payload["ligand_receptor_pairs"]:
        require_keys(pair, ["pair_id", "ligand", "receptor", "evidence"], name="ligand_receptor_pairs[]")
        if pair["ligand"] not in genes:
            raise AssertionError(f"Missing ligand gene in toy input: {pair['ligand']}")
        if pair["receptor"] not in genes:
            raise AssertionError(f"Missing receptor gene in toy input: {pair['receptor']}")

    deg_entries = payload["degs"]
    if not deg_entries:
        raise AssertionError("At least one DEG entry is required for the starter DEG mode.")
    for entry in deg_entries:
        require_keys(entry, ["group", "gene"], name="degs[]")
        if entry["group"] not in groups:
            raise AssertionError(f"DEG entry references unknown group: {entry['group']}")
        if entry["gene"] not in genes:
            raise AssertionError(f"DEG entry references unknown gene: {entry['gene']}")

    adjacency = payload.get("adjacency", [])
    known_cells = set(cell_ids)
    for edge in adjacency:
        if len(edge) != 2:
            raise AssertionError(f"Adjacency edges must contain exactly two cell IDs: {edge}")
        if edge[0] not in known_cells or edge[1] not in known_cells:
            raise AssertionError(f"Adjacency references unknown cells: {edge}")


def cell_expression(cell: dict[str, Any], gene: str) -> float:
    return normalize_float(cell["expression"].get(gene, 0.0))


def assignment_from_cells(cells: list[dict[str, Any]]) -> dict[str, str]:
    return {cell["cell_id"]: str(cell["group"]) for cell in cells}


def compute_group_stats(
    cells: list[dict[str, Any]],
    assignment: dict[str, str],
    genes: list[str],
) -> dict[str, dict[str, dict[str, float]]]:
    grouped_cells: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        grouped_cells.setdefault(assignment[cell["cell_id"]], []).append(cell)

    stats: dict[str, dict[str, dict[str, float]]] = {}
    for group, group_cells in sorted(grouped_cells.items()):
        stats[group] = {}
        cell_count = len(group_cells)
        for gene in genes:
            values = [cell_expression(cell, gene) for cell in group_cells]
            stats[group][gene] = {
                "cell_count": float(cell_count),
                "mean_expression": sum(values) / float(cell_count),
                "expression_fraction": sum(value > 0 for value in values) / float(cell_count),
            }
    return stats


def build_group_rows(
    group_stats: dict[str, dict[str, dict[str, float]]],
    genes: list[str],
) -> list[dict[str, Any]]:
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


def deg_lookup(payload: dict[str, Any]) -> dict[str, set[str]]:
    lookup: dict[str, set[str]] = {}
    for entry in payload["degs"]:
        lookup.setdefault(str(entry["group"]), set()).add(str(entry["gene"]))
    return lookup


def pair_row_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (str(row["source_group"]), str(row["target_group"]), str(row["interacting_pair"]))


def compute_pair_rows(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    assignment: dict[str, str],
    genes: list[str],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, dict[str, float]]]]:
    threshold = normalize_float(payload["expression_threshold"])
    group_stats = compute_group_stats(cells, assignment, genes)
    deg_sets = deg_lookup(payload)
    rows: list[dict[str, Any]] = []
    groups = sorted(group_stats)

    for pair in payload["ligand_receptor_pairs"]:
        ligand = str(pair["ligand"])
        receptor = str(pair["receptor"])
        interacting_pair = f"{ligand}_{receptor}"
        for source_group in groups:
            ligand_stats = group_stats[source_group][ligand]
            for target_group in groups:
                receptor_stats = group_stats[target_group][receptor]
                expression_pass = (
                    ligand_stats["expression_fraction"] >= threshold
                    and receptor_stats["expression_fraction"] >= threshold
                )
                mean_expression = (
                    (ligand_stats["mean_expression"] + receptor_stats["mean_expression"]) / 2.0
                    if expression_pass
                    else 0.0
                )
                min_expression = (
                    min(ligand_stats["mean_expression"], receptor_stats["mean_expression"])
                    if expression_pass
                    else 0.0
                )
                geometric_mean = (
                    math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"])
                    if expression_pass
                    else 0.0
                )
                ligand_is_deg = ligand in deg_sets.get(source_group, set())
                receptor_is_deg = receptor in deg_sets.get(target_group, set())
                deg_support_score = (0.5 if ligand_is_deg else 0.0) + (0.5 if receptor_is_deg else 0.0)
                rows.append(
                    {
                        "pair_id": pair["pair_id"],
                        "interacting_pair": interacting_pair,
                        "ligand": ligand,
                        "receptor": receptor,
                        "source_group": source_group,
                        "target_group": target_group,
                        "evidence": pair["evidence"],
                        "expression_pass": expression_pass,
                        "ligand_mean": ligand_stats["mean_expression"],
                        "receptor_mean": receptor_stats["mean_expression"],
                        "ligand_fraction": ligand_stats["expression_fraction"],
                        "receptor_fraction": receptor_stats["expression_fraction"],
                        "mean_expression": mean_expression,
                        "min_expression": min_expression,
                        "geometric_mean": geometric_mean,
                        "ligand_is_deg": ligand_is_deg,
                        "receptor_is_deg": receptor_is_deg,
                        "deg_relevant": expression_pass and (ligand_is_deg or receptor_is_deg),
                        "deg_support_score": deg_support_score if expression_pass else 0.0,
                    }
                )
    return rows, group_stats


def assign_rank_scores(rows: list[dict[str, Any]], field: str, output_field: str) -> None:
    ranked_rows = [row for row in rows if row["expression_pass"]]
    ranked_rows.sort(key=lambda row: (-normalize_float(row[field]), pair_row_key(row)))
    if not ranked_rows:
        raise AssertionError("No interactions passed the expression threshold.")

    for row in rows:
        row[output_field] = 0.0
    if len(ranked_rows) == 1:
        ranked_rows[0][output_field] = 1.0
        return

    total = len(ranked_rows)
    for position, row in enumerate(ranked_rows, start=1):
        row[output_field] = 1.0 - ((position - 1) / float(total - 1))


def generate_permuted_assignments(cells: list[dict[str, Any]], permutations: int, seed: int) -> list[dict[str, str]]:
    labels = [str(cell["group"]) for cell in cells]
    baseline = tuple(labels)
    rng = random.Random(seed)
    seen: set[tuple[str, ...]] = set()
    assignments: list[dict[str, str]] = []
    attempts = 0
    max_attempts = max(100, permutations * 40)

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
    null_stats: dict[tuple[str, str, str], list[float]] = {pair_row_key(row): [] for row in observed_rows}
    for assignment in generate_permuted_assignments(cells, int(payload["null_permutations"]), int(payload["random_seed"])):
        permuted_rows, _ = compute_pair_rows(payload, cells, assignment, genes)
        for row in permuted_rows:
            null_stats[pair_row_key(row)].append(normalize_float(row["mean_expression"]))

    permutation_count = int(payload["null_permutations"])
    summary_rows: list[dict[str, Any]] = []
    for row in observed_rows:
        observed = normalize_float(row["mean_expression"])
        stats = null_stats[pair_row_key(row)]
        ge_count = sum(value >= observed for value in stats)
        empirical_pvalue = (ge_count + 1) / float(permutation_count + 1)
        row["null_mean_expression"] = sum(stats) / float(len(stats))
        row["null_max_expression"] = max(stats)
        row["null_ge_count"] = ge_count
        row["empirical_pvalue"] = empirical_pvalue
        row["is_significant"] = row["expression_pass"] and empirical_pvalue <= normalize_float(payload["pvalue_cutoff"])
        summary_rows.append(
            {
                "interacting_pair": row["interacting_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "observed_mean_expression": observed,
                "null_mean_expression": row["null_mean_expression"],
                "null_max_expression": row["null_max_expression"],
                "null_ge_count": ge_count,
                "empirical_pvalue": empirical_pvalue,
            }
        )

    summary_rows.sort(
        key=lambda row: (
            normalize_float(row["empirical_pvalue"]),
            -normalize_float(row["observed_mean_expression"]),
            pair_row_key(row),
        )
    )
    return summary_rows


def summarize_deg_support(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary_rows: list[dict[str, Any]] = []
    for row in rows:
        summary_rows.append(
            {
                "interacting_pair": row["interacting_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand_is_deg": row["ligand_is_deg"],
                "receptor_is_deg": row["receptor_is_deg"],
                "deg_relevant": row["deg_relevant"],
                "deg_support_score": row["deg_support_score"],
            }
        )
    summary_rows.sort(key=pair_row_key)
    return summary_rows


def compute_adjacency_support(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    assignment: dict[str, str],
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
        return summarize_adjacency_rows(rows)

    for row in rows:
        group_edges = 0
        supporting_edges = 0
        for first_cell_id, second_cell_id in adjacency:
            for source_id, target_id in ((first_cell_id, second_cell_id), (second_cell_id, first_cell_id)):
                if assignment[source_id] != row["source_group"] or assignment[target_id] != row["target_group"]:
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

    summary_rows.extend(summarize_adjacency_rows(rows))
    return summary_rows


def summarize_adjacency_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary_rows: list[dict[str, Any]] = []
    for row in rows:
        summary_rows.append(
            {
                "interacting_pair": row["interacting_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "supporting_edges": row["supporting_edges"],
                "group_edges": row["group_edges"],
                "adjacency_fraction": row["adjacency_fraction"],
                "spatial_support": row["spatial_support"],
            }
        )
    summary_rows.sort(
        key=lambda row: (
            -normalize_float(row["adjacency_fraction"]),
            pair_row_key(row),
        )
    )
    return summary_rows


def compute_specificity_rank(rows: list[dict[str, Any]]) -> None:
    total_group_pairs = len({(row["source_group"], row["target_group"]) for row in rows})
    significant_counts: dict[str, int] = {}
    for row in rows:
        if row["is_significant"]:
            significant_counts[row["interacting_pair"]] = significant_counts.get(row["interacting_pair"], 0) + 1
    for row in rows:
        row["cpdb_specificity_rank"] = significant_counts.get(row["interacting_pair"], 0) / float(total_group_pairs)


def finalize_scores(rows: list[dict[str, Any]]) -> None:
    assign_rank_scores(rows, "min_expression", "min_expression_rank")
    assign_rank_scores(rows, "geometric_mean", "geometric_mean_rank")
    compute_specificity_rank(rows)

    for row in rows:
        if not row["expression_pass"]:
            row["score"] = 0.0
            row["method"] = METHOD_LABEL
            continue
        row["score"] = (
            0.30 * normalize_float(row["min_expression_rank"])
            + 0.20 * normalize_float(row["geometric_mean_rank"])
            + 0.20 * (1.0 - normalize_float(row["empirical_pvalue"]))
            + 0.10 * normalize_float(row["deg_support_score"])
            + 0.10 * normalize_float(row["adjacency_fraction"])
            + 0.10 * (1.0 - normalize_float(row["cpdb_specificity_rank"]))
        )
        row["method"] = METHOD_LABEL


def sort_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        [row for row in rows if row["expression_pass"]],
        key=lambda row: (
            -int(bool(row["is_significant"])),
            -normalize_float(row["score"]),
            normalize_float(row["empirical_pvalue"]),
            normalize_float(row["cpdb_specificity_rank"]),
            pair_row_key(row),
        ),
    )


def evidence_string(row: dict[str, Any]) -> str:
    reasons: list[str] = []
    if row["is_significant"]:
        reasons.append("significant_mean")
    else:
        reasons.append("expression_pass")
    if row["ligand_is_deg"] and row["receptor_is_deg"]:
        reasons.append("dual_deg_support")
    elif row["deg_relevant"]:
        reasons.append("deg_supported")
    if row["spatial_support"] == "supported":
        reasons.append("adjacency_supported")
    elif row["spatial_support"] == "limited":
        reasons.append("adjacency_limited")
    reasons.append(f"empirical_p={row['empirical_pvalue']:.3f}")
    return " + ".join(reasons)


def build_significant_means(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    significant_rows = [row for row in rows if row["is_significant"]]
    significant_rows.sort(
        key=lambda row: (
            normalize_float(row["empirical_pvalue"]),
            -normalize_float(row["mean_expression"]),
            pair_row_key(row),
        )
    )
    deliverable_rows: list[dict[str, Any]] = []
    for row in significant_rows:
        deliverable_rows.append(
            {
                "interacting_pair": row["interacting_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "mean_expression": row["mean_expression"],
                "pvalue": row["empirical_pvalue"],
                "rank": row["cpdb_specificity_rank"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "deg_relevant": row["deg_relevant"],
                "spatial_support": row["spatial_support"],
                "adjacency_fraction": row["adjacency_fraction"],
            }
        )
    return deliverable_rows


def build_ranked_interactions(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, Any]]:
    ranked_rows = sort_rows(rows)[:top_n]
    deliverable_rows: list[dict[str, Any]] = []
    for rank, row in enumerate(ranked_rows, start=1):
        deliverable_rows.append(
            {
                "rank": rank,
                "interacting_pair": row["interacting_pair"],
                "evidence": evidence_string(row),
                "supporting_group_pair": f"{row['source_group']}->{row['target_group']}",
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "mean_expression": row["mean_expression"],
                "pvalue": row["empirical_pvalue"],
                "score": row["score"],
                "deg_relevant": row["deg_relevant"],
                "spatial_support": row["spatial_support"],
                "cpdb_specificity_rank": row["cpdb_specificity_rank"],
                "method": row["method"],
            }
        )
    return deliverable_rows


def build_report_sections(
    payload: dict[str, Any],
    significant_rows: list[dict[str, Any]],
    ranked_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    top_ranked = ranked_rows[0]
    significant_preview = significant_rows[:3]
    if significant_preview:
        significant_bullets = [
            (
                f"{row['interacting_pair']} in {row['source_group']}->{row['target_group']} "
                f"with mean {row['mean_expression']:.3f}, empirical p-value {row['pvalue']:.3f}, "
                f"and spatial status {row['spatial_support']}."
            )
            for row in significant_preview
        ]
    else:
        significant_bullets = ["No rows crossed the configured empirical p-value cutoff in this toy run."]

    return [
        {
            "name": "Run context",
            "bullets": [
                f"Starter label: {payload['run_label']}.",
                f"Toy cells: {len(payload['cells'])} across {len({cell['group'] for cell in payload['cells']})} groups.",
                (
                    f"Ligand-receptor catalog size: {len(payload['ligand_receptor_pairs'])}; "
                    f"DEG entries: {len(payload['degs'])}; deterministic shuffles: {payload['null_permutations']}."
                ),
            ],
        },
        {
            "name": "Significant interactions",
            "bullets": significant_bullets,
        },
        {
            "name": "Ranking logic",
            "bullets": [
                "Rows first pass an expression-fraction filter based on ligand and receptor detection across the source and target groups.",
                "Empirical p-values come from a deterministic cell-label shuffle null over the toy groups, mirroring CellPhoneDB statistical mode at tiny scale.",
                "Ranking combines min-expression, geometric mean, empirical p-value, DEG relevance, adjacency support, and interaction specificity.",
                (
                    f"Top ranked interaction: {top_ranked['interacting_pair']} in {top_ranked['supporting_group_pair']} "
                    f"with score {top_ranked['score']:.3f} and empirical p-value {top_ranked['pvalue']:.3f}."
                ),
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter uses synthetic expression values, a tiny ligand-receptor catalog, and a toy adjacency graph.",
                "It does not reproduce the full CellPhoneDB database, heteromer handling, production permutation depth, or FDR procedures.",
                "Use the outputs only as a deterministic shell-friendly starter path and upgrade to the real public workflow for biological interpretation.",
            ],
        },
    ]


def write_outputs(
    payload: dict[str, Any],
    input_path: Path,
    outdir: Path,
    significant_rows: list[dict[str, Any]],
    ranked_rows: list[dict[str, Any]],
    group_rows: list[dict[str, Any]],
    null_rows: list[dict[str, Any]],
    deg_rows: list[dict[str, Any]],
    adjacency_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    significant_path = outdir / "cellphonedb_significant_means.tsv"
    ranked_path = outdir / "interaction_ranked.tsv"
    report_path = outdir / "cellphonedb_report.md"
    group_path = outdir / "qc_group_expression.tsv"
    null_path = outdir / "qc_statistical_null.tsv"
    deg_path = outdir / "qc_deg_support.tsv"
    adjacency_path = outdir / "qc_adjacency_support.tsv"

    write_tsv(
        significant_path,
        significant_rows,
        [
            "interacting_pair",
            "source_group",
            "target_group",
            "mean_expression",
            "pvalue",
            "rank",
            "ligand",
            "receptor",
            "deg_relevant",
            "spatial_support",
            "adjacency_fraction",
        ],
    )
    write_tsv(
        ranked_path,
        ranked_rows,
        [
            "rank",
            "interacting_pair",
            "evidence",
            "supporting_group_pair",
            "source_group",
            "target_group",
            "mean_expression",
            "pvalue",
            "score",
            "deg_relevant",
            "spatial_support",
            "cpdb_specificity_rank",
            "method",
        ],
    )
    write_tsv(group_path, group_rows, ["group", "gene", "cell_count", "mean_expression", "expression_fraction"])
    write_tsv(
        null_path,
        null_rows,
        [
            "interacting_pair",
            "source_group",
            "target_group",
            "observed_mean_expression",
            "null_mean_expression",
            "null_max_expression",
            "null_ge_count",
            "empirical_pvalue",
        ],
    )
    write_tsv(
        deg_path,
        deg_rows,
        [
            "interacting_pair",
            "source_group",
            "target_group",
            "ligand_is_deg",
            "receptor_is_deg",
            "deg_relevant",
            "deg_support_score",
        ],
    )
    write_tsv(
        adjacency_path,
        adjacency_rows,
        [
            "interacting_pair",
            "source_group",
            "target_group",
            "supporting_edges",
            "group_edges",
            "adjacency_fraction",
            "spatial_support",
        ],
    )
    write_markdown(report_path, REPORT_TITLE, build_report_sections(payload, significant_rows, ranked_rows))

    top_interaction = ranked_rows[0]
    summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path),
        "configuration": {
            "expression_threshold": payload["expression_threshold"],
            "null_permutations": payload["null_permutations"],
            "pvalue_cutoff": payload["pvalue_cutoff"],
            "random_seed": payload["random_seed"],
            "ranking_top_n": payload["ranking_top_n"],
        },
        "counts": {
            "cell_count": len(payload["cells"]),
            "group_count": len({cell["group"] for cell in payload["cells"]}),
            "catalog_pair_count": len(payload["ligand_receptor_pairs"]),
            "significant_rows": len(significant_rows),
            "ranked_rows": len(ranked_rows),
        },
        "top_interaction": {
            "interacting_pair": top_interaction["interacting_pair"],
            "source_group": top_interaction["source_group"],
            "target_group": top_interaction["target_group"],
            "score": round(normalize_float(top_interaction["score"]), 6),
            "pvalue": round(normalize_float(top_interaction["pvalue"]), 6),
        },
        "written_files": [],
    }
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["written_files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_outputs(payload: dict[str, Any], *, input_path: Path, outdir: Path) -> dict[str, Any]:
    validate_input(payload)
    cells = payload["cells"]
    genes = collect_genes(cells)
    assignment = assignment_from_cells(cells)
    pair_rows, group_stats = compute_pair_rows(payload, cells, assignment, genes)
    null_rows = summarize_null_statistics(payload, cells, genes, pair_rows)
    deg_rows = summarize_deg_support(pair_rows)
    adjacency_rows = compute_adjacency_support(payload, cells, assignment, pair_rows)
    finalize_scores(pair_rows)
    significant_rows = build_significant_means(pair_rows)
    ranked_rows = build_ranked_interactions(pair_rows, int(payload["ranking_top_n"]))
    if not significant_rows:
        raise AssertionError("The starter did not produce any significant means rows.")
    if not ranked_rows:
        raise AssertionError("The starter did not produce any ranked interaction rows.")
    group_rows = build_group_rows(group_stats, genes)
    return write_outputs(
        payload,
        input_path,
        outdir,
        significant_rows,
        ranked_rows,
        group_rows,
        null_rows,
        deg_rows,
        adjacency_rows,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the portable CellPhoneDB ligand-receptor starter.")
    parser.add_argument("--input", type=Path, required=True, help="Path to the raw toy or project input JSON.")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory where starter outputs will be written.")
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    summary = build_outputs(payload, input_path=args.input, outdir=args.outdir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
