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
REPORT_TITLE = "LIANA Portable Starter Report"
ALLOWED_METHODS = ("min_expr", "geometric_mean", "specificity", "shuffle_pvalue")
METHOD_SORT_DIRECTION = {
    "min_expr": "desc",
    "geometric_mean": "desc",
    "specificity": "desc",
    "shuffle_pvalue": "asc",
}
METHOD_FIELD_NAME = {
    "min_expr": "min_expr",
    "geometric_mean": "geometric_mean",
    "specificity": "specificity",
    "shuffle_pvalue": "empirical_pvalue",
}
METHOD_RANK_FIELD_NAME = {
    "min_expr": "min_expr_rank",
    "geometric_mean": "geometric_mean_rank",
    "specificity": "specificity_method_rank",
    "shuffle_pvalue": "shuffle_pvalue_rank",
}


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
            writer.writerow({field: format_value(row.get(field, "")) for field in fieldnames})


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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


def canonicalize_alias(value: str, aliases: dict[str, str]) -> str:
    normalized = value.strip()
    lowered = normalized.lower()
    if normalized in aliases:
        return str(aliases[normalized])
    if lowered in aliases:
        return str(aliases[lowered])
    raise AssertionError(f"Missing harmonization alias for value: {value}")


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
            "min_cells",
            "null_permutations",
            "random_seed",
            "ranking_top_n",
            "selected_methods",
            "aggregate_method",
            "harmonization",
            "cells",
            "ligand_receptor_catalog",
            "expected_invariants",
        ],
        name="toy input",
    )
    if payload["aggregate_method"] != "mean_rank":
        raise AssertionError("This starter currently supports only aggregate_method='mean_rank'.")
    selected_methods = payload["selected_methods"]
    if not 2 <= len(selected_methods) <= 4:
        raise AssertionError("The starter expects 2 to 4 selected methods.")
    for method in selected_methods:
        if method not in ALLOWED_METHODS:
            raise AssertionError(f"Unsupported selected method: {method}")

    cells = payload["cells"]
    if len(cells) < 4:
        raise AssertionError("The starter expects at least four cells.")
    cell_ids = [str(cell["cell_id"]) for cell in cells]
    if len(cell_ids) != len(set(cell_ids)):
        raise AssertionError("Cell IDs must be unique.")

    harmonization = payload["harmonization"]
    require_keys(
        harmonization,
        ["group_aliases", "condition_aliases", "resource_aliases", "condition_order"],
        name="harmonization",
    )
    if len(harmonization["condition_order"]) != 2:
        raise AssertionError("The starter expects exactly two ordered conditions for condition-aware summaries.")

    genes = collect_genes(cells)
    for cell in cells:
        require_keys(cell, ["cell_id", "group", "condition", "expression"], name=f"cell {cell.get('cell_id', '<unknown>')}")
        if not isinstance(cell["expression"], dict) or not cell["expression"]:
            raise AssertionError(f"Cell {cell['cell_id']} must contain a non-empty expression mapping.")
        canonicalize_alias(str(cell["group"]), harmonization["group_aliases"])
        canonicalize_alias(str(cell["condition"]), harmonization["condition_aliases"])

    pairs = payload["ligand_receptor_catalog"]
    if not pairs:
        raise AssertionError("At least one ligand-receptor pair is required.")
    for pair in pairs:
        require_keys(pair, ["pair_id", "ligand", "receptor", "resources", "evidence"], name="ligand_receptor_catalog[]")
        if pair["ligand"] not in genes:
            raise AssertionError(f"Missing ligand gene in toy input: {pair['ligand']}")
        if pair["receptor"] not in genes:
            raise AssertionError(f"Missing receptor gene in toy input: {pair['receptor']}")
        for resource in pair["resources"]:
            canonicalize_alias(str(resource), harmonization["resource_aliases"])

    adjacency = payload.get("adjacency", [])
    known_cells = set(cell_ids)
    for edge in adjacency:
        if len(edge) != 2:
            raise AssertionError(f"Adjacency edges must contain exactly two cell IDs: {edge}")
        if edge[0] not in known_cells or edge[1] not in known_cells:
            raise AssertionError(f"Adjacency references unknown cells: {edge}")


def harmonize_cells(payload: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str], list[str]]:
    harmonization = payload["harmonization"]
    cells: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    groups: set[str] = set()
    conditions: set[str] = set()

    for cell in payload["cells"]:
        normalized_group = canonicalize_alias(str(cell["group"]), harmonization["group_aliases"])
        normalized_condition = canonicalize_alias(str(cell["condition"]), harmonization["condition_aliases"])
        copied = dict(cell)
        copied["normalized_group"] = normalized_group
        copied["normalized_condition"] = normalized_condition
        cells.append(copied)
        groups.add(normalized_group)
        conditions.add(normalized_condition)
        qc_rows.append(
            {
                "cell_id": cell["cell_id"],
                "raw_group": cell["group"],
                "normalized_group": normalized_group,
                "raw_condition": cell["condition"],
                "normalized_condition": normalized_condition,
                "harmonized": normalized_group != cell["group"] or normalized_condition != cell["condition"],
            }
        )

    ordered_conditions = [condition for condition in harmonization["condition_order"] if condition in conditions]
    if len(ordered_conditions) != 2:
        raise AssertionError("The toy input must contain both ordered conditions after harmonization.")
    ordered_groups = sorted(groups)
    return cells, qc_rows, ordered_groups, ordered_conditions


def harmonize_catalog(payload: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    resource_aliases = payload["harmonization"]["resource_aliases"]
    pair_specs: list[dict[str, Any]] = []
    resource_rows: list[dict[str, Any]] = []

    for pair in payload["ligand_receptor_catalog"]:
        canonical_resources = sorted({canonicalize_alias(str(resource), resource_aliases) for resource in pair["resources"]})
        interacting_pair = f"{pair['ligand']}_{pair['receptor']}"
        pair_spec = dict(pair)
        pair_spec["pair"] = interacting_pair
        pair_spec["canonical_resources"] = canonical_resources
        pair_spec["resource_count"] = len(canonical_resources)
        pair_specs.append(pair_spec)
        resource_rows.append(
            {
                "pair": interacting_pair,
                "resource_count": len(canonical_resources),
                "resources": ";".join(canonical_resources),
            }
        )

    resource_rows.sort(key=lambda row: (-int(row["resource_count"]), row["pair"]))
    return pair_specs, resource_rows


def expression_value(cell: dict[str, Any], gene: str) -> float:
    return normalize_float(cell["expression"].get(gene, 0.0))


def assignment_from_cells(cells: list[dict[str, Any]]) -> dict[str, str]:
    return {str(cell["cell_id"]): str(cell["normalized_group"]) for cell in cells}


def compute_group_stats(
    cells: list[dict[str, Any]],
    assignment: dict[str, str],
    genes: list[str],
) -> dict[str, dict[str, dict[str, float]]]:
    grouped_cells: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        grouped_cells.setdefault(assignment[str(cell["cell_id"])], []).append(cell)

    stats: dict[str, dict[str, dict[str, float]]] = {}
    for group, group_cells in sorted(grouped_cells.items()):
        stats[group] = {}
        cell_count = len(group_cells)
        for gene in genes:
            values = [expression_value(cell, gene) for cell in group_cells]
            stats[group][gene] = {
                "cell_count": float(cell_count),
                "mean_expression": sum(values) / float(cell_count),
                "expression_fraction": sum(value > 0 for value in values) / float(cell_count),
            }
    return stats


def compute_condition_stats(
    cells: list[dict[str, Any]],
    genes: list[str],
    conditions: list[str],
) -> dict[str, dict[str, dict[str, dict[str, float]]]]:
    stats: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
    for condition in conditions:
        condition_cells = [cell for cell in cells if cell["normalized_condition"] == condition]
        stats[condition] = compute_group_stats(condition_cells, assignment_from_cells(condition_cells), genes)
    return stats


def build_group_rows(
    overall_stats: dict[str, dict[str, dict[str, float]]],
    condition_stats: dict[str, dict[str, dict[str, dict[str, float]]]],
    genes: list[str],
    conditions: list[str],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    for group in sorted(overall_stats):
        for gene in genes:
            stats = overall_stats[group][gene]
            rows.append(
                {
                    "group": group,
                    "condition": "all",
                    "gene": gene,
                    "cell_count": int(stats["cell_count"]),
                    "mean_expression": stats["mean_expression"],
                    "expression_fraction": stats["expression_fraction"],
                }
            )

    for condition in conditions:
        for group in sorted(condition_stats[condition]):
            for gene in genes:
                stats = condition_stats[condition][group][gene]
                rows.append(
                    {
                        "group": group,
                        "condition": condition,
                        "gene": gene,
                        "cell_count": int(stats["cell_count"]),
                        "mean_expression": stats["mean_expression"],
                        "expression_fraction": stats["expression_fraction"],
                    }
                )
    return rows


def safe_group_gene_stats(
    stats: dict[str, dict[str, dict[str, float]]],
    group: str,
    gene: str,
) -> dict[str, float]:
    if group not in stats:
        return {"cell_count": 0.0, "mean_expression": 0.0, "expression_fraction": 0.0}
    return stats[group].get(gene, {"cell_count": 0.0, "mean_expression": 0.0, "expression_fraction": 0.0})


def specificity_component(
    group_stats: dict[str, dict[str, dict[str, float]]],
    gene: str,
    target_group: str,
) -> float:
    target_mean = safe_group_gene_stats(group_stats, target_group, gene)["mean_expression"]
    other_means = [
        safe_group_gene_stats(group_stats, group_name, gene)["mean_expression"]
        for group_name in sorted(group_stats)
        if group_name != target_group
    ]
    baseline = sum(other_means) / float(len(other_means)) if other_means else 0.0
    return (target_mean + 0.5) / (baseline + 0.5)


def row_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (str(row["pair"]), str(row["source_group"]), str(row["target_group"]))


def compute_candidate_rows(
    payload: dict[str, Any],
    pair_specs: list[dict[str, Any]],
    groups: list[str],
    group_stats: dict[str, dict[str, dict[str, float]]],
) -> list[dict[str, Any]]:
    threshold = normalize_float(payload["expression_threshold"])
    min_cells = int(payload["min_cells"])
    rows: list[dict[str, Any]] = []

    for pair in pair_specs:
        ligand = str(pair["ligand"])
        receptor = str(pair["receptor"])
        for source_group in groups:
            ligand_stats = safe_group_gene_stats(group_stats, source_group, ligand)
            for target_group in groups:
                receptor_stats = safe_group_gene_stats(group_stats, target_group, receptor)
                expression_pass = (
                    int(ligand_stats["cell_count"]) >= min_cells
                    and int(receptor_stats["cell_count"]) >= min_cells
                    and ligand_stats["expression_fraction"] >= threshold
                    and receptor_stats["expression_fraction"] >= threshold
                )
                min_expr = (
                    min(ligand_stats["mean_expression"], receptor_stats["mean_expression"])
                    if expression_pass
                    else 0.0
                )
                geometric_mean = (
                    math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"])
                    if expression_pass
                    else 0.0
                )
                specificity = (
                    math.sqrt(
                        specificity_component(group_stats, ligand, source_group)
                        * specificity_component(group_stats, receptor, target_group)
                    )
                    if expression_pass
                    else 0.0
                )
                rows.append(
                    {
                        "pair_id": pair["pair_id"],
                        "pair": pair["pair"],
                        "ligand": ligand,
                        "receptor": receptor,
                        "evidence": pair["evidence"],
                        "source_group": source_group,
                        "target_group": target_group,
                        "resources": ";".join(pair["canonical_resources"]),
                        "resource_count": pair["resource_count"],
                        "expression_pass": expression_pass,
                        "ligand_mean": ligand_stats["mean_expression"],
                        "receptor_mean": receptor_stats["mean_expression"],
                        "ligand_fraction": ligand_stats["expression_fraction"],
                        "receptor_fraction": receptor_stats["expression_fraction"],
                        "min_expr": min_expr,
                        "geometric_mean": geometric_mean,
                        "specificity": specificity,
                    }
                )
    return rows


def compute_null_support(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    genes: list[str],
    rows: list[dict[str, Any]],
) -> None:
    threshold = normalize_float(payload["expression_threshold"])
    min_cells = int(payload["min_cells"])
    iterations = int(payload["null_permutations"])
    rng = random.Random(int(payload["random_seed"]))
    base_groups = [str(cell["normalized_group"]) for cell in cells]
    row_lookup = {row_key(row): row for row in rows}
    null_summary = {
        key: {"ge_or_equal": 0, "score_sum": 0.0, "score_max": 0.0}
        for key, row in row_lookup.items()
        if row["expression_pass"]
    }

    for _ in range(iterations):
        shuffled_groups = list(base_groups)
        rng.shuffle(shuffled_groups)
        perm_assignment = {
            str(cell["cell_id"]): shuffled_groups[index]
            for index, cell in enumerate(cells)
        }
        perm_stats = compute_group_stats(cells, perm_assignment, genes)
        for key, row in row_lookup.items():
            if not row["expression_pass"]:
                continue
            ligand_stats = safe_group_gene_stats(perm_stats, str(row["source_group"]), str(row["ligand"]))
            receptor_stats = safe_group_gene_stats(perm_stats, str(row["target_group"]), str(row["receptor"]))
            if (
                int(ligand_stats["cell_count"]) >= min_cells
                and int(receptor_stats["cell_count"]) >= min_cells
                and ligand_stats["expression_fraction"] >= threshold
                and receptor_stats["expression_fraction"] >= threshold
            ):
                perm_score = math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"])
            else:
                perm_score = 0.0
            null_summary[key]["score_sum"] += perm_score
            null_summary[key]["score_max"] = max(null_summary[key]["score_max"], perm_score)
            if perm_score >= normalize_float(row["geometric_mean"]) - 1e-12:
                null_summary[key]["ge_or_equal"] += 1

    for row in rows:
        key = row_key(row)
        if not row["expression_pass"]:
            row["null_mean_score"] = 0.0
            row["null_max_score"] = 0.0
            row["empirical_pvalue"] = 1.0
            continue
        summary = null_summary[key]
        row["null_mean_score"] = summary["score_sum"] / float(iterations)
        row["null_max_score"] = summary["score_max"]
        row["empirical_pvalue"] = (summary["ge_or_equal"] + 1.0) / float(iterations + 1)


def assign_method_ranks(rows: list[dict[str, Any]], selected_methods: list[str]) -> None:
    eligible_rows = [row for row in rows if row["expression_pass"]]
    max_rank = len(eligible_rows) + 1

    for method in selected_methods:
        reverse = METHOD_SORT_DIRECTION[method] == "desc"
        field_name = METHOD_FIELD_NAME[method]
        rank_field = METHOD_RANK_FIELD_NAME[method]
        sorted_rows = sorted(
            eligible_rows,
            key=lambda row: (
                (-normalize_float(row[field_name]) if reverse else normalize_float(row[field_name])),
                row_key(row),
            ),
        )
        for rank, row in enumerate(sorted_rows, start=1):
            row[rank_field] = float(rank)
        for row in rows:
            row.setdefault(rank_field, float(max_rank))

    magnitude_methods = [method for method in selected_methods if method in {"min_expr", "geometric_mean"}]
    specificity_methods = [method for method in selected_methods if method in {"specificity", "shuffle_pvalue"}]
    if not magnitude_methods:
        magnitude_methods = list(selected_methods)
    if not specificity_methods:
        specificity_methods = list(selected_methods)

    for row in rows:
        aggregate_rank = sum(normalize_float(row[METHOD_RANK_FIELD_NAME[method]]) for method in selected_methods) / float(len(selected_methods))
        row["magnitude_rank"] = sum(normalize_float(row[METHOD_RANK_FIELD_NAME[method]]) for method in magnitude_methods) / float(len(magnitude_methods))
        row["specificity_rank"] = sum(normalize_float(row[METHOD_RANK_FIELD_NAME[method]]) for method in specificity_methods) / float(len(specificity_methods))
        row["aggregate_rank"] = aggregate_rank
        row["aggregate_score"] = (1.0 / aggregate_rank) if row["expression_pass"] else 0.0


def build_condition_summary(
    payload: dict[str, Any],
    pair_specs: list[dict[str, Any]],
    groups: list[str],
    condition_stats: dict[str, dict[str, dict[str, dict[str, float]]]],
    conditions: list[str],
) -> list[dict[str, Any]]:
    threshold = normalize_float(payload["expression_threshold"])
    min_cells = int(payload["min_cells"])
    baseline_condition, comparison_condition = conditions
    rows: list[dict[str, Any]] = []

    for pair in pair_specs:
        ligand = str(pair["ligand"])
        receptor = str(pair["receptor"])
        for source_group in groups:
            for target_group in groups:
                baseline_ligand = safe_group_gene_stats(condition_stats[baseline_condition], source_group, ligand)
                baseline_receptor = safe_group_gene_stats(condition_stats[baseline_condition], target_group, receptor)
                comparison_ligand = safe_group_gene_stats(condition_stats[comparison_condition], source_group, ligand)
                comparison_receptor = safe_group_gene_stats(condition_stats[comparison_condition], target_group, receptor)

                def condition_score(ligand_stats: dict[str, float], receptor_stats: dict[str, float]) -> tuple[float, float]:
                    passes = (
                        int(ligand_stats["cell_count"]) >= min_cells
                        and int(receptor_stats["cell_count"]) >= min_cells
                        and ligand_stats["expression_fraction"] >= threshold
                        and receptor_stats["expression_fraction"] >= threshold
                    )
                    if not passes:
                        return 0.0, 0.0
                    return (
                        min(ligand_stats["mean_expression"], receptor_stats["mean_expression"]),
                        math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"]),
                    )

                baseline_min, baseline_geo = condition_score(baseline_ligand, baseline_receptor)
                comparison_min, comparison_geo = condition_score(comparison_ligand, comparison_receptor)
                delta_geo = comparison_geo - baseline_geo
                if abs(delta_geo) < 1e-9:
                    dominant_condition = "matched"
                elif delta_geo > 0:
                    dominant_condition = comparison_condition
                else:
                    dominant_condition = baseline_condition
                rows.append(
                    {
                        "pair": pair["pair"],
                        "source_group": source_group,
                        "target_group": target_group,
                        "baseline_condition": baseline_condition,
                        "comparison_condition": comparison_condition,
                        "baseline_min_expr": baseline_min,
                        "comparison_min_expr": comparison_min,
                        "baseline_geometric_mean": baseline_geo,
                        "comparison_geometric_mean": comparison_geo,
                        "delta_geometric_mean": delta_geo,
                        "dominant_condition": dominant_condition,
                    }
                )
    return rows


def build_adjacency_rows(
    pair_specs: list[dict[str, Any]],
    cells: list[dict[str, Any]],
    groups: list[str],
    adjacency: list[list[str]],
) -> list[dict[str, Any]]:
    adjacency_set = {tuple(sorted((str(edge[0]), str(edge[1])))) for edge in adjacency}
    group_to_cells: dict[str, list[str]] = {}
    for cell in cells:
        group_to_cells.setdefault(str(cell["normalized_group"]), []).append(str(cell["cell_id"]))

    summary_by_group_pair: dict[tuple[str, str], dict[str, Any]] = {}
    for source_group in groups:
        source_ids = group_to_cells.get(source_group, [])
        for target_group in groups:
            target_ids = group_to_cells.get(target_group, [])
            if source_group == target_group:
                possible_pairs = {
                    tuple(sorted((source_ids[i], source_ids[j])))
                    for i in range(len(source_ids))
                    for j in range(i + 1, len(source_ids))
                }
            else:
                possible_pairs = {
                    tuple(sorted((source_id, target_id)))
                    for source_id in source_ids
                    for target_id in target_ids
                    if source_id != target_id
                }
            supporting_edges = sum(1 for edge in possible_pairs if edge in adjacency_set)
            possible_edges = len(possible_pairs)
            fraction = supporting_edges / float(possible_edges) if possible_edges else 0.0
            summary_by_group_pair[(source_group, target_group)] = {
                "supporting_edges": supporting_edges,
                "possible_edges": possible_edges,
                "adjacency_fraction": fraction,
                "spatial_support": "supported" if supporting_edges > 0 else "unsupported",
            }

    rows: list[dict[str, Any]] = []
    for pair in pair_specs:
        for source_group in groups:
            for target_group in groups:
                summary = summary_by_group_pair[(source_group, target_group)]
                rows.append(
                    {
                        "pair": pair["pair"],
                        "source_group": source_group,
                        "target_group": target_group,
                        "supporting_edges": summary["supporting_edges"],
                        "possible_edges": summary["possible_edges"],
                        "adjacency_fraction": summary["adjacency_fraction"],
                        "spatial_support": summary["spatial_support"],
                    }
                )
    return rows


def lookup_by_triplet(rows: list[dict[str, Any]]) -> dict[tuple[str, str, str], dict[str, Any]]:
    return {
        (str(row["pair"]), str(row["source_group"]), str(row["target_group"])): row
        for row in rows
    }


def build_rankings(
    payload: dict[str, Any],
    rows: list[dict[str, Any]],
    condition_rows: list[dict[str, Any]],
    adjacency_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    condition_lookup = lookup_by_triplet(condition_rows)
    adjacency_lookup = lookup_by_triplet(adjacency_rows)
    selected_methods = ",".join(payload["selected_methods"])
    ranked = [row for row in rows if row["expression_pass"]]
    ranked.sort(
        key=lambda row: (
            normalize_float(row["aggregate_rank"]),
            -int(row["resource_count"]),
            row_key(row),
        )
    )

    output_rows: list[dict[str, Any]] = []
    for rank, row in enumerate(ranked[: int(payload["ranking_top_n"])], start=1):
        condition = condition_lookup[row_key(row)]
        adjacency = adjacency_lookup[row_key(row)]
        output_rows.append(
            {
                "rank": rank,
                "pair": row["pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "aggregate_score": row["aggregate_score"],
                "aggregate_rank": row["aggregate_rank"],
                "magnitude_rank": row["magnitude_rank"],
                "specificity_rank": row["specificity_rank"],
                "resource_count": row["resource_count"],
                "resources": row["resources"],
                "empirical_pvalue": row["empirical_pvalue"],
                "selected_methods": selected_methods,
                "dominant_condition": condition["dominant_condition"],
                "delta_geometric_mean": condition["delta_geometric_mean"],
                "adjacency_fraction": adjacency["adjacency_fraction"],
                "spatial_support": adjacency["spatial_support"],
            }
        )
    return output_rows


def build_method_qc_rows(
    payload: dict[str, Any],
    rows: list[dict[str, Any]],
    condition_rows: list[dict[str, Any]],
    adjacency_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    condition_lookup = lookup_by_triplet(condition_rows)
    adjacency_lookup = lookup_by_triplet(adjacency_rows)
    qc_rows: list[dict[str, Any]] = []

    for row in sorted(rows, key=row_key):
        condition = condition_lookup[row_key(row)]
        adjacency = adjacency_lookup[row_key(row)]
        qc_row = dict(row)
        qc_row["selected_methods"] = ",".join(payload["selected_methods"])
        qc_row["delta_geometric_mean"] = condition["delta_geometric_mean"]
        qc_row["dominant_condition"] = condition["dominant_condition"]
        qc_row["adjacency_fraction"] = adjacency["adjacency_fraction"]
        qc_row["spatial_support"] = adjacency["spatial_support"]
        qc_rows.append(qc_row)
    return qc_rows


def build_report_sections(
    payload: dict[str, Any],
    rankings: list[dict[str, Any]],
    resource_rows: list[dict[str, Any]],
    summary: dict[str, Any],
) -> list[dict[str, Any]]:
    top_interaction = summary["top_interaction"]
    top_condition_delta = summary["top_condition_delta"]
    sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Selected methods: {', '.join(payload['selected_methods'])}.",
                f"Expression threshold: {payload['expression_threshold']}; null permutations: {payload['null_permutations']}.",
                f"Harmonized groups: {', '.join(summary['counts']['groups'])}; conditions: {', '.join(summary['counts']['conditions'])}.",
            ],
        },
        {
            "name": "Rankings",
            "bullets": [
                f"Top interaction: {top_interaction['pair']} from {top_interaction['source_group']} to {top_interaction['target_group']} with aggregate score {top_interaction['aggregate_score']:.3f}.",
                f"Top interaction empirical p-value: {top_interaction['empirical_pvalue']:.3f}; resource count: {top_interaction['resource_count']}.",
                f"Top condition-aware gain: {top_condition_delta['pair']} with delta geometric mean {top_condition_delta['delta_geometric_mean']:.3f}.",
            ],
        },
        {
            "name": "Resource overlap",
            "bullets": [
                f"Highest overlap pair: {resource_rows[0]['pair']} supported by {resource_rows[0]['resource_count']} canonical resources.",
                f"Canonical resources seen: {', '.join(summary['counts']['resources'])}.",
            ],
        },
        {
            "name": "Condition-aware summary",
            "bullets": [
                f"Comparison order: {payload['harmonization']['condition_order'][0]} to {payload['harmonization']['condition_order'][1]}.",
                f"Top ranked interaction dominant condition: {rankings[0]['dominant_condition']}.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter uses a deterministic mean-rank approximation instead of the full LIANA or LIANA+ implementation.",
                "The shuffle null is intentionally tiny and the catalog plus adjacency graph are synthetic.",
                "Use the full public stack for real resource snapshots, curated complexes, and biological interpretation.",
            ],
        },
    ]
    return sections


def validate_tsv(path: Path, required_columns: list[str]) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AssertionError(f"Missing TSV header in {path}")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise AssertionError(f"Missing TSV columns in {path}: {missing}")


def validate_markdown(path: Path, required_sections: list[str]) -> None:
    text = path.read_text(encoding="utf-8")
    for section in required_sections:
        heading = f"## {section}"
        if heading not in text:
            raise AssertionError(f"Missing markdown section {heading} in {path}")


def validate_json(path: Path, required_keys: list[str]) -> None:
    payload = load_json(path)
    missing = [key for key in required_keys if key not in payload]
    if missing:
        raise AssertionError(f"Missing JSON keys in {path}: {missing}")


def validate_outputs(skill_dir: Path, outdir: Path) -> None:
    metadata = load_json(skill_dir / "metadata.yaml")
    for artifact in metadata["deliverables"] + metadata.get("qc_artifacts", []):
        path = outdir / artifact["path"]
        if not path.exists():
            raise AssertionError(f"Missing artifact: {path}")
        kind = artifact["kind"]
        if kind == "tsv":
            validate_tsv(path, artifact.get("required_columns", []))
        elif kind == "md":
            validate_markdown(path, artifact.get("required_sections", []))
        elif kind == "json":
            validate_json(path, artifact.get("required_keys", []))
        else:
            raise AssertionError(f"Unsupported artifact kind in metadata: {kind}")


def run_analysis(skill_dir: Path, *, input_path: Path, outdir: Path) -> dict[str, Any]:
    payload = load_json(input_path)
    validate_input(payload)
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cells, harmonization_rows, groups, conditions = harmonize_cells(payload)
    pair_specs, resource_rows = harmonize_catalog(payload)
    genes = collect_genes(cells)
    overall_stats = compute_group_stats(cells, assignment_from_cells(cells), genes)
    condition_stats = compute_condition_stats(cells, genes, conditions)
    group_rows = build_group_rows(overall_stats, condition_stats, genes, conditions)
    candidate_rows = compute_candidate_rows(payload, pair_specs, groups, overall_stats)
    compute_null_support(payload, cells, genes, candidate_rows)
    assign_method_ranks(candidate_rows, payload["selected_methods"])
    condition_rows = build_condition_summary(payload, pair_specs, groups, condition_stats, conditions)
    adjacency_rows = build_adjacency_rows(pair_specs, cells, groups, payload.get("adjacency", []))
    rankings = build_rankings(payload, candidate_rows, condition_rows, adjacency_rows)
    method_qc_rows = build_method_qc_rows(payload, candidate_rows, condition_rows, adjacency_rows)

    if not rankings:
        raise AssertionError("The starter produced no ranked interactions after filtering.")

    counts = {
        "cells": len(cells),
        "pairs": len(pair_specs),
        "groups": groups,
        "conditions": conditions,
        "resources": sorted({resource for row in resource_rows for resource in str(row["resources"]).split(";") if resource}),
        "candidate_rows": len(candidate_rows),
        "expression_pass_rows": sum(1 for row in candidate_rows if row["expression_pass"]),
    }
    candidate_lookup = {row_key(row): row for row in candidate_rows}
    eligible_condition_rows = [
        row
        for row in condition_rows
        if candidate_lookup.get(row_key(row), {}).get("expression_pass")
    ]
    top_condition_delta = max(
        eligible_condition_rows or condition_rows,
        key=lambda row: (normalize_float(row["delta_geometric_mean"]), row_key(row)),
    )
    summary = {
        "run_label": payload["run_label"],
        "selected_methods": payload["selected_methods"],
        "counts": counts,
        "top_interaction": {
            "pair": rankings[0]["pair"],
            "source_group": rankings[0]["source_group"],
            "target_group": rankings[0]["target_group"],
            "aggregate_score": rankings[0]["aggregate_score"],
            "resource_count": rankings[0]["resource_count"],
            "empirical_pvalue": rankings[0]["empirical_pvalue"],
        },
        "top_condition_delta": {
            "pair": top_condition_delta["pair"],
            "source_group": top_condition_delta["source_group"],
            "target_group": top_condition_delta["target_group"],
            "delta_geometric_mean": top_condition_delta["delta_geometric_mean"],
            "dominant_condition": top_condition_delta["dominant_condition"],
        },
    }
    report_sections = build_report_sections(payload, rankings, resource_rows, summary)

    write_tsv(
        outdir / "liana_rankings.tsv",
        rankings,
        [
            "rank",
            "pair",
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "aggregate_score",
            "aggregate_rank",
            "magnitude_rank",
            "specificity_rank",
            "resource_count",
            "resources",
            "empirical_pvalue",
            "selected_methods",
            "dominant_condition",
            "delta_geometric_mean",
            "adjacency_fraction",
            "spatial_support",
        ],
    )
    write_tsv(outdir / "resource_overlap.tsv", resource_rows, ["pair", "resource_count", "resources"])
    write_markdown(outdir / "liana_report.md", REPORT_TITLE, report_sections)
    write_tsv(
        outdir / "qc_metadata_harmonization.tsv",
        harmonization_rows,
        ["cell_id", "raw_group", "normalized_group", "raw_condition", "normalized_condition", "harmonized"],
    )
    write_tsv(
        outdir / "qc_group_expression.tsv",
        group_rows,
        ["group", "condition", "gene", "cell_count", "mean_expression", "expression_fraction"],
    )
    write_tsv(
        outdir / "qc_method_scores.tsv",
        method_qc_rows,
        [
            "pair",
            "source_group",
            "target_group",
            "expression_pass",
            "ligand_mean",
            "receptor_mean",
            "ligand_fraction",
            "receptor_fraction",
            "min_expr",
            "geometric_mean",
            "specificity",
            "null_mean_score",
            "null_max_score",
            "empirical_pvalue",
            "min_expr_rank",
            "geometric_mean_rank",
            "specificity_method_rank",
            "specificity_rank",
            "shuffle_pvalue_rank",
            "aggregate_rank",
            "aggregate_score",
            "delta_geometric_mean",
            "dominant_condition",
            "adjacency_fraction",
            "spatial_support",
            "resources",
            "resource_count",
            "selected_methods",
        ],
    )
    write_tsv(
        outdir / "qc_condition_summary.tsv",
        condition_rows,
        [
            "pair",
            "source_group",
            "target_group",
            "baseline_condition",
            "comparison_condition",
            "baseline_min_expr",
            "comparison_min_expr",
            "baseline_geometric_mean",
            "comparison_geometric_mean",
            "delta_geometric_mean",
            "dominant_condition",
        ],
    )
    write_tsv(
        outdir / "qc_adjacency_support.tsv",
        adjacency_rows,
        ["pair", "source_group", "target_group", "supporting_edges", "possible_edges", "adjacency_fraction", "spatial_support"],
    )

    summary["written_files"] = sorted(str(path.relative_to(outdir)) for path in outdir.rglob("*") if path.is_file())
    summary["written_files"].append("run_summary.json")
    summary["written_files"] = sorted(set(summary["written_files"]))
    write_json(outdir / "run_summary.json", summary)
    validate_outputs(skill_dir, outdir)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"Run the LIANA portable starter for {SKILL_DIR.name}.")
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run_analysis(SKILL_DIR, input_path=args.input, outdir=args.outdir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
