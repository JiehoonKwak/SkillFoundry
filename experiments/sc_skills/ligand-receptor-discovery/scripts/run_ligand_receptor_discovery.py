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
REPORT_TITLE = "Ligand-Receptor Discovery Starter Report"


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
        return str(value).lower()
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


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_tsv(path: Path, required_columns: list[str]) -> None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
    missing = [column for column in required_columns if column not in fieldnames]
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


def collect_genes(cells: list[dict[str, Any]]) -> list[str]:
    genes: set[str] = set()
    for cell in cells:
        genes.update(cell["expression"].keys())
    return sorted(genes)


def cell_expression(cell: dict[str, Any], gene: str) -> float:
    return normalize_float(cell["expression"].get(gene, 0.0))


def assignment_from_cells(cells: list[dict[str, Any]]) -> dict[str, str]:
    return {str(cell["cell_id"]): str(cell["group"]) for cell in cells}


def validate_input(payload: dict[str, Any]) -> None:
    require_keys(
        payload,
        [
            "run_label",
            "expression_threshold",
            "active_expression_threshold",
            "null_permutations",
            "random_seed",
            "ranking_top_n",
            "spatial_support_threshold",
            "resource_aliases",
            "resource_weights",
            "cells",
            "ligand_receptor_catalog",
        ],
        name="toy input",
    )
    cells = payload["cells"]
    if len(cells) < 4:
        raise AssertionError("The starter expects at least four cells.")
    cell_ids = [str(cell["cell_id"]) for cell in cells]
    if len(cell_ids) != len(set(cell_ids)):
        raise AssertionError("Cell IDs must be unique.")
    groups = {str(cell["group"]) for cell in cells}
    if len(groups) < 2:
        raise AssertionError("The starter expects at least two groups.")
    for cell in cells:
        require_keys(cell, ["cell_id", "group", "expression"], name=f"cell {cell.get('cell_id', '<unknown>')}")
        if not isinstance(cell["expression"], dict) or not cell["expression"]:
            raise AssertionError(f"Cell {cell['cell_id']} must have a non-empty expression mapping.")
    genes = collect_genes(cells)
    for pair in payload["ligand_receptor_catalog"]:
        require_keys(pair, ["pair_id", "ligand", "receptor", "resources", "evidence"], name="ligand_receptor_catalog[]")
        if pair["ligand"] not in genes:
            raise AssertionError(f"Missing ligand gene in toy input: {pair['ligand']}")
        if pair["receptor"] not in genes:
            raise AssertionError(f"Missing receptor gene in toy input: {pair['receptor']}")
        if not pair["resources"]:
            raise AssertionError(f"Catalog entry {pair['pair_id']} must declare at least one resource.")
    known_cells = set(cell_ids)
    for edge in payload.get("adjacency", []):
        if len(edge) != 2:
            raise AssertionError(f"Adjacency edges must contain exactly two cell IDs: {edge}")
        if edge[0] not in known_cells or edge[1] not in known_cells:
            raise AssertionError(f"Adjacency references unknown cells: {edge}")


def normalize_resource_name(name: str, aliases: dict[str, str]) -> str:
    return aliases.get(name.lower(), name)


def harmonize_catalog(payload: dict[str, Any]) -> list[dict[str, Any]]:
    aliases = {str(key).lower(): str(value) for key, value in payload["resource_aliases"].items()}
    resource_weights = {str(key): normalize_float(value) for key, value in payload["resource_weights"].items()}
    pair_specs: list[dict[str, Any]] = []
    max_weight_sum = 0.0
    for pair in payload["ligand_receptor_catalog"]:
        resources = sorted({normalize_resource_name(str(resource), aliases) for resource in pair["resources"]})
        weight_sum = sum(resource_weights.get(resource, 0.5) for resource in resources)
        max_weight_sum = max(max_weight_sum, weight_sum)
        pair_specs.append(
            {
                "pair_id": str(pair["pair_id"]),
                "candidate_pair": f"{pair['ligand']}_{pair['receptor']}",
                "ligand": str(pair["ligand"]),
                "receptor": str(pair["receptor"]),
                "resources": resources,
                "resource_count": len(resources),
                "resource_weight_sum": weight_sum,
                "pathway": str(pair.get("pathway", "unspecified")),
                "evidence": str(pair["evidence"]),
            }
        )
    max_weight_sum = max(max_weight_sum, 1.0)
    for spec in pair_specs:
        spec["knowledge_score"] = spec["resource_weight_sum"] / max_weight_sum
    return pair_specs


def compute_group_stats(
    cells: list[dict[str, Any]],
    group_assignment: dict[str, str],
    genes: list[str],
) -> dict[str, dict[str, dict[str, float]]]:
    grouped_cells: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        grouped_cells.setdefault(group_assignment[str(cell["cell_id"])], []).append(cell)

    stats: dict[str, dict[str, dict[str, float]]] = {}
    for group, group_cells in sorted(grouped_cells.items()):
        stats[group] = {}
        cell_count = float(len(group_cells))
        for gene in genes:
            values = [cell_expression(cell, gene) for cell in group_cells]
            stats[group][gene] = {
                "cell_count": cell_count,
                "mean_expression": sum(values) / cell_count,
                "expression_fraction": sum(value > 0 for value in values) / cell_count,
            }
    return stats


def build_group_expression_rows(
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


def pair_row_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (str(row["candidate_pair"]), str(row["source_group"]), str(row["target_group"]))


def compute_pair_scores(
    payload: dict[str, Any],
    pair_specs: list[dict[str, Any]],
    group_stats: dict[str, dict[str, dict[str, float]]],
) -> list[dict[str, Any]]:
    threshold = normalize_float(payload["expression_threshold"])
    groups = sorted(group_stats)
    rows: list[dict[str, Any]] = []
    for spec in pair_specs:
        for source_group in groups:
            ligand_stats = group_stats[source_group][spec["ligand"]]
            for target_group in groups:
                receptor_stats = group_stats[target_group][spec["receptor"]]
                expression_pass = (
                    ligand_stats["expression_fraction"] >= threshold
                    and receptor_stats["expression_fraction"] >= threshold
                )
                min_expr = min(ligand_stats["mean_expression"], receptor_stats["mean_expression"]) if expression_pass else 0.0
                geometric_mean = (
                    math.sqrt(ligand_stats["mean_expression"] * receptor_stats["mean_expression"]) if expression_pass else 0.0
                )
                expression_score = 0.5 * (min_expr + geometric_mean)
                rows.append(
                    {
                        "candidate_pair": spec["candidate_pair"],
                        "pair_id": spec["pair_id"],
                        "source_group": source_group,
                        "target_group": target_group,
                        "ligand": spec["ligand"],
                        "receptor": spec["receptor"],
                        "pathway": spec["pathway"],
                        "evidence": spec["evidence"],
                        "resources": ";".join(spec["resources"]),
                        "resource_count": spec["resource_count"],
                        "knowledge_score": spec["knowledge_score"],
                        "expression_pass": expression_pass,
                        "ligand_mean": ligand_stats["mean_expression"],
                        "receptor_mean": receptor_stats["mean_expression"],
                        "ligand_fraction": ligand_stats["expression_fraction"],
                        "receptor_fraction": receptor_stats["expression_fraction"],
                        "min_expr": min_expr,
                        "geometric_mean": geometric_mean,
                        "expression_score": expression_score,
                    }
                )
    return rows


def generate_permuted_assignments(cells: list[dict[str, Any]], permutations: int, seed: int) -> list[dict[str, str]]:
    labels = [str(cell["group"]) for cell in cells]
    baseline = tuple(labels)
    rng = random.Random(seed)
    seen: set[tuple[str, ...]] = set()
    assignments: list[dict[str, str]] = []
    attempts = 0
    max_attempts = max(100, permutations * 30)
    while len(assignments) < permutations and attempts < max_attempts:
        attempts += 1
        shuffled = labels[:]
        rng.shuffle(shuffled)
        candidate = tuple(shuffled)
        if candidate == baseline or candidate in seen:
            continue
        seen.add(candidate)
        assignments.append({str(cells[index]["cell_id"]): shuffled[index] for index in range(len(cells))})
    if len(assignments) != permutations:
        raise AssertionError("Could not generate the requested number of unique deterministic shuffles.")
    return assignments


def summarize_null_statistics(
    payload: dict[str, Any],
    cells: list[dict[str, Any]],
    genes: list[str],
    pair_specs: list[dict[str, Any]],
    observed_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    permutations = int(payload["null_permutations"])
    null_scores: dict[tuple[str, str, str], list[float]] = {pair_row_key(row): [] for row in observed_rows}
    for assignment in generate_permuted_assignments(cells, permutations, int(payload["random_seed"])):
        permuted_stats = compute_group_stats(cells, assignment, genes)
        permuted_rows = compute_pair_scores(payload, pair_specs, permuted_stats)
        for row in permuted_rows:
            null_scores[pair_row_key(row)].append(normalize_float(row["expression_score"]))

    null_rows: list[dict[str, Any]] = []
    for row in observed_rows:
        stats = null_scores[pair_row_key(row)]
        observed_score = normalize_float(row["expression_score"])
        ge_count = sum(value >= observed_score for value in stats)
        empirical_pvalue = (ge_count + 1) / float(permutations + 1)
        row["null_mean_score"] = sum(stats) / float(len(stats))
        row["null_max_score"] = max(stats)
        row["null_ge_count"] = ge_count
        row["empirical_pvalue"] = empirical_pvalue
        null_rows.append(
            {
                "candidate_pair": row["candidate_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "observed_expression_score": observed_score,
                "null_mean_score": row["null_mean_score"],
                "null_max_score": row["null_max_score"],
                "null_ge_count": ge_count,
                "empirical_pvalue": empirical_pvalue,
            }
        )
    return null_rows


def build_spatial_support_rows(payload: dict[str, Any], cells: list[dict[str, Any]], pair_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cell_lookup = {str(cell["cell_id"]): cell for cell in cells}
    group_assignment = assignment_from_cells(cells)
    threshold = normalize_float(payload["active_expression_threshold"])
    support_threshold = normalize_float(payload["spatial_support_threshold"])
    spatial_rows: list[dict[str, Any]] = []
    for row in pair_rows:
        possible_edges = 0
        supporting_edges = 0
        for cell_a, cell_b in payload.get("adjacency", []):
            for source_cell, target_cell in ((cell_a, cell_b), (cell_b, cell_a)):
                if group_assignment[source_cell] != row["source_group"] or group_assignment[target_cell] != row["target_group"]:
                    continue
                possible_edges += 1
                if (
                    cell_expression(cell_lookup[source_cell], str(row["ligand"])) >= threshold
                    and cell_expression(cell_lookup[target_cell], str(row["receptor"])) >= threshold
                ):
                    supporting_edges += 1
        adjacency_fraction = supporting_edges / float(possible_edges) if possible_edges else 0.0
        spatial_support = "supported" if possible_edges and adjacency_fraction >= support_threshold else "not_supported"
        if possible_edges == 0:
            spatial_support = "no_edges"
        row["supporting_edges"] = supporting_edges
        row["possible_edges"] = possible_edges
        row["adjacency_fraction"] = adjacency_fraction
        row["spatial_support"] = spatial_support
        spatial_rows.append(
            {
                "candidate_pair": row["candidate_pair"],
                "source_group": row["source_group"],
                "target_group": row["target_group"],
                "ligand": row["ligand"],
                "receptor": row["receptor"],
                "supporting_edges": supporting_edges,
                "possible_edges": possible_edges,
                "adjacency_fraction": adjacency_fraction,
                "spatial_support": spatial_support,
            }
        )
    return spatial_rows


def assign_support_scores(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    passing_rows = [row for row in rows if row["expression_pass"]]
    if not passing_rows:
        raise AssertionError("No ligand-receptor candidates passed the expression filter.")
    max_expression = max(normalize_float(row["expression_score"]) for row in passing_rows) or 1.0
    for row in rows:
        if not row["expression_pass"]:
            row["support_score"] = 0.0
            continue
        expression_component = normalize_float(row["expression_score"]) / max_expression
        knowledge_component = normalize_float(row["knowledge_score"])
        null_component = 1.0 - normalize_float(row["empirical_pvalue"])
        spatial_component = normalize_float(row["adjacency_fraction"])
        row["support_score"] = (
            0.55 * expression_component
            + 0.20 * knowledge_component
            + 0.15 * null_component
            + 0.10 * spatial_component
        )
    ranked_rows = sorted(
        passing_rows,
        key=lambda row: (
            -normalize_float(row["support_score"]),
            normalize_float(row["empirical_pvalue"]),
            -normalize_float(row["knowledge_score"]),
            pair_row_key(row),
        ),
    )
    for index, row in enumerate(ranked_rows, start=1):
        row["rank"] = index
    return ranked_rows


def build_candidate_pairs(ranked_rows: list[dict[str, Any]], top_n: int) -> list[dict[str, Any]]:
    selected = ranked_rows[:top_n]
    return [
        {
            "rank": row["rank"],
            "candidate_pair": row["candidate_pair"],
            "source_group": row["source_group"],
            "target_group": row["target_group"],
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "support_score": row["support_score"],
            "expression_score": row["expression_score"],
            "knowledge_score": row["knowledge_score"],
            "resource_count": row["resource_count"],
            "empirical_pvalue": row["empirical_pvalue"],
            "adjacency_fraction": row["adjacency_fraction"],
            "spatial_support": row["spatial_support"],
            "pathway": row["pathway"],
        }
        for row in selected
    ]


def build_evidence_table(candidate_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    evidence_rows: list[dict[str, Any]] = []
    for row in candidate_rows:
        for evidence_type, evidence_value in (
            ("resources", row["resources"]),
            ("knowledge_score", format_value(row["knowledge_score"])),
            ("expression_means", f"ligand={format_value(row['ligand_mean'])};receptor={format_value(row['receptor_mean'])}"),
            ("expression_filter", f"ligand_fraction={format_value(row['ligand_fraction'])};receptor_fraction={format_value(row['receptor_fraction'])}"),
            ("geometric_mean", format_value(row["geometric_mean"])),
            ("empirical_pvalue", format_value(row["empirical_pvalue"])),
            ("spatial_support", f"{row['spatial_support']}:{format_value(row['adjacency_fraction'])}"),
            ("catalog_evidence", row["evidence"]),
        ):
            evidence_rows.append(
                {
                    "candidate_pair": row["candidate_pair"],
                    "source_group": row["source_group"],
                    "target_group": row["target_group"],
                    "evidence_type": evidence_type,
                    "evidence_value": evidence_value,
                }
            )
    return evidence_rows


def build_report_sections(payload: dict[str, Any], candidate_rows: list[dict[str, Any]], summary: dict[str, Any]) -> list[dict[str, Any]]:
    top_row = candidate_rows[0]
    prioritized_bullets = [
        (
            f"Rank {row['rank']}: {row['candidate_pair']} from {row['source_group']} to {row['target_group']} "
            f"with support_score={format_value(row['support_score'])}, empirical_pvalue={format_value(row['empirical_pvalue'])}, "
            f"and spatial_support={row['spatial_support']}."
        )
        for row in candidate_rows
    ]
    return [
        {
            "name": "Run context",
            "bullets": [
                f"Cells: {summary['counts']['cells']}; groups: {summary['counts']['groups']}; catalog pairs: {summary['counts']['catalog_pairs']}.",
                (
                    f"Expression threshold: {format_value(payload['expression_threshold'])}; active edge threshold: "
                    f"{format_value(payload['active_expression_threshold'])}; deterministic shuffles: {payload['null_permutations']}."
                ),
                (
                    f"Top candidate after ranking: {top_row['candidate_pair']} from {top_row['source_group']} to "
                    f"{top_row['target_group']}."
                ),
            ],
        },
        {
            "name": "Prioritized candidates",
            "bullets": prioritized_bullets,
        },
        {
            "name": "Evidence summary",
            "bullets": [
                (
                    f"{top_row['candidate_pair']} is backed by {top_row['resource_count']} catalog resources "
                    f"({top_row['resources']}) and keeps a knowledge_score of {format_value(top_row['knowledge_score'])}."
                ),
                (
                    f"The top pair has observed expression_score={format_value(top_row['expression_score'])} versus "
                    f"null_mean_score={format_value(top_row['null_mean_score'])} and empirical_pvalue={format_value(top_row['empirical_pvalue'])}."
                ),
                (
                    f"Spatial support is {top_row['spatial_support']} with {top_row['supporting_edges']} supporting edges out of "
                    f"{top_row['possible_edges']} observed adjacency opportunities."
                ),
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "This starter uses a tiny synthetic catalog and raw toy matrix rather than a real OmniPath or CellPhoneDB resource snapshot.",
                "The null model is intentionally small and deterministic, so the empirical p-values are only starter-scale heuristics.",
                "Adjacency support is a graph-edge fraction surrogate, not a full Squidpy or LIANA+ spatial statistic.",
            ],
        },
    ]


def run_analysis(skill_dir: Path, *, input_path: Path, outdir: Path) -> dict[str, Any]:
    payload = load_json(input_path)
    validate_input(payload)
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    cells = payload["cells"]
    genes = collect_genes(cells)
    pair_specs = harmonize_catalog(payload)
    baseline_assignment = assignment_from_cells(cells)
    group_stats = compute_group_stats(cells, baseline_assignment, genes)
    group_rows = build_group_expression_rows(group_stats, genes)
    pair_rows = compute_pair_scores(payload, pair_specs, group_stats)
    null_rows = summarize_null_statistics(payload, cells, genes, pair_specs, pair_rows)
    spatial_rows = build_spatial_support_rows(payload, cells, pair_rows)
    ranked_rows = assign_support_scores(pair_rows)
    candidate_rows = build_candidate_pairs(ranked_rows, int(payload["ranking_top_n"]))
    evidence_rows = build_evidence_table(ranked_rows[: int(payload["ranking_top_n"])])

    summary = {
        "run_label": payload["run_label"],
        "counts": {
            "cells": len(cells),
            "groups": sorted(group_stats),
            "catalog_pairs": len(pair_specs),
            "ranked_candidates": len(ranked_rows),
            "spatially_supported": sum(1 for row in ranked_rows if row["spatial_support"] == "supported"),
        },
        "top_candidate": {
            "candidate_pair": candidate_rows[0]["candidate_pair"],
            "source_group": candidate_rows[0]["source_group"],
            "target_group": candidate_rows[0]["target_group"],
            "support_score": candidate_rows[0]["support_score"],
            "empirical_pvalue": candidate_rows[0]["empirical_pvalue"],
        },
    }
    report_sections = build_report_sections(payload, ranked_rows[: int(payload["ranking_top_n"])], summary)

    write_tsv(
        outdir / "candidate_pairs.tsv",
        candidate_rows,
        [
            "rank",
            "candidate_pair",
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "support_score",
            "expression_score",
            "knowledge_score",
            "resource_count",
            "empirical_pvalue",
            "adjacency_fraction",
            "spatial_support",
            "pathway",
        ],
    )
    write_tsv(
        outdir / "evidence_table.tsv",
        evidence_rows,
        ["candidate_pair", "source_group", "target_group", "evidence_type", "evidence_value"],
    )
    write_markdown(outdir / "discovery_summary.md", REPORT_TITLE, report_sections)
    write_tsv(
        outdir / "qc_group_expression.tsv",
        group_rows,
        ["group", "gene", "cell_count", "mean_expression", "expression_fraction"],
    )
    write_tsv(
        outdir / "qc_pair_scores.tsv",
        pair_rows,
        [
            "candidate_pair",
            "source_group",
            "target_group",
            "expression_pass",
            "ligand",
            "receptor",
            "ligand_mean",
            "receptor_mean",
            "ligand_fraction",
            "receptor_fraction",
            "min_expr",
            "geometric_mean",
            "expression_score",
            "resource_count",
            "knowledge_score",
            "resources",
            "pathway",
            "evidence",
        ],
    )
    write_tsv(
        outdir / "qc_shuffle_null.tsv",
        null_rows,
        [
            "candidate_pair",
            "source_group",
            "target_group",
            "observed_expression_score",
            "null_mean_score",
            "null_max_score",
            "null_ge_count",
            "empirical_pvalue",
        ],
    )
    write_tsv(
        outdir / "qc_spatial_support.tsv",
        spatial_rows,
        [
            "candidate_pair",
            "source_group",
            "target_group",
            "ligand",
            "receptor",
            "supporting_edges",
            "possible_edges",
            "adjacency_fraction",
            "spatial_support",
        ],
    )

    summary["written_files"] = sorted(str(path.relative_to(outdir)) for path in outdir.rglob("*") if path.is_file())
    summary["written_files"].append("run_summary.json")
    summary["written_files"] = sorted(set(summary["written_files"]))
    write_json(outdir / "run_summary.json", summary)
    validate_outputs(skill_dir, outdir)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=f"Run the ligand-receptor portable starter for {SKILL_DIR.name}.")
    parser.add_argument("--input", type=Path, default=SKILL_DIR / "examples" / "toy_input.json")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = run_analysis(SKILL_DIR, input_path=args.input, outdir=args.outdir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
