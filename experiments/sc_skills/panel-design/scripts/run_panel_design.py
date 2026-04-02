#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = SKILL_DIR / "examples" / "toy_input.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata() -> dict:
    return load_json(SKILL_DIR / "metadata.yaml")


def clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def round_float(value: float) -> float:
    return round(value, 4)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def write_markdown(path: Path, title: str, sections: list[dict]) -> None:
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_markdown_sections(path: Path) -> list[str]:
    sections: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            sections.append(line[3:].strip())
    return sections


def cosine_similarity(left: dict[str, float], right: dict[str, float], cell_types: list[str]) -> float:
    numerator = sum(left[cell_type] * right[cell_type] for cell_type in cell_types)
    left_norm = math.sqrt(sum(left[cell_type] ** 2 for cell_type in cell_types))
    right_norm = math.sqrt(sum(right[cell_type] ** 2 for cell_type in cell_types))
    if left_norm == 0.0 or right_norm == 0.0:
        return 0.0
    return numerator / (left_norm * right_norm)


def validate_payload(payload: dict) -> None:
    if "deliverables" in payload:
        raise ValueError("Toy input must contain raw inputs and expected invariants, not precomputed deliverables.")

    required_keys = {
        "candidate_markers",
        "cell_types",
        "expected_invariants",
        "panel_name",
        "panel_size",
        "platform_context",
        "run_label",
        "selection_policy",
        "source_weights",
    }
    missing = sorted(required_keys - set(payload))
    if missing:
        raise ValueError(f"Toy input is missing required keys: {', '.join(missing)}")

    if payload["panel_size"] < 1:
        raise ValueError("Panel size must be positive.")

    cell_type_names = [item["name"] for item in payload["cell_types"]]
    if len(cell_type_names) != len(set(cell_type_names)):
        raise ValueError("Cell type names must be unique.")

    minimum_quota_total = sum(int(item["quota"]) for item in payload["cell_types"])
    if payload["panel_size"] < minimum_quota_total:
        raise ValueError("Panel size is smaller than the required cell-type quotas.")

    if not payload["candidate_markers"]:
        raise ValueError("At least one candidate marker is required.")

    seen_genes: set[str] = set()
    for candidate in payload["candidate_markers"]:
        gene = candidate["gene"]
        if gene in seen_genes:
            raise ValueError(f"Candidate genes must be unique; repeated gene: {gene}")
        seen_genes.add(gene)
        signals = candidate["cell_type_signals"]
        missing_cell_types = sorted(set(cell_type_names) - set(signals))
        extra_cell_types = sorted(set(signals) - set(cell_type_names))
        if missing_cell_types or extra_cell_types:
            raise ValueError(
                f"Candidate {gene} has mismatched cell type signals; "
                f"missing={missing_cell_types}, extra={extra_cell_types}"
            )


def aggregate_markers(payload: dict) -> list[dict]:
    policy = payload["selection_policy"]
    score_weights = policy["score_weights"]
    source_weights = payload["source_weights"]
    coverage_floor = float(policy["coverage_signal_floor"])
    total_source_weight = sum(source_weights.values())
    cell_type_order = [item["name"] for item in payload["cell_types"]]
    rows: list[dict] = []

    for candidate in payload["candidate_markers"]:
        signals = candidate["cell_type_signals"]
        ordered_signals = sorted(signals.items(), key=lambda item: (-item[1], item[0]))
        best_cell_type, best_signal = ordered_signals[0]
        second_signal = ordered_signals[1][1] if len(ordered_signals) > 1 else 0.0
        evidence_score = sum(
            candidate["evidence"][source_name] * source_weight
            for source_name, source_weight in source_weights.items()
        ) / total_source_weight
        specificity_gap = max(0.0, best_signal - second_signal)
        coverage_breadth_count = sum(1 for value in signals.values() if value >= coverage_floor)
        coverage_breadth = coverage_breadth_count / len(cell_type_order)
        feasibility_score = clamp(
            candidate["probe_design_score"] - (0.5 * candidate["crowding_penalty"])
        )
        aggregate_score = (
            score_weights["evidence"] * evidence_score
            + score_weights["specificity"] * specificity_gap
            + score_weights["feasibility"] * feasibility_score
            + score_weights["breadth"] * coverage_breadth
        )
        evidence_sources = sorted(
            [
                source_name
                for source_name, source_value in candidate["evidence"].items()
                if source_value >= 0.8
            ]
        )
        row = {
            "gene": candidate["gene"],
            "best_cell_type": best_cell_type,
            "best_signal": round_float(best_signal),
            "second_signal": round_float(second_signal),
            "specificity_gap": round_float(specificity_gap),
            "coverage_breadth": round_float(coverage_breadth),
            "coverage_breadth_count": coverage_breadth_count,
            "evidence_score": round_float(evidence_score),
            "evidence_source_count": len(evidence_sources),
            "evidence_sources": ";".join(evidence_sources),
            "feasibility_score": round_float(feasibility_score),
            "aggregate_score": round_float(aggregate_score),
            "probe_design_score": round_float(candidate["probe_design_score"]),
            "crowding_penalty": round_float(candidate["crowding_penalty"]),
            "redundancy_group": candidate.get("redundancy_group", ""),
            "signals": {cell_type: round_float(signals[cell_type]) for cell_type in cell_type_order},
        }
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -row["aggregate_score"],
            -row["best_signal"],
            -row["evidence_score"],
            row["gene"],
        )
    )
    for rank, row in enumerate(rows, start=1):
        row["aggregate_rank"] = rank
    return rows


def build_redundancy_qc(aggregated_rows: list[dict], payload: dict) -> tuple[list[dict], dict[str, dict]]:
    cell_types = [item["name"] for item in payload["cell_types"]]
    threshold = float(payload["selection_policy"]["redundancy_threshold"])
    adjacency: dict[str, set[str]] = {row["gene"]: set() for row in aggregated_rows}
    pair_triggers: dict[tuple[str, str], list[str]] = {}
    max_similarity: dict[str, float] = {row["gene"]: 0.0 for row in aggregated_rows}
    gene_lookup = {row["gene"]: row for row in aggregated_rows}

    for index, left_row in enumerate(aggregated_rows):
        for right_row in aggregated_rows[index + 1 :]:
            similarity = cosine_similarity(left_row["signals"], right_row["signals"], cell_types)
            key = tuple(sorted((left_row["gene"], right_row["gene"])))
            triggers: list[str] = []
            if (
                left_row["redundancy_group"]
                and left_row["redundancy_group"] == right_row["redundancy_group"]
            ):
                triggers.append("shared_group")
            if (
                similarity >= threshold
                and left_row["best_cell_type"] == right_row["best_cell_type"]
            ):
                triggers.append("high_similarity")
            if not triggers:
                continue

            adjacency[left_row["gene"]].add(right_row["gene"])
            adjacency[right_row["gene"]].add(left_row["gene"])
            pair_triggers[key] = triggers
            max_similarity[left_row["gene"]] = max(max_similarity[left_row["gene"]], similarity)
            max_similarity[right_row["gene"]] = max(max_similarity[right_row["gene"]], similarity)

    components: list[list[str]] = []
    visited: set[str] = set()
    for row in aggregated_rows:
        gene = row["gene"]
        if gene in visited:
            continue
        stack = [gene]
        component: list[str] = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    stack.append(neighbor)
        components.append(sorted(component))

    qc_rows: list[dict] = []
    redundancy_lookup: dict[str, dict] = {}
    components.sort(key=lambda component: component[0])
    for component_index, component in enumerate(components, start=1):
        ranked_component = sorted(
            component,
            key=lambda gene: (
                -gene_lookup[gene]["aggregate_score"],
                -gene_lookup[gene]["best_signal"],
                gene,
            ),
        )
        representative_gene = ranked_component[0]
        for gene in ranked_component:
            partners = sorted(adjacency[gene])
            trigger_values = []
            for partner in partners:
                pair_key = tuple(sorted((gene, partner)))
                trigger_values.extend(pair_triggers.get(pair_key, []))
            trigger_summary = ";".join(sorted(set(trigger_values)))
            retained = gene == representative_gene
            qc_rows.append(
                {
                    "gene": gene,
                    "best_cell_type": gene_lookup[gene]["best_cell_type"],
                    "aggregate_score": gene_lookup[gene]["aggregate_score"],
                    "redundancy_component": f"component_{component_index:02d}",
                    "representative_gene": representative_gene,
                    "retained_after_pruning": "yes" if retained else "no",
                    "paired_with": ";".join(partners),
                    "trigger": trigger_summary,
                    "max_similarity": round_float(max_similarity[gene]),
                }
            )
            redundancy_lookup[gene] = {
                "component": f"component_{component_index:02d}",
                "representative_gene": representative_gene,
                "retained": retained,
                "partners": partners,
                "trigger": trigger_summary,
            }
    qc_rows.sort(
        key=lambda row: (
            row["redundancy_component"],
            0 if row["retained_after_pruning"] == "yes" else 1,
            row["gene"],
        )
    )
    return qc_rows, redundancy_lookup


def marginal_balance_gain(
    signals: dict[str, float],
    current_coverage: dict[str, float],
    targets: dict[str, float],
    signal_floor: float,
) -> float:
    gain = 0.0
    for cell_type, target in targets.items():
        signal = signals[cell_type]
        if signal < signal_floor:
            continue
        capped_before = min(current_coverage[cell_type], target)
        capped_after = min(current_coverage[cell_type] + signal, target)
        gain += max(0.0, capped_after - capped_before)
    return gain


def update_coverage(current_coverage: dict[str, float], signals: dict[str, float]) -> None:
    for cell_type, signal in signals.items():
        current_coverage[cell_type] = round_float(current_coverage[cell_type] + signal)


def infer_fill_role(candidate: dict, bridge_signal_floor: float) -> str:
    ordered = sorted(candidate["signals"].items(), key=lambda item: (-item[1], item[0]))
    second_signal = ordered[1][1] if len(ordered) > 1 else 0.0
    if second_signal >= bridge_signal_floor:
        return "immune bridge"
    return "coverage booster"


def build_rationale(candidate: dict, redundancy_lookup: dict[str, dict]) -> str:
    evidence_summary = candidate["evidence_sources"] or "weighted evidence"
    redundancy = redundancy_lookup[candidate["gene"]]
    if candidate["selection_stage"] == "quota":
        stage_reason = f"{candidate['assigned_cell_type']} quota"
    else:
        stage_reason = f"balance fill gain {candidate['marginal_gain']:.2f}"
    if redundancy["partners"]:
        redundancy_note = f"survived {redundancy['component']}"
    else:
        redundancy_note = "singleton after pruning"
    return (
        f"{stage_reason}; aggregate {candidate['aggregate_score']:.3f}; "
        f"{evidence_summary}; {redundancy_note}"
    )


def select_panel(
    aggregated_rows: list[dict],
    redundancy_lookup: dict[str, dict],
    payload: dict,
) -> tuple[list[dict], list[dict], dict[str, float]]:
    policy = payload["selection_policy"]
    retained_candidates = [
        row for row in aggregated_rows if redundancy_lookup[row["gene"]]["retained"]
    ]
    targets = {
        item["name"]: float(item["target_coverage"])
        for item in payload["cell_types"]
    }
    role_lookup = {item["name"]: item["role_label"] for item in payload["cell_types"]}
    current_coverage = {item["name"]: 0.0 for item in payload["cell_types"]}
    selected_rows: list[dict] = []
    trace_rows: list[dict] = []
    selected_genes: set[str] = set()
    balance_signal_floor = float(policy["balance_signal_floor"])
    balance_bonus_weight = float(policy["balance_bonus_weight"])
    breadth_bonus_weight = float(policy["breadth_bonus_weight"])
    bridge_signal_floor = float(policy["bridge_signal_floor"])

    def append_selection(
        candidate: dict,
        selection_stage: str,
        assigned_cell_type: str,
        panel_role: str,
        selection_priority: float,
        marginal_gain: float,
    ) -> None:
        selected_genes.add(candidate["gene"])
        update_coverage(current_coverage, candidate["signals"])
        selected_entry = dict(candidate)
        selected_entry["selection_stage"] = selection_stage
        selected_entry["assigned_cell_type"] = assigned_cell_type
        selected_entry["panel_role"] = panel_role
        selected_entry["selection_priority"] = round_float(selection_priority)
        selected_entry["marginal_gain"] = round_float(marginal_gain)
        selected_rows.append(selected_entry)

        trace_row = {
            "step": len(selected_rows),
            "selected_gene": candidate["gene"],
            "selection_stage": selection_stage,
            "assigned_cell_type": assigned_cell_type,
            "panel_role": panel_role,
            "aggregate_score": candidate["aggregate_score"],
            "marginal_gain": round_float(marginal_gain),
            "selection_priority": round_float(selection_priority),
        }
        for cell_type in current_coverage:
            trace_row[f"{cell_type}_coverage"] = round_float(current_coverage[cell_type])
        trace_rows.append(trace_row)

    for cell_type in payload["cell_types"]:
        target_name = cell_type["name"]
        for _ in range(int(cell_type["quota"])):
            eligible = [
                row
                for row in retained_candidates
                if row["gene"] not in selected_genes and row["best_cell_type"] == target_name
            ]
            if not eligible:
                eligible = [
                    row
                    for row in retained_candidates
                    if row["gene"] not in selected_genes and row["signals"][target_name] >= balance_signal_floor
                ]
            if not eligible:
                raise ValueError(f"No retained candidate can satisfy quota for cell type '{target_name}'.")

            best_candidate = max(
                eligible,
                key=lambda row: (
                    row["aggregate_score"],
                    row["signals"][target_name],
                    row["gene"],
                ),
            )
            marginal_gain = marginal_balance_gain(
                best_candidate["signals"],
                current_coverage,
                targets,
                balance_signal_floor,
            )
            selection_priority = best_candidate["aggregate_score"] + (0.1 * marginal_gain)
            append_selection(
                best_candidate,
                selection_stage="quota",
                assigned_cell_type=target_name,
                panel_role=role_lookup[target_name],
                selection_priority=selection_priority,
                marginal_gain=marginal_gain,
            )

    while len(selected_rows) < int(payload["panel_size"]):
        eligible = [
            row for row in retained_candidates if row["gene"] not in selected_genes
        ]
        if not eligible:
            raise ValueError("Not enough retained candidates to fill the requested panel size.")
        scored_candidates: list[tuple[float, dict, float]] = []
        for candidate in eligible:
            marginal_gain = marginal_balance_gain(
                candidate["signals"],
                current_coverage,
                targets,
                balance_signal_floor,
            )
            selection_priority = (
                candidate["aggregate_score"]
                + (balance_bonus_weight * marginal_gain)
                + (breadth_bonus_weight * candidate["coverage_breadth"])
            )
            scored_candidates.append((selection_priority, candidate, marginal_gain))

        selection_priority, best_candidate, marginal_gain = max(
            scored_candidates,
            key=lambda item: (
                item[0],
                item[1]["aggregate_score"],
                item[1]["gene"],
            ),
        )
        append_selection(
            best_candidate,
            selection_stage="balance_fill",
            assigned_cell_type=best_candidate["best_cell_type"],
            panel_role=infer_fill_role(best_candidate, bridge_signal_floor),
            selection_priority=selection_priority,
            marginal_gain=marginal_gain,
        )

    for row in selected_rows:
        row["rationale"] = build_rationale(row, redundancy_lookup)

    ranked_rows = sorted(
        selected_rows,
        key=lambda row: (
            -row["aggregate_score"],
            -row["selection_priority"],
            row["gene"],
        ),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row["rank"] = rank
    return ranked_rows, trace_rows, current_coverage


def write_outputs(
    outdir: Path,
    metadata: dict,
    aggregated_rows: list[dict],
    redundancy_rows: list[dict],
    selected_rows: list[dict],
    trace_rows: list[dict],
    final_coverage: dict[str, float],
    payload: dict,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    panel_columns = [
        "rank",
        "target_gene",
        "panel_role",
        "rationale",
        "selected_stage",
        "best_cell_type",
        "aggregate_score",
        "selection_priority",
        "marginal_gain",
        "evidence_source_count",
        "redundancy_component",
    ]
    panel_rows = [
        {
            "rank": row["rank"],
            "target_gene": row["gene"],
            "panel_role": row["panel_role"],
            "rationale": row["rationale"],
            "selected_stage": row["selection_stage"],
            "best_cell_type": row["best_cell_type"],
            "aggregate_score": f"{row['aggregate_score']:.4f}",
            "selection_priority": f"{row['selection_priority']:.4f}",
            "marginal_gain": f"{row['marginal_gain']:.4f}",
            "evidence_source_count": row["evidence_source_count"],
            "redundancy_component": row["redundancy_component"],
        }
        for row in selected_rows
    ]
    write_tsv(outdir / "panel_candidates.tsv", panel_rows, panel_columns)

    aggregation_columns = [
        "aggregate_rank",
        "gene",
        "best_cell_type",
        "aggregate_score",
        "evidence_score",
        "specificity_gap",
        "coverage_breadth",
        "feasibility_score",
        "evidence_sources",
    ]
    aggregation_rows = [
        {
            "aggregate_rank": row["aggregate_rank"],
            "gene": row["gene"],
            "best_cell_type": row["best_cell_type"],
            "aggregate_score": f"{row['aggregate_score']:.4f}",
            "evidence_score": f"{row['evidence_score']:.4f}",
            "specificity_gap": f"{row['specificity_gap']:.4f}",
            "coverage_breadth": f"{row['coverage_breadth']:.4f}",
            "feasibility_score": f"{row['feasibility_score']:.4f}",
            "evidence_sources": row["evidence_sources"],
        }
        for row in aggregated_rows
    ]
    write_tsv(outdir / "aggregated_marker_qc.tsv", aggregation_rows, aggregation_columns)

    redundancy_columns = [
        "gene",
        "best_cell_type",
        "aggregate_score",
        "redundancy_component",
        "representative_gene",
        "retained_after_pruning",
        "paired_with",
        "trigger",
        "max_similarity",
    ]
    redundancy_output_rows = [
        {
            "gene": row["gene"],
            "best_cell_type": row["best_cell_type"],
            "aggregate_score": f"{row['aggregate_score']:.4f}",
            "redundancy_component": row["redundancy_component"],
            "representative_gene": row["representative_gene"],
            "retained_after_pruning": row["retained_after_pruning"],
            "paired_with": row["paired_with"],
            "trigger": row["trigger"],
            "max_similarity": f"{row['max_similarity']:.4f}",
        }
        for row in redundancy_rows
    ]
    write_tsv(outdir / "redundancy_qc.tsv", redundancy_output_rows, redundancy_columns)

    trace_columns = [
        "step",
        "selected_gene",
        "selection_stage",
        "assigned_cell_type",
        "panel_role",
        "aggregate_score",
        "marginal_gain",
        "selection_priority",
    ] + [f"{item['name']}_coverage" for item in payload["cell_types"]]
    trace_output_rows = []
    for row in trace_rows:
        output_row = {
            "step": row["step"],
            "selected_gene": row["selected_gene"],
            "selection_stage": row["selection_stage"],
            "assigned_cell_type": row["assigned_cell_type"],
            "panel_role": row["panel_role"],
            "aggregate_score": f"{row['aggregate_score']:.4f}",
            "marginal_gain": f"{row['marginal_gain']:.4f}",
            "selection_priority": f"{row['selection_priority']:.4f}",
        }
        for cell_type in payload["cell_types"]:
            output_row[f"{cell_type['name']}_coverage"] = f"{row[f'{cell_type['name']}_coverage']:.4f}"
        trace_output_rows.append(output_row)
    write_tsv(outdir / "coverage_balance_qc.tsv", trace_output_rows, trace_columns)

    selected_gene_list = [row["gene"] for row in selected_rows]
    pruned_gene_list = [
        row["gene"] for row in redundancy_rows if row["retained_after_pruning"] == "no"
    ]
    panel_rationale_sections = [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}",
                f"Panel name: {payload['panel_name']}",
                f"Starter selected {len(selected_gene_list)} genes from {len(aggregated_rows)} synthetic candidates for {payload['platform_context']['assay']}.",
            ],
        },
        {
            "name": "Coverage strategy",
            "bullets": [
                f"Quota anchors: {', '.join(row['gene'] for row in selected_rows if row['selection_stage'] == 'quota')}.",
                f"Balance fill: {', '.join(row['gene'] for row in selected_rows if row['selection_stage'] == 'balance_fill') or 'none'}.",
                f"Final coverage totals: "
                + ", ".join(
                    f"{cell_type}={final_coverage[cell_type]:.2f}"
                    for cell_type in final_coverage
                ),
            ],
        },
        {
            "name": "Trade-offs",
            "bullets": [
                f"Redundancy pruning removed {len(pruned_gene_list)} genes: {', '.join(pruned_gene_list)}.",
                "The starter rewards high specificity for quota anchors, then uses a balance bonus to reinforce still-shallow compartments.",
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "Evidence values, probe-design scores, and crowding penalties are synthetic surrogates rather than direct database or sequence outputs.",
                "This starter does not implement sequence-level probe design, classifier retraining, or real atlas retrieval.",
            ],
        },
    ]
    write_markdown(outdir / "panel_rationale.md", "Panel Rationale", panel_rationale_sections)

    platform_notes_sections = [
        {
            "name": "Platform assumptions",
            "bullets": [
                f"Species: {payload['platform_context']['species']}.",
                f"Reference tissue: {payload['platform_context']['reference_tissue']}.",
                f"Panel size request: {payload['panel_size']} genes inside a nominal {payload['platform_context']['max_recommended_panel_size']}-target platform budget.",
            ],
        },
        {
            "name": "Probe constraints",
            "bullets": payload["platform_context"]["approximation_flags"],
        },
        {
            "name": "Follow-up",
            "bullets": payload["platform_context"]["design_notes"]
            + [
                "Replace the toy scores with atlas-derived marker prevalence, platform-specific probe filters, and tissue-matched references before real use.",
            ],
        },
    ]
    write_markdown(outdir / "platform_notes.md", "Platform Notes", platform_notes_sections)

    summary = {
        "aggregated_candidate_count": len(aggregated_rows),
        "coverage_by_cell_type": {cell_type: round_float(value) for cell_type, value in final_coverage.items()},
        "panel_name": payload["panel_name"],
        "pruned_gene_count": len(pruned_gene_list),
        "pruned_genes": pruned_gene_list,
        "retained_candidate_count": len(aggregated_rows) - len(pruned_gene_list),
        "run_label": payload["run_label"],
        "selected_gene_count": len(selected_gene_list),
        "selected_genes": selected_gene_list,
        "selection_trace_steps": len(trace_rows),
        "top_ranked_gene": panel_rows[0]["target_gene"] if panel_rows else "",
    }
    write_json(outdir / "run_summary.json", summary)
    validate_outputs(SKILL_DIR, outdir, metadata)


def validate_outputs(skill_dir: Path, outdir: Path, metadata: dict | None = None) -> None:
    metadata = metadata or load_json(skill_dir / "metadata.yaml")
    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise ValueError(f"Missing deliverable: {deliverable['path']}")
        if deliverable["kind"] == "tsv":
            rows = read_tsv(path)
            if not rows:
                raise ValueError(f"Deliverable is empty: {deliverable['path']}")
            required_columns = set(deliverable["required_columns"])
            observed_columns = set(rows[0])
            missing_columns = sorted(required_columns - observed_columns)
            if missing_columns:
                raise ValueError(
                    f"Deliverable {deliverable['path']} is missing columns: {', '.join(missing_columns)}"
                )
        if deliverable["kind"] == "md":
            sections = set(parse_markdown_sections(path))
            missing_sections = sorted(set(deliverable["required_sections"]) - sections)
            if missing_sections:
                raise ValueError(
                    f"Deliverable {deliverable['path']} is missing sections: {', '.join(missing_sections)}"
                )

    for extra_output in [
        "aggregated_marker_qc.tsv",
        "coverage_balance_qc.tsv",
        "redundancy_qc.tsv",
        "run_summary.json",
    ]:
        if not (outdir / extra_output).exists():
            raise ValueError(f"Missing QC artifact: {extra_output}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the portable panel-design starter.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    try:
        payload = load_json(args.input)
        validate_payload(payload)
        metadata = load_metadata()
        aggregated_rows = aggregate_markers(payload)
        redundancy_rows, redundancy_lookup = build_redundancy_qc(aggregated_rows, payload)
        for row in aggregated_rows:
            row["redundancy_component"] = redundancy_lookup[row["gene"]]["component"]
        selected_rows, trace_rows, final_coverage = select_panel(
            aggregated_rows,
            redundancy_lookup,
            payload,
        )
        write_outputs(
            args.outdir,
            metadata,
            aggregated_rows,
            redundancy_rows,
            selected_rows,
            trace_rows,
            final_coverage,
            payload,
        )
    except Exception as exc:  # pragma: no cover - exercised through CLI tests
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
