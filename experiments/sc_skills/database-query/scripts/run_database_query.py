#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = SKILL_DIR / "examples" / "toy_input.json"
ENTITY_TYPES = ("gene", "protein", "variant", "disease")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata() -> dict:
    return load_json(SKILL_DIR / "metadata.yaml")


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def ordered_unique(items: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def normalize_lookup_key(text: str) -> str:
    lowered = text.strip().lower()
    lowered = lowered.replace("_", " ")
    lowered = re.sub(r"[-/]+", " ", lowered)
    lowered = re.sub(r"\s+", " ", lowered)
    lowered = re.sub(r"[^a-z0-9.:() +]+", "", lowered)
    return re.sub(r"\s+", " ", lowered).strip()


def candidate_keys(text: str) -> list[str]:
    stripped = text.strip()
    primary = normalize_lookup_key(stripped)
    variants = [primary]
    compact = primary.replace(" ", "")
    if compact and compact != primary:
        variants.append(compact)
    upper_key = normalize_lookup_key(stripped.upper())
    if upper_key and upper_key not in variants:
        variants.append(upper_key)
    bare = normalize_lookup_key(re.sub(r"\s+", "", stripped))
    if bare and bare not in variants:
        variants.append(bare)
    return [value for value in ordered_unique(variants) if value]


def is_uniprot_accession(text: str) -> bool:
    token = text.strip().upper()
    patterns = (
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]",
        r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}",
    )
    return any(re.fullmatch(pattern, token) for pattern in patterns)


def sort_sources(sources: list[str], source_weights: dict[str, float]) -> list[str]:
    return sorted(ordered_unique(sources), key=lambda item: (-source_weights.get(item, 0.0), item))


def build_concept_index(payload: dict) -> dict[str, dict]:
    return {
        concept["concept_id"]: concept
        for concept in payload["local_lookup_tables"]["concepts"]
    }


def build_alias_index(payload: dict, concepts: dict[str, dict]) -> dict[str, list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    for alias in payload["local_lookup_tables"]["aliases"]:
        concept = concepts[alias["concept_id"]]
        entry = {
            "concept_id": alias["concept_id"],
            "entity_type": alias["entity_type"],
            "alias_source": alias["alias_source"],
            "label": concept["canonical_label"],
            "normalized_id": concept["preferred_id"],
        }
        for key in candidate_keys(alias["alias"]):
            index[key].append(entry)
    return index


def build_source_indexes(payload: dict) -> dict[str, dict[str, list[dict]]]:
    indexes: dict[str, dict[str, list[dict]]] = {}
    for source_name, source_payload in payload["mock_api_payloads"].items():
        source_index: dict[str, list[dict]] = defaultdict(list)
        for record in source_payload.get("records", []):
            normalized_query_keys = ordered_unique(
                [
                    key
                    for value in record.get("query_keys", [])
                    for key in candidate_keys(value)
                ]
            )
            normalized_source_id = normalize_lookup_key(record["source_id"])
            entry = dict(record)
            entry["source"] = source_name
            entry["normalized_query_keys"] = normalized_query_keys
            entry["normalized_source_id"] = normalized_source_id
            for key in ordered_unique(normalized_query_keys + [normalized_source_id]):
                if key:
                    source_index[key].append(entry)
        indexes[source_name] = source_index
    return indexes


def collect_hits(
    query: str,
    keys: list[str],
    alias_index: dict[str, list[dict]],
    source_indexes: dict[str, dict[str, list[dict]]],
    source_weights: dict[str, float],
) -> list[dict]:
    hits: list[dict] = []
    seen: set[tuple[str, str, str]] = set()

    for key in keys:
        for alias in alias_index.get(key, []):
            dedupe_key = ("local_alias", alias["concept_id"], alias["alias_source"])
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            hits.append(
                {
                    "query": query,
                    "source": "local_alias",
                    "concept_id": alias["concept_id"],
                    "entity_type": alias["entity_type"],
                    "source_id": alias["alias_source"],
                    "matched_key": key,
                    "match_kind": "alias",
                    "base_score": float(source_weights.get("local_alias", 0.0)),
                    "label": alias["label"],
                }
            )

    for source_name, source_index in source_indexes.items():
        for key in keys:
            for record in source_index.get(key, []):
                dedupe_key = (source_name, record["concept_id"], record["source_id"])
                if dedupe_key in seen:
                    continue
                seen.add(dedupe_key)
                hits.append(
                    {
                        "query": query,
                        "source": source_name,
                        "concept_id": record["concept_id"],
                        "entity_type": record["entity_type"],
                        "source_id": record["source_id"],
                        "matched_key": key,
                        "match_kind": "source_id" if key == record["normalized_source_id"] else "query_key",
                        "base_score": float(source_weights.get(source_name, 0.0)),
                        "label": record["label"],
                        "clinical_significance": record.get("clinical_significance", ""),
                    }
                )
    return hits


def infer_entity_type(
    query: str,
    keys: list[str],
    hits: list[dict],
    payload: dict,
) -> tuple[str, dict[str, float], int]:
    scores = {entity_type: 0.0 for entity_type in ENTITY_TYPES}
    alias_hits = sum(1 for hit in hits if hit["source"] == "local_alias")
    lookup = keys[0] if keys else normalize_lookup_key(query)
    compact = lookup.replace(" ", "")
    keywords = payload["query_type_keywords"]

    if compact.startswith("rs") and compact[2:].isdigit():
        scores["variant"] += 7.0
    if any(marker in lookup for marker in keywords.get("variant_markers", []) if marker != "rs"):
        scores["variant"] += 4.0
    if is_uniprot_accession(query):
        scores["protein"] += 7.0
    if any(term in lookup for term in keywords.get("protein", [])):
        scores["protein"] += 2.0
    if any(term in lookup for term in keywords.get("disease", [])):
        scores["disease"] += 4.0
    if " " not in lookup and not is_uniprot_accession(query) and not compact.startswith("rs"):
        scores["gene"] += 1.5

    for hit in hits:
        weight = 1.25 if hit["source"] == "local_alias" else 1.0
        scores[hit["entity_type"]] += weight

    inferred = sorted(ENTITY_TYPES, key=lambda entity: (-scores[entity], entity))[0]
    return inferred, {key: round(value, 3) for key, value in scores.items()}, alias_hits


def reconcile_hits(
    query: str,
    inferred_entity_type: str,
    hits: list[dict],
    concepts: dict[str, dict],
    payload: dict,
) -> tuple[dict, list[dict]]:
    if not hits:
        raise ValueError(f"No hits were found for query: {query!r}")

    grouped: dict[str, list[dict]] = defaultdict(list)
    for hit in hits:
        grouped[hit["concept_id"]].append(hit)

    scored_groups: list[dict] = []
    preferred_source = payload["preferred_sources"][inferred_entity_type]
    source_weights = payload["source_weights"]

    for concept_id, group_hits in grouped.items():
        concept = concepts[concept_id]
        public_sources = sort_sources(
            [hit["source"] for hit in group_hits if hit["source"] != "local_alias"],
            source_weights,
        )
        base_score = sum(hit["base_score"] for hit in group_hits)
        type_bonus = 3.0 if concept["entity_type"] == inferred_entity_type else -1.5
        preferred_bonus = 1.5 if preferred_source in public_sources else 0.0
        cross_source_bonus = 0.75 * max(len(public_sources) - 1, 0)
        local_alias_bonus = 0.5 if any(hit["source"] == "local_alias" for hit in group_hits) else 0.0
        score = round(base_score + type_bonus + preferred_bonus + cross_source_bonus + local_alias_bonus, 3)
        scored_groups.append(
            {
                "concept_id": concept_id,
                "entity_type": concept["entity_type"],
                "resolved_id": concept["preferred_id"],
                "match_label": concept["canonical_label"],
                "supporting_sources": public_sources,
                "source_hit_count": len(public_sources),
                "winning_score": score,
                "all_hits": group_hits,
                "local_alias_count": sum(1 for hit in group_hits if hit["source"] == "local_alias"),
            }
        )

    scored_groups.sort(key=lambda item: (-item["winning_score"], item["concept_id"]))
    winner = dict(scored_groups[0])
    runner_up = scored_groups[1] if len(scored_groups) > 1 else None
    winner["runner_up_concept_id"] = runner_up["concept_id"] if runner_up else ""
    winner["score_margin"] = round(
        winner["winning_score"] - (runner_up["winning_score"] if runner_up else 0.0),
        3,
    )
    winner["ambiguity_flag"] = "yes" if runner_up else "no"
    winner["cross_source_support"] = "yes" if winner["source_hit_count"] >= 2 else "no"
    return winner, scored_groups


def build_query_result(query: str, winner: dict) -> dict:
    evidence_parts = [
        f"entity={winner['entity_type']}",
        f"supporting_sources={winner['source_hit_count']}",
        f"score_margin={winner['score_margin']:.3f}",
        f"local_alias_hits={winner['local_alias_count']}",
    ]
    if winner["ambiguity_flag"] == "yes":
        evidence_parts.append(f"runner_up={winner['runner_up_concept_id']}")
    clinical_significance = next(
        (
            hit.get("clinical_significance", "")
            for hit in winner["all_hits"]
            if hit.get("clinical_significance")
        ),
        "",
    )
    if clinical_significance:
        evidence_parts.append(f"clinical_significance={clinical_significance}")
    sources = winner["supporting_sources"] or ["local_alias"]
    return {
        "query": query,
        "normalized_id": winner["resolved_id"],
        "source": ";".join(sources),
        "match_label": winner["match_label"],
        "evidence": "|".join(evidence_parts),
    }


def build_source_provenance_sections(
    payload: dict,
    normalization_rows: list[dict],
    fanout_rows: list[dict],
) -> list[dict]:
    source_hit_counts: dict[str, int] = defaultdict(int)
    for row in fanout_rows:
        if row["source"] != "local_alias":
            source_hit_counts[row["source"]] += 1
    api_bullets = [
        f"{source}: {source_hit_counts.get(source, 0)} fan-out hit(s) across the bundled mock payload."
        for source in payload["mock_api_payloads"].keys()
        if source_hit_counts.get(source, 0) > 0
    ]
    return [
        {
            "name": "Run context",
            "bullets": [
                f"Run label: {payload['run_label']}",
                f"Queries processed: {len(payload['queries'])}",
                "Execution stayed local and deterministic with no network access.",
                "The starter wrote machine-readable QC tables for normalization, fan-out, and reconciliation.",
            ],
        },
        {
            "name": "APIs consulted",
            "bullets": api_bullets or ["No mock public-source hits were observed."],
        },
        {
            "name": "Normalization",
            "bullets": [
                f"{row['query']} -> {row['normalized_query']} [{row['inferred_entity_type']}]"
                for row in normalization_rows
            ],
        },
        {
            "name": "Caveats",
            "bullets": [
                "Source payloads are tiny bundled surrogates, not live API responses.",
                "Identifier syntax may resemble public accessions while still remaining toy fixtures.",
                "Full HGVS parsing, clinical evidence review, and live pagination are outside starter scope.",
            ],
        },
    ]


def build_resolution_sections(results_rows: list[dict], reconciliation_rows: list[dict]) -> list[dict]:
    resolved_bullets = [
        f"{row['query']} -> {row['normalized_id']} via {row['source']} ({row['match_label']})"
        for row in results_rows
    ]
    ambiguous_rows = [row for row in reconciliation_rows if row["ambiguity_flag"] == "yes"]
    ambiguity_bullets = [
        f"{row['query']} kept a runner-up {row['runner_up_concept_id']} but resolved to {row['winning_concept_id']} with margin {row['score_margin']}."
        for row in ambiguous_rows
    ] or ["No runner-up concepts survived scoring for this toy run."]
    return [
        {
            "name": "Resolved identifiers",
            "bullets": resolved_bullets,
        },
        {
            "name": "Ambiguity handling",
            "bullets": ambiguity_bullets,
        },
        {
            "name": "Follow-up",
            "bullets": [
                "Replace mock payloads with pinned official API snapshots before relying on the output outside the starter.",
                "Extend the variant path with a real HGVS parser and transcript selection when richer inputs appear.",
                "Add disease-ontology or MONDO-style crosswalks if multiple disease vocabularies must be merged.",
            ],
        },
    ]


def validate_metadata_outputs(metadata: dict, outdir: Path) -> None:
    for spec in metadata["deliverables"] + metadata.get("auxiliary_outputs", []):
        path = outdir / spec["path"]
        if not path.exists():
            raise AssertionError(f"Missing output: {path}")
        if spec["kind"] == "tsv":
            rows = read_tsv(path)
            missing = [column for column in spec.get("required_columns", []) if column not in (rows[0].keys() if rows else [])]
            if missing:
                raise AssertionError(f"Missing TSV columns in {path}: {missing}")
        elif spec["kind"] == "md":
            text = path.read_text(encoding="utf-8")
            for section in spec.get("required_sections", []):
                heading = f"## {section}"
                if heading not in text:
                    raise AssertionError(f"Missing section {heading} in {path}")
        elif spec["kind"] == "json":
            payload = load_json(path)
            missing = [key for key in spec.get("required_keys", []) if key not in payload]
            if missing:
                raise AssertionError(f"Missing JSON keys in {path}: {missing}")


def validate_expected_invariants(
    payload: dict,
    results_rows: list[dict],
    normalization_rows: list[dict],
    reconciliation_rows: list[dict],
    run_summary: dict,
) -> None:
    expected = payload.get("expected_invariants")
    if not expected:
        return

    if len(results_rows) != expected["query_count"]:
        raise AssertionError("Resolved query count does not match expected query count.")
    if run_summary["resolved_count"] != expected["resolved_count"]:
        raise AssertionError("Resolved count in run summary does not match expected value.")
    if run_summary["cross_source_supported_query_count"] < expected["minimum_cross_source_supported_queries"]:
        raise AssertionError("Cross-source support fell below the expected minimum.")

    expected_ambiguous = set(expected.get("ambiguous_queries", []))
    observed_ambiguous = {
        row["query"]
        for row in reconciliation_rows
        if row["ambiguity_flag"] == "yes"
    }
    if expected_ambiguous != observed_ambiguous:
        raise AssertionError(
            f"Ambiguous queries mismatch. expected={sorted(expected_ambiguous)} observed={sorted(observed_ambiguous)}"
        )

    result_index = {row["query"]: row for row in results_rows}
    normalization_index = {row["query"]: row for row in normalization_rows}
    reconciliation_index = {row["query"]: row for row in reconciliation_rows}

    for item in expected["expected_results"]:
        result_row = result_index[item["query"]]
        normalization_row = normalization_index[item["query"]]
        reconciliation_row = reconciliation_index[item["query"]]
        if result_row["normalized_id"] != item["normalized_id"]:
            raise AssertionError(f"Unexpected normalized identifier for {item['query']}.")
        if normalization_row["inferred_entity_type"] != item["inferred_entity_type"]:
            raise AssertionError(f"Unexpected entity type for {item['query']}.")
        observed_sources = set(filter(None, reconciliation_row["supporting_sources"].split(";")))
        if not set(item["supporting_sources"]).issubset(observed_sources):
            raise AssertionError(f"Missing supporting sources for {item['query']}.")


def run_skill(input_path: Path, outdir: Path) -> dict:
    payload = load_json(input_path)
    metadata = load_metadata()
    concepts = build_concept_index(payload)
    alias_index = build_alias_index(payload, concepts)
    source_indexes = build_source_indexes(payload)
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    normalization_rows: list[dict] = []
    fanout_rows: list[dict] = []
    reconciliation_rows: list[dict] = []
    results_rows: list[dict] = []

    for query in payload["queries"]:
        keys = candidate_keys(query)
        hits = collect_hits(query, keys, alias_index, source_indexes, payload["source_weights"])
        inferred_entity_type, type_scores, local_alias_hits = infer_entity_type(query, keys, hits, payload)
        winner, _ = reconcile_hits(query, inferred_entity_type, hits, concepts, payload)

        normalization_rows.append(
            {
                "query": query,
                "normalized_query": keys[0] if keys else "",
                "inferred_entity_type": inferred_entity_type,
                "candidate_keys": ";".join(keys),
                "local_alias_hits": str(local_alias_hits),
                "gene_type_score": f"{type_scores['gene']:.3f}",
                "protein_type_score": f"{type_scores['protein']:.3f}",
                "variant_type_score": f"{type_scores['variant']:.3f}",
                "disease_type_score": f"{type_scores['disease']:.3f}",
            }
        )

        for hit in hits:
            fanout_rows.append(
                {
                    "query": query,
                    "source": hit["source"],
                    "concept_id": hit["concept_id"],
                    "entity_type": hit["entity_type"],
                    "source_id": hit["source_id"],
                    "matched_key": hit["matched_key"],
                    "match_kind": hit["match_kind"],
                    "base_score": f"{float(hit['base_score']):.3f}",
                }
            )

        reconciliation_rows.append(
            {
                "query": query,
                "winning_concept_id": winner["concept_id"],
                "inferred_entity_type": inferred_entity_type,
                "resolved_id": winner["resolved_id"],
                "supporting_sources": ";".join(winner["supporting_sources"]),
                "source_hit_count": str(winner["source_hit_count"]),
                "winning_score": f"{winner['winning_score']:.3f}",
                "runner_up_concept_id": winner["runner_up_concept_id"],
                "score_margin": f"{winner['score_margin']:.3f}",
                "ambiguity_flag": winner["ambiguity_flag"],
            }
        )
        results_rows.append(build_query_result(query, winner))

    query_columns = metadata["deliverables"][0]["required_columns"]
    normalization_columns = metadata["auxiliary_outputs"][0]["required_columns"]
    fanout_columns = metadata["auxiliary_outputs"][1]["required_columns"]
    reconciliation_columns = metadata["auxiliary_outputs"][2]["required_columns"]

    write_tsv(outdir / "query_results.tsv", results_rows, query_columns)
    write_tsv(outdir / "normalization_qc.tsv", normalization_rows, normalization_columns)
    write_tsv(outdir / "fanout_qc.tsv", fanout_rows, fanout_columns)
    write_tsv(outdir / "reconciliation_qc.tsv", reconciliation_rows, reconciliation_columns)
    write_markdown(
        outdir / "source_provenance.md",
        "Source Provenance",
        build_source_provenance_sections(payload, normalization_rows, fanout_rows),
    )
    write_markdown(
        outdir / "resolution_notes.md",
        "Resolution Notes",
        build_resolution_sections(results_rows, reconciliation_rows),
    )

    source_hit_counts: dict[str, int] = defaultdict(int)
    for row in fanout_rows:
        source_hit_counts[row["source"]] += 1

    run_summary = {
        "run_label": payload["run_label"],
        "input_path": str(input_path.resolve()),
        "outdir": str(outdir),
        "query_count": len(payload["queries"]),
        "resolved_count": len(results_rows),
        "ambiguous_query_count": sum(1 for row in reconciliation_rows if row["ambiguity_flag"] == "yes"),
        "cross_source_supported_query_count": sum(
            1 for row in reconciliation_rows if int(row["source_hit_count"]) >= 2
        ),
        "source_hit_counts": dict(sorted(source_hit_counts.items())),
        "written_files": [],
    }
    write_json(outdir / "run_summary.json", run_summary)
    run_summary["written_files"] = sorted(
        str(path.relative_to(outdir))
        for path in outdir.rglob("*")
        if path.is_file()
    )
    write_json(outdir / "run_summary.json", run_summary)

    validate_metadata_outputs(metadata, outdir)
    validate_expected_invariants(payload, results_rows, normalization_rows, reconciliation_rows, run_summary)
    return run_summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the deterministic database-query starter.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        summary = run_skill(args.input, args.outdir)
    except Exception as exc:  # pragma: no cover - exercised through subprocess tests
        print(f"database-query starter failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
