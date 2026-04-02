#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = SKILL_DIR / "examples" / "toy_input.json"
VALID_BASES = set("ACGT")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_metadata() -> dict:
    return load_json(SKILL_DIR / "metadata.yaml")


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def gc_pct(sequence: str) -> float:
    gc_count = sequence.count("G") + sequence.count("C")
    return round((100.0 * gc_count) / len(sequence), 3)


def wallace_tm(sequence: str) -> float:
    return float((2 * (sequence.count("A") + sequence.count("T"))) + (4 * (sequence.count("G") + sequence.count("C"))))


def max_homopolymer_run(sequence: str) -> int:
    best = 0
    run = 0
    previous = ""
    for base in sequence:
        if base == previous:
            run += 1
        else:
            run = 1
            previous = base
        best = max(best, run)
    return best


def longest_exact_run(left: str, right: str) -> int:
    best = 0
    run = 0
    for left_base, right_base in zip(left, right):
        if left_base == right_base:
            run += 1
            best = max(best, run)
        else:
            run = 0
    return best


def max_end_self_complement(sequence: str, tail_length: int = 6) -> int:
    tail = sequence[-tail_length:]
    rc = reverse_complement(sequence)
    best = 0
    for start in range(0, len(rc) - len(tail) + 1):
        best = max(best, longest_exact_run(tail, rc[start : start + len(tail)]))
    return best


def max_cross_end_complement(left_primer: str, right_primer: str, tail_length: int = 6) -> int:
    left_tail = left_primer[-tail_length:]
    right_tail = right_primer[-tail_length:]
    right_rc = reverse_complement(right_primer)
    left_rc = reverse_complement(left_primer)
    best = 0
    for start in range(0, len(right_rc) - len(left_tail) + 1):
        best = max(best, longest_exact_run(left_tail, right_rc[start : start + len(left_tail)]))
    for start in range(0, len(left_rc) - len(right_tail) + 1):
        best = max(best, longest_exact_run(right_tail, left_rc[start : start + len(right_tail)]))
    return best


def normalize_sequence(sequence: str, *, label: str) -> str:
    normalized = sequence.strip().upper()
    invalid = sorted(set(normalized) - VALID_BASES)
    if invalid:
        raise ValueError(f"{label} contains unsupported bases: {', '.join(invalid)}")
    if not normalized:
        raise ValueError(f"{label} cannot be empty")
    return normalized


def validate_payload(payload: dict) -> None:
    if "deliverables" in payload:
        raise ValueError("Toy input must contain raw inputs and expected invariants, not precomputed deliverables.")

    required_keys = {
        "run_label",
        "lookup_policy",
        "primer_constraints",
        "queries",
        "references",
        "expected_invariants",
    }
    missing = sorted(required_keys - set(payload))
    if missing:
        raise ValueError(f"Toy input is missing required keys: {', '.join(missing)}")

    if not payload["queries"]:
        raise ValueError("At least one query is required.")
    if not payload["references"]:
        raise ValueError("At least one reference is required.")

    query_ids: set[str] = set()
    for query in payload["queries"]:
        query_id = query["query_id"]
        if query_id in query_ids:
            raise ValueError(f"Duplicate query_id: {query_id}")
        query_ids.add(query_id)
        normalize_sequence(query["sequence"], label=f"query {query_id}")

    reference_ids: set[str] = set()
    for reference in payload["references"]:
        reference_id = reference["reference_id"]
        if reference_id in reference_ids:
            raise ValueError(f"Duplicate reference_id: {reference_id}")
        reference_ids.add(reference_id)
        sequence = normalize_sequence(reference["sequence"], label=f"reference {reference_id}")
        for feature in reference["features"]:
            if feature["start"] < 1 or feature["end"] > len(sequence) or feature["start"] > feature["end"]:
                raise ValueError(f"Feature {feature['feature_id']} is out of bounds for {reference_id}")
        for variant in reference["variants"]:
            if variant["position"] < 1 or variant["position"] > len(sequence):
                raise ValueError(f"Variant {variant['variant_label']} is out of bounds for {reference_id}")
            if sequence[variant["position"] - 1] != variant["ref_base"]:
                raise ValueError(
                    f"Variant {variant['variant_label']} ref_base does not match reference sequence for {reference_id}"
                )
            if variant["ref_base"] == variant["alt_base"]:
                raise ValueError(f"Variant {variant['variant_label']} must change the base in {reference_id}")

    constraints = payload["primer_constraints"]
    min_length, max_length = constraints["primer_length_range"]
    if min_length > max_length:
        raise ValueError("primer_length_range must be ordered as [min, max].")
    min_product, max_product = constraints["product_size_range"]
    if min_product > max_product:
        raise ValueError("product_size_range must be ordered as [min, max].")


def find_best_alignment(query_sequence: str, reference_sequence: str, *, min_overlap: int) -> dict:
    orientations = {
        "+": query_sequence,
        "-": reverse_complement(query_sequence),
    }
    best: dict | None = None

    for orientation, oriented_query in orientations.items():
        for offset in range(-(len(oriented_query) - min_overlap), len(reference_sequence) - min_overlap + 1):
            query_start = max(0, -offset)
            reference_start = max(0, offset)
            overlap = min(len(oriented_query) - query_start, len(reference_sequence) - reference_start)
            if overlap < min_overlap:
                continue

            query_slice = oriented_query[query_start : query_start + overlap]
            reference_slice = reference_sequence[reference_start : reference_start + overlap]
            matches = sum(
                1 for query_base, reference_base in zip(query_slice, reference_slice) if query_base == reference_base
            )
            mismatches = overlap - matches
            longest_seed = longest_exact_run(query_slice, reference_slice)
            identity_pct = round((100.0 * matches) / overlap, 3)
            coverage_pct = round((100.0 * overlap) / len(oriented_query), 3)
            score = round((3.0 * matches) - (2.0 * mismatches) + longest_seed + (0.1 * overlap), 3)

            candidate = {
                "orientation": orientation,
                "oriented_query": oriented_query,
                "query_start": query_start,
                "query_end": query_start + overlap,
                "reference_start": reference_start,
                "reference_end": reference_start + overlap,
                "overlap_length": overlap,
                "matches": matches,
                "mismatches": mismatches,
                "identity_pct": identity_pct,
                "coverage_pct": coverage_pct,
                "longest_seed": longest_seed,
                "score": score,
            }

            ordering = (
                candidate["score"],
                candidate["identity_pct"],
                candidate["coverage_pct"],
                candidate["longest_seed"],
                -candidate["reference_start"],
                candidate["orientation"] == "+",
            )
            if best is None or ordering > (
                best["score"],
                best["identity_pct"],
                best["coverage_pct"],
                best["longest_seed"],
                -best["reference_start"],
                best["orientation"] == "+",
            ):
                best = candidate

    if best is None:
        raise ValueError("Failed to compute an alignment candidate.")
    return best


def overlap_rows_for_query(query: dict, reference: dict, alignment: dict) -> tuple[list[dict], list[dict]]:
    feature_rows: list[dict] = []
    variant_rows: list[dict] = []

    aligned_reference_start = alignment["reference_start"] + 1
    aligned_reference_end = alignment["reference_end"]
    oriented_query = alignment["oriented_query"]

    for feature in reference["features"]:
        overlap_start = max(feature["start"], aligned_reference_start)
        overlap_end = min(feature["end"], aligned_reference_end)
        if overlap_start > overlap_end:
            continue

        feature_rows.append(
            {
                "query_id": query["query_id"],
                "feature_source": feature["source"],
                "feature_type": feature["type"],
                "feature_id": feature["feature_id"],
                "start": feature["start"],
                "end": feature["end"],
                "label": feature["label"],
                "reference_id": reference["reference_id"],
                "reference_symbol": reference["symbol"],
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                "overlap_length": overlap_end - overlap_start + 1,
            }
        )

    for variant in reference["variants"]:
        if not aligned_reference_start <= variant["position"] <= aligned_reference_end:
            continue

        query_index = alignment["query_start"] + (variant["position"] - aligned_reference_start)
        observed_base = oriented_query[query_index]
        if observed_base == variant["alt_base"]:
            observed_state = "mutant"
        elif observed_base == variant["ref_base"]:
            observed_state = "reference"
        else:
            observed_state = "other"

        context_start = max(0, query_index - 4)
        context_end = min(len(oriented_query), query_index + 5)
        context_sequence = oriented_query[context_start:context_end]
        variant_row = {
            "query_id": query["query_id"],
            "matched_symbol": reference["symbol"],
            "variant_label": variant["variant_label"],
            "position": variant["position"],
            "reference_base": variant["ref_base"],
            "alternate_base": variant["alt_base"],
            "observed_base": observed_base,
            "observed_state": observed_state,
            "query_coordinate": query_index + 1,
            "context_sequence": context_sequence,
            "amplicon_target_label": variant["amplicon_target_label"],
        }
        variant_rows.append(variant_row)
        feature_rows.append(
            {
                "query_id": query["query_id"],
                "feature_source": "computed_variant",
                "feature_type": "variant_site",
                "feature_id": variant["variant_label"],
                "start": variant["position"],
                "end": variant["position"],
                "label": f"{variant['variant_label']} ({observed_state})",
                "reference_id": reference["reference_id"],
                "reference_symbol": reference["symbol"],
                "overlap_start": variant["position"],
                "overlap_end": variant["position"],
                "overlap_length": 1,
            }
        )

    return feature_rows, variant_rows


def evaluate_primer(sequence: str, *, optimal_tm: float) -> dict:
    gc_value = gc_pct(sequence)
    tm_value = wallace_tm(sequence)
    gc_clamp = sum(1 for base in sequence[-5:] if base in {"G", "C"})
    homopolymer = max_homopolymer_run(sequence)
    end_self_complement = max_end_self_complement(sequence)
    penalty = round(
        abs(tm_value - optimal_tm)
        + (abs(gc_value - 50.0) / 5.0)
        + (abs(len(sequence) - 20) * 0.5)
        + max(0, homopolymer - 4) * 2.0
        + max(0, end_self_complement - 3) * 2.5,
        3,
    )
    return {
        "sequence": sequence,
        "length": len(sequence),
        "gc_pct": gc_value,
        "tm": round(tm_value, 3),
        "gc_clamp": gc_clamp,
        "homopolymer": homopolymer,
        "end_self_complement": end_self_complement,
        "penalty": penalty,
    }


def design_primers_for_variant(
    query: dict,
    aligned_query: str,
    variant_row: dict,
    constraints: dict,
) -> tuple[list[dict], list[dict]]:
    min_length, max_length = constraints["primer_length_range"]
    gc_min, gc_max = constraints["gc_range"]
    tm_min, tm_max = constraints["tm_range"]
    product_min, product_max = constraints["product_size_range"]
    clamp_min, clamp_max = constraints["gc_clamp_range"]
    optimal_tm = float(constraints["optimal_tm"])
    flank_min = int(constraints["variant_flank_min"])
    max_homopolymer = int(constraints["max_homopolymer"])
    max_end_self_complement_allowed = int(constraints["max_end_self_complement"])
    max_tm_gap = float(constraints["max_tm_gap"])
    optimal_product = (product_min + product_max) / 2.0
    variant_index = int(variant_row["query_coordinate"]) - 1

    left_candidates: list[dict] = []
    right_candidates: list[dict] = []
    qc_rows: list[dict] = []
    accepted_pairs: list[dict] = []

    def primer_passes(primer: dict) -> bool:
        return (
            gc_min <= primer["gc_pct"] <= gc_max
            and tm_min <= primer["tm"] <= tm_max
            and clamp_min <= primer["gc_clamp"] <= clamp_max
            and primer["homopolymer"] <= max_homopolymer
            and primer["end_self_complement"] <= max_end_self_complement_allowed
        )

    for start in range(0, len(aligned_query)):
        for length in range(min_length, max_length + 1):
            end = start + length
            if end > len(aligned_query):
                continue
            if variant_index - (end - 1) < flank_min:
                continue
            forward_primer = evaluate_primer(aligned_query[start:end], optimal_tm=optimal_tm)
            forward_primer["start"] = start + 1
            forward_primer["end"] = end
            if primer_passes(forward_primer):
                left_candidates.append(forward_primer)

    for start in range(0, len(aligned_query)):
        for length in range(min_length, max_length + 1):
            end = start + length
            if end > len(aligned_query):
                continue
            if start - variant_index < flank_min:
                continue
            reverse_primer = evaluate_primer(reverse_complement(aligned_query[start:end]), optimal_tm=optimal_tm)
            reverse_primer["start"] = start + 1
            reverse_primer["end"] = end
            if primer_passes(reverse_primer):
                right_candidates.append(reverse_primer)

    for left in left_candidates:
        for right in right_candidates:
            product_size = right["end"] - left["start"] + 1
            if not product_min <= product_size <= product_max:
                continue

            tm_gap = round(abs(left["tm"] - right["tm"]), 3)
            if tm_gap > max_tm_gap:
                continue

            cross_end = max_cross_end_complement(left["sequence"], right["sequence"])
            pair_score = round(
                left["penalty"]
                + right["penalty"]
                + tm_gap
                + (abs(product_size - optimal_product) / 6.0)
                + (cross_end * 1.5),
                3,
            )
            pair = {
                "query_id": query["query_id"],
                "amplicon_target_label": variant_row["amplicon_target_label"],
                "variant_label": variant_row["variant_label"],
                "left_primer": left["sequence"],
                "right_primer": right["sequence"],
                "left_start": left["start"],
                "left_end": left["end"],
                "right_start": right["start"],
                "right_end": right["end"],
                "left_tm": left["tm"],
                "right_tm": right["tm"],
                "left_gc_pct": left["gc_pct"],
                "right_gc_pct": right["gc_pct"],
                "product_size": product_size,
                "tm_gap": tm_gap,
                "cross_end_complement": cross_end,
                "pair_score": pair_score,
                "variant_offset_from_left_end": variant_index - (left["end"] - 1),
                "variant_offset_from_right_start": right["start"] - (variant_index + 1),
            }
            accepted_pairs.append(pair)

    accepted_pairs.sort(
        key=lambda row: (
            row["pair_score"],
            row["product_size"],
            row["left_start"],
            row["right_start"],
            row["left_primer"],
            row["right_primer"],
        )
    )

    for rank, pair in enumerate(accepted_pairs, start=1):
        qc_row = dict(pair)
        qc_row["candidate_rank"] = rank
        qc_rows.append(qc_row)

    if not accepted_pairs:
        raise ValueError(f"No acceptable primer pairs found for {query['query_id']}")

    deliverable_rows: list[dict] = []
    for pair_rank, pair in enumerate(accepted_pairs[: constraints["num_return"]], start=1):
        deliverable_rows.append(
            {
                "query_id": pair["query_id"],
                "pair_rank": pair_rank,
                "left_primer": pair["left_primer"],
                "right_primer": pair["right_primer"],
                "product_size": pair["product_size"],
                "amplicon_target_label": pair["amplicon_target_label"],
                "left_start": pair["left_start"],
                "right_end": pair["right_end"],
                "left_tm": pair["left_tm"],
                "right_tm": pair["right_tm"],
                "left_gc_pct": pair["left_gc_pct"],
                "right_gc_pct": pair["right_gc_pct"],
                "pair_score": pair["pair_score"],
            }
        )

    return deliverable_rows, qc_rows


def analyze_queries(payload: dict) -> dict:
    metadata = load_metadata()
    lookup_rows: list[dict] = []
    sequence_rows: list[dict] = []
    feature_rows: list[dict] = []
    variant_rows: list[dict] = []
    primer_rows: list[dict] = []
    primer_qc_rows: list[dict] = []
    matched_oriented_queries: dict[str, str] = {}

    min_overlap = int(payload["lookup_policy"]["min_overlap"])

    for query in payload["queries"]:
        query_sequence = normalize_sequence(query["sequence"], label=f"query {query['query_id']}")
        candidate_rows: list[dict] = []

        for reference in payload["references"]:
            reference_sequence = normalize_sequence(reference["sequence"], label=f"reference {reference['reference_id']}")
            alignment = find_best_alignment(query_sequence, reference_sequence, min_overlap=min_overlap)
            candidate_rows.append(
                {
                    "query_id": query["query_id"],
                    "reference_id": reference["reference_id"],
                    "matched_symbol": reference["symbol"],
                    "matched_accession": reference["accession"],
                    "matched_species": reference["species"],
                    **alignment,
                }
            )

        candidate_rows.sort(
            key=lambda row: (
                -row["score"],
                -row["identity_pct"],
                -row["coverage_pct"],
                -row["longest_seed"],
                row["matched_symbol"],
            )
        )
        for rank, candidate in enumerate(candidate_rows, start=1):
            lookup_rows.append(
                {
                    "query_id": candidate["query_id"],
                    "candidate_rank": rank,
                    "reference_id": candidate["reference_id"],
                    "matched_symbol": candidate["matched_symbol"],
                    "matched_accession": candidate["matched_accession"],
                    "matched_species": candidate["matched_species"],
                    "orientation": candidate["orientation"],
                    "overlap_length": candidate["overlap_length"],
                    "matches": candidate["matches"],
                    "mismatches": candidate["mismatches"],
                    "identity_pct": candidate["identity_pct"],
                    "coverage_pct": candidate["coverage_pct"],
                    "longest_seed": candidate["longest_seed"],
                    "aligned_reference_start": candidate["reference_start"] + 1,
                    "aligned_reference_end": candidate["reference_end"],
                    "score": candidate["score"],
                }
            )

        best = candidate_rows[0]
        if best["identity_pct"] < float(payload["lookup_policy"]["minimum_identity_pct"]):
            raise ValueError(f"Best hit for {query['query_id']} fell below the minimum identity threshold")

        reference = next(item for item in payload["references"] if item["reference_id"] == best["reference_id"])
        aligned_feature_rows, aligned_variant_rows = overlap_rows_for_query(query, reference, best)
        if not aligned_variant_rows:
            raise ValueError(f"No tracked variant was covered by the winning alignment for {query['query_id']}")

        matched_oriented_queries[query["query_id"]] = best["oriented_query"]
        feature_rows.extend(aligned_feature_rows)
        variant_rows.extend(aligned_variant_rows)

        primary_variant = sorted(
            aligned_variant_rows,
            key=lambda row: (
                {"mutant": 0, "reference": 1, "other": 2}[row["observed_state"]],
                row["position"],
                row["variant_label"],
            ),
        )[0]

        sequence_rows.append(
            {
                "query_id": query["query_id"],
                "matched_symbol": best["matched_symbol"],
                "matched_accession": best["matched_accession"],
                "matched_species": best["matched_species"],
                "variant_label": primary_variant["variant_label"],
                "orientation": best["orientation"],
                "identity_pct": best["identity_pct"],
                "coverage_pct": best["coverage_pct"],
                "aligned_reference_start": best["reference_start"] + 1,
                "aligned_reference_end": best["reference_end"],
                "variant_state": primary_variant["observed_state"],
                "query_length": len(query_sequence),
            }
        )

    variant_rows.sort(key=lambda row: (row["query_id"], row["position"], row["variant_label"]))
    for variant_row in variant_rows:
        aligned_query = matched_oriented_queries[variant_row["query_id"]]
        query = next(item for item in payload["queries"] if item["query_id"] == variant_row["query_id"])
        deliverable_pairs, qc_pairs = design_primers_for_variant(
            query,
            aligned_query,
            variant_row,
            payload["primer_constraints"],
        )
        primer_rows.extend(deliverable_pairs)
        primer_qc_rows.extend(qc_pairs)

    output = {
        "metadata": metadata,
        "sequence_rows": sorted(sequence_rows, key=lambda row: row["query_id"]),
        "feature_rows": sorted(
            feature_rows,
            key=lambda row: (row["query_id"], int(row["start"]), int(row["end"]), row["feature_id"]),
        ),
        "primer_rows": sorted(
            primer_rows,
            key=lambda row: (row["query_id"], int(row["pair_rank"]), row["left_primer"], row["right_primer"]),
        ),
        "lookup_rows": sorted(
            lookup_rows,
            key=lambda row: (row["query_id"], int(row["candidate_rank"]), row["matched_symbol"]),
        ),
        "variant_rows": variant_rows,
        "primer_qc_rows": sorted(
            primer_qc_rows,
            key=lambda row: (row["query_id"], float(row["pair_score"]), int(row["candidate_rank"])),
        ),
    }
    return output


def write_outputs(payload: dict, outdir: Path) -> dict:
    analysis = analyze_queries(payload)
    metadata = analysis["metadata"]
    outdir.mkdir(parents=True, exist_ok=True)

    deliverable_columns = {
        item["path"]: list(item["required_columns"])
        for item in metadata["deliverables"]
    }
    deliverable_columns["sequence_summary.tsv"].extend(
        [
            "orientation",
            "identity_pct",
            "coverage_pct",
            "aligned_reference_start",
            "aligned_reference_end",
            "variant_state",
            "query_length",
        ]
    )
    deliverable_columns["feature_annotations.tsv"].extend(
        [
            "reference_id",
            "reference_symbol",
            "overlap_start",
            "overlap_end",
            "overlap_length",
        ]
    )
    deliverable_columns["primer_candidates.tsv"].extend(
        [
            "left_start",
            "right_end",
            "left_tm",
            "right_tm",
            "left_gc_pct",
            "right_gc_pct",
            "pair_score",
        ]
    )

    write_tsv(outdir / "sequence_summary.tsv", analysis["sequence_rows"], deliverable_columns["sequence_summary.tsv"])
    write_tsv(outdir / "feature_annotations.tsv", analysis["feature_rows"], deliverable_columns["feature_annotations.tsv"])
    write_tsv(outdir / "primer_candidates.tsv", analysis["primer_rows"], deliverable_columns["primer_candidates.tsv"])

    write_tsv(
        outdir / "lookup_qc.tsv",
        analysis["lookup_rows"],
        [
            "query_id",
            "candidate_rank",
            "reference_id",
            "matched_symbol",
            "matched_accession",
            "matched_species",
            "orientation",
            "overlap_length",
            "matches",
            "mismatches",
            "identity_pct",
            "coverage_pct",
            "longest_seed",
            "aligned_reference_start",
            "aligned_reference_end",
            "score",
        ],
    )
    write_tsv(
        outdir / "variant_qc.tsv",
        analysis["variant_rows"],
        [
            "query_id",
            "matched_symbol",
            "variant_label",
            "position",
            "reference_base",
            "alternate_base",
            "observed_base",
            "observed_state",
            "query_coordinate",
            "context_sequence",
            "amplicon_target_label",
        ],
    )
    write_tsv(
        outdir / "primer_qc.tsv",
        analysis["primer_qc_rows"],
        [
            "query_id",
            "candidate_rank",
            "amplicon_target_label",
            "variant_label",
            "left_primer",
            "right_primer",
            "left_start",
            "left_end",
            "right_start",
            "right_end",
            "left_tm",
            "right_tm",
            "left_gc_pct",
            "right_gc_pct",
            "product_size",
            "tm_gap",
            "cross_end_complement",
            "pair_score",
            "variant_offset_from_left_end",
            "variant_offset_from_right_start",
        ],
    )

    written_files = sorted(
        [
            "feature_annotations.tsv",
            "lookup_qc.tsv",
            "primer_candidates.tsv",
            "primer_qc.tsv",
            "run_summary.json",
            "sequence_summary.tsv",
            "variant_qc.tsv",
        ]
    )
    summary = {
        "run_label": payload["run_label"],
        "query_count": len(payload["queries"]),
        "resolved_query_count": len(analysis["sequence_rows"]),
        "reverse_complement_query_count": sum(1 for row in analysis["sequence_rows"] if row["orientation"] == "-"),
        "mutant_variant_count": sum(1 for row in analysis["variant_rows"] if row["observed_state"] == "mutant"),
        "feature_annotation_count": len(analysis["feature_rows"]),
        "primer_pair_count": len(analysis["primer_rows"]),
        "written_files": written_files,
    }
    write_json(outdir / "run_summary.json", summary)
    return summary


def validate_outputs(outdir: Path, *, input_path: Path = DEFAULT_INPUT) -> dict:
    payload = load_json(input_path)
    metadata = load_metadata()
    expected = payload["expected_invariants"]

    for deliverable in metadata["deliverables"]:
        path = outdir / deliverable["path"]
        if not path.exists():
            raise AssertionError(f"Missing deliverable: {path}")
        rows = read_tsv(path)
        if not rows:
            raise AssertionError(f"Deliverable is empty: {path}")
        missing_columns = [column for column in deliverable["required_columns"] if column not in rows[0]]
        if missing_columns:
            raise AssertionError(f"Missing required columns in {path}: {missing_columns}")

    for extra_name in ("lookup_qc.tsv", "variant_qc.tsv", "primer_qc.tsv", "run_summary.json"):
        if not (outdir / extra_name).exists():
            raise AssertionError(f"Missing QC artifact: {outdir / extra_name}")

    sequence_rows = read_tsv(outdir / "sequence_summary.tsv")
    feature_rows = read_tsv(outdir / "feature_annotations.tsv")
    primer_rows = read_tsv(outdir / "primer_candidates.tsv")
    lookup_rows = read_tsv(outdir / "lookup_qc.tsv")
    variant_rows = read_tsv(outdir / "variant_qc.tsv")
    primer_qc_rows = read_tsv(outdir / "primer_qc.tsv")
    summary = load_json(outdir / "run_summary.json")

    if len(sequence_rows) != expected["query_count"]:
        raise AssertionError("Unexpected number of sequence_summary rows")
    if summary["resolved_query_count"] != expected["query_count"]:
        raise AssertionError("Resolved query count does not match expected query count")
    if summary["reverse_complement_query_count"] != expected["reverse_complement_query_count"]:
        raise AssertionError("Unexpected reverse complement match count")
    if summary["mutant_variant_count"] != expected["mutant_variant_count"]:
        raise AssertionError("Unexpected mutant variant count")
    if len(feature_rows) < expected["minimum_feature_annotation_count"]:
        raise AssertionError("Feature annotations did not meet the minimum row count")
    if len(primer_rows) != expected["query_count"] * expected["primer_pairs_per_query"]:
        raise AssertionError("Unexpected number of primer candidate rows")
    if len(lookup_rows) != expected["query_count"] * expected["lookup_candidates_per_query"]:
        raise AssertionError("Unexpected number of lookup QC rows")
    if len(variant_rows) != expected["query_count"]:
        raise AssertionError("Unexpected number of variant QC rows")
    if len(primer_qc_rows) < len(primer_rows):
        raise AssertionError("Primer QC rows should cover at least the reported primer candidates")

    sequence_index = {row["query_id"]: row for row in sequence_rows}
    variant_index = {row["query_id"]: row for row in variant_rows}
    top_primer = {row["query_id"]: row for row in primer_rows if row["pair_rank"] == "1"}

    product_min, product_max = payload["primer_constraints"]["product_size_range"]

    for query_id, symbol in expected["matched_symbols"].items():
        if sequence_index[query_id]["matched_symbol"] != symbol:
            raise AssertionError(f"Unexpected matched symbol for {query_id}")
        if sequence_index[query_id]["matched_accession"] != expected["matched_accessions"][query_id]:
            raise AssertionError(f"Unexpected matched accession for {query_id}")
        if sequence_index[query_id]["orientation"] != expected["orientation_by_query"][query_id]:
            raise AssertionError(f"Unexpected orientation for {query_id}")
        if float(sequence_index[query_id]["identity_pct"]) < expected["minimum_identity_pct"]:
            raise AssertionError(f"Identity below minimum threshold for {query_id}")
        if variant_index[query_id]["observed_state"] != expected["variant_state_by_query"][query_id]:
            raise AssertionError(f"Unexpected variant state for {query_id}")
        if int(top_primer[query_id]["product_size"]) != expected["top_product_size_by_query"][query_id]:
            raise AssertionError(f"Unexpected top primer product size for {query_id}")
        if top_primer[query_id]["amplicon_target_label"] != expected["amplicon_target_by_query"][query_id]:
            raise AssertionError(f"Unexpected amplicon target label for {query_id}")

    for row in primer_rows:
        product_size = int(row["product_size"])
        if not product_min <= product_size <= product_max:
            raise AssertionError(f"Primer product size fell outside the configured range: {product_size}")

    return {
        "validated_queries": sorted(sequence_index),
        "feature_annotation_count": len(feature_rows),
        "primer_pair_count": len(primer_rows),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the sequence-analysis portable starter.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    payload = load_json(args.input)
    validate_payload(payload)
    summary = write_outputs(payload, args.outdir.resolve())
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - surfaced in CLI tests
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
