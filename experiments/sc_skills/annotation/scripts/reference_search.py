#!/usr/bin/env python3
from __future__ import annotations

import math
from typing import Any

import pandas as pd


NORMAL_LIKE_TERMS = {
    "adjacent-normal",
    "adjacent normal",
    "baseline",
    "control",
    "healthy",
    "homeostatic",
    "normal",
    "non-diseased",
    "steady state",
    "steady-state",
    "unaffected",
}

SPECIES_ALIASES = {
    "human": {"human", "homo sapiens"},
    "mouse": {"mouse", "mus musculus"},
}


def normalize_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    return " ".join(str(value).strip().lower().replace("_", " ").replace("-", " ").split())


def tokenize(*values: Any) -> set[str]:
    tokens: set[str] = set()
    for value in values:
        normalized = normalize_text(value)
        if normalized:
            tokens.update(normalized.split())
    return tokens


def build_reference_query_string(species: str, tissue: str, condition: str) -> str:
    return " ".join(part for part in [normalize_text(species), normalize_text(tissue), normalize_text(condition)] if part)


def build_reference_query_variants(species: str, tissue: str, condition: str) -> list[str]:
    base = build_reference_query_string(species, tissue, condition)
    variants = [base]
    normalized_condition = normalize_text(condition)
    if normalized_condition in NORMAL_LIKE_TERMS:
        variants.extend(
            [
                build_reference_query_string(species, tissue, "healthy"),
                build_reference_query_string(species, tissue, "control"),
            ]
        )
    variants.append(build_reference_query_string(species, tissue, "atlas"))
    deduped: list[str] = []
    for item in variants:
        if item and item not in deduped:
            deduped.append(item)
    return deduped


def _coerce_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    normalized = normalize_text(value)
    return normalized in {"1", "true", "t", "yes", "y", "present", "available"}


def _coerce_int(value: Any) -> int:
    if value in {None, ""}:
        return 0
    try:
        if pd.isna(value):
            return 0
    except TypeError:
        pass
    return int(float(value))


def _species_match_score(query_species: str, candidate_species: str) -> tuple[bool, float, list[str], list[str]]:
    notes: list[str] = []
    rejections: list[str] = []
    query_norm = normalize_text(query_species)
    candidate_norm = normalize_text(candidate_species)
    if not query_norm or not candidate_norm:
        notes.append("species_missing")
        return False, 0.0, notes, rejections

    query_aliases = SPECIES_ALIASES.get(query_norm, {query_norm})
    candidate_aliases = SPECIES_ALIASES.get(candidate_norm, {candidate_norm})
    if query_aliases & candidate_aliases:
        notes.append("species_exact")
        return True, 6.0, notes, rejections

    rejections.append(f"species_mismatch:{candidate_species}")
    return False, -6.0, notes, rejections


def _tissue_match_score(query_tissue: str, candidate_tissue: str, title: str) -> tuple[bool, float, list[str], list[str]]:
    notes: list[str] = []
    rejections: list[str] = []
    query_norm = normalize_text(query_tissue)
    candidate_norm = normalize_text(candidate_tissue)
    if not query_norm:
        notes.append("tissue_missing_from_query")
        return False, 0.0, notes, rejections
    if candidate_norm == query_norm:
        notes.append("tissue_exact")
        return True, 5.0, notes, rejections

    overlap = tokenize(query_norm) & tokenize(candidate_norm, title)
    if overlap:
        notes.append("tissue_partial")
        return True, 2.0, notes, rejections

    rejections.append(f"tissue_mismatch:{candidate_tissue}")
    return False, -4.0, notes, rejections


def _condition_match_score(query_condition: str, candidate_disease: str) -> tuple[bool, float, list[str], list[str]]:
    notes: list[str] = []
    rejections: list[str] = []
    query_norm = normalize_text(query_condition)
    candidate_norm = normalize_text(candidate_disease)
    if not query_norm:
        notes.append("condition_unspecified")
        return False, 0.0, notes, rejections

    if query_norm == candidate_norm:
        notes.append("condition_exact")
        return True, 4.0, notes, rejections

    if query_norm in NORMAL_LIKE_TERMS:
        if candidate_norm in NORMAL_LIKE_TERMS:
            notes.append("healthy_like")
            return True, 3.5, notes, rejections
        if not candidate_norm:
            notes.append("disease_missing")
            return False, 0.0, notes, rejections
        rejections.append(f"disease_not_normal:{candidate_disease}")
        return False, -3.5, notes, rejections

    if query_norm and query_norm in candidate_norm:
        notes.append("condition_partial")
        return True, 2.0, notes, rejections

    if candidate_norm:
        notes.append("condition_mismatch")
        return False, -1.5, notes, rejections
    return False, 0.0, notes, rejections


def _annotation_score(record: dict[str, Any]) -> tuple[bool, bool, float, list[str], list[str]]:
    notes: list[str] = []
    rejections: list[str] = []
    broad_column = normalize_text(record.get("cell_type_obs_key"))
    if not broad_column:
        broad_column = normalize_text(record.get("cell_type_column"))
    subtype_column = normalize_text(record.get("subtype_obs_key"))
    if not subtype_column:
        subtype_column = normalize_text(record.get("subtype_column"))
    has_broad = _coerce_bool(record.get("has_broad_labels")) or bool(broad_column)
    has_subtype = _coerce_bool(record.get("has_subtype_labels")) or bool(subtype_column)
    score = 0.0
    if has_broad:
        score += 2.0
        notes.append("broad_labels_present")
    else:
        rejections.append("missing_broad_labels")
    if has_subtype:
        score += 1.5
        notes.append("subtype_labels_present")
    else:
        notes.append("subtype_labels_missing")
    return has_broad, has_subtype, score, notes, rejections


def _cell_count_score(cell_count: int, min_cells: int) -> tuple[float, list[str], list[str]]:
    notes: list[str] = []
    rejections: list[str] = []
    if cell_count < int(min_cells):
        rejections.append(f"too_few_cells:{cell_count}")
    if cell_count > 0:
        notes.append("sufficient_scale" if cell_count >= int(min_cells) else "low_scale")
    score = min(float(cell_count) / 10000.0, 2.0)
    return score, notes, rejections


def rank_reference_candidates(records: list[dict[str, Any]], query_context: dict[str, Any]) -> pd.DataFrame:
    query_string = build_reference_query_string(
        str(query_context.get("species", "")),
        str(query_context.get("tissue", "")),
        str(query_context.get("condition", "")),
    )
    query_tokens = tokenize(query_string)
    min_cells = int(query_context.get("min_cells", 1000))
    rows: list[dict[str, Any]] = []
    for original in records:
        record = dict(original)
        dataset_id = str(record.get("dataset_id", ""))
        dataset_title = str(record.get("dataset_title", ""))
        collection_name = str(record.get("collection_name", ""))
        species = str(record.get("species") or record.get("organism") or "")
        tissue = str(record.get("tissue") or record.get("tissue_general") or "")
        disease = str(record.get("disease") or record.get("condition") or "")
        assay = str(record.get("assay") or record.get("platform_name") or "")
        suspension_type = str(record.get("suspension_type") or "")
        cell_count = _coerce_int(record.get("dataset_total_cell_count", record.get("cell_count", 0)))
        searchable_tokens = tokenize(dataset_title, collection_name, species, tissue, disease, assay)
        overlap = len(query_tokens & searchable_tokens)

        species_match, species_score, species_notes, species_rejections = _species_match_score(
            str(query_context.get("species", "")),
            species,
        )
        tissue_match, tissue_score, tissue_notes, tissue_rejections = _tissue_match_score(
            str(query_context.get("tissue", "")),
            tissue,
            dataset_title,
        )
        condition_match, condition_score, condition_notes, condition_rejections = _condition_match_score(
            str(query_context.get("condition", "")),
            disease,
        )
        has_broad, has_subtype, annotation_score, annotation_notes, annotation_rejections = _annotation_score(record)
        cell_count_score, cell_count_notes, cell_count_rejections = _cell_count_score(cell_count, min_cells)

        notes = species_notes + tissue_notes + condition_notes + annotation_notes + cell_count_notes
        rejections = species_rejections + tissue_rejections + condition_rejections + annotation_rejections + cell_count_rejections

        suspension_score = 0.5 if normalize_text(suspension_type) == "cell" else -1.5
        if normalize_text(suspension_type) not in {"", "cell"}:
            rejections.append(f"suspension_not_cell:{suspension_type}")

        primary_data_score = 0.25 if _coerce_bool(record.get("is_primary_data", True)) else 0.0
        token_overlap_score = overlap * 0.15
        total_score = (
            species_score
            + tissue_score
            + condition_score
            + annotation_score
            + cell_count_score
            + suspension_score
            + primary_data_score
            + token_overlap_score
        )

        row = dict(record)
        row.update(
            {
                "dataset_id": dataset_id,
                "dataset_title": dataset_title,
                "collection_name": collection_name,
                "species": species,
                "tissue": tissue,
                "disease": disease,
                "assay": assay,
                "suspension_type": suspension_type,
                "cell_count": cell_count,
                "dataset_total_cell_count": cell_count,
                "query_string": query_string,
                "query_token_overlap": overlap,
                "species_match": species_match,
                "tissue_match": tissue_match,
                "condition_match": condition_match,
                "has_broad_labels": has_broad,
                "has_subtype_labels": has_subtype,
                "species_score": round(species_score, 6),
                "tissue_score": round(tissue_score, 6),
                "condition_score": round(condition_score, 6),
                "annotation_score": round(annotation_score, 6),
                "cell_count_score": round(cell_count_score, 6),
                "token_overlap_score": round(token_overlap_score, 6),
                "suspension_score": round(suspension_score, 6),
                "primary_data_score": round(primary_data_score, 6),
                "score": round(total_score, 6),
                "accept_for_transfer": len(rejections) == 0,
                "selection_notes": ";".join(notes) if notes else "none",
                "rejection_reasons": ";".join(sorted(dict.fromkeys(rejections))) if rejections else "none",
            }
        )
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame = frame.sort_values(
        by=["accept_for_transfer", "score", "has_subtype_labels", "dataset_total_cell_count", "dataset_id"],
        ascending=[False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    frame.insert(0, "rank", list(range(1, len(frame) + 1)))
    frame["selected"] = False
    accepted = frame.index[frame["accept_for_transfer"]].tolist()
    if accepted:
        frame.loc[accepted[0], "selected"] = True
    return frame


def build_reference_selection_summary(frame: pd.DataFrame, query_context: dict[str, Any]) -> dict[str, Any]:
    accepted = frame.loc[frame["accept_for_transfer"]].copy()
    selected = frame.loc[frame["selected"]].copy()
    if selected.empty:
        return {
            "query_string": build_reference_query_string(
                str(query_context.get("species", "")),
                str(query_context.get("tissue", "")),
                str(query_context.get("condition", "")),
            ),
            "query_variants": build_reference_query_variants(
                str(query_context.get("species", "")),
                str(query_context.get("tissue", "")),
                str(query_context.get("condition", "")),
            ),
            "candidate_count": int(frame.shape[0]),
            "accepted_candidate_count": int(accepted.shape[0]),
            "selected_dataset_id": None,
            "selected_dataset_title": None,
            "selected_score": None,
            "stop_reason": "No candidate passed the hard reference gates. Stop and curate a better reference.",
            "top_candidate_rejections": [
                {
                    "dataset_id": str(row["dataset_id"]),
                    "rejection_reasons": str(row["rejection_reasons"]),
                }
                for _, row in frame.head(3).iterrows()
                if str(row["rejection_reasons"]) != "none"
            ],
        }

    row = selected.iloc[0]
    return {
        "query_string": str(row["query_string"]),
        "query_variants": build_reference_query_variants(
            str(query_context.get("species", "")),
            str(query_context.get("tissue", "")),
            str(query_context.get("condition", "")),
        ),
        "candidate_count": int(frame.shape[0]),
        "accepted_candidate_count": int(accepted.shape[0]),
        "selected_dataset_id": str(row["dataset_id"]),
        "selected_dataset_title": str(row["dataset_title"]),
        "selected_score": float(row["score"]),
        "selected_notes": str(row["selection_notes"]),
        "selected_rejection_reasons": str(row["rejection_reasons"]),
        "selected_rank": int(row["rank"]),
        "stop_reason": None,
    }
