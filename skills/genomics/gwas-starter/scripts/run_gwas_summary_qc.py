#!/usr/bin/env python3
"""Run deterministic GWAS summary-statistics QC and interpretation triage."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from statistics import median


DEFAULT_CONFIG = {
    "min_info": 0.8,
    "min_minor_allele_frequency": 0.01,
    "low_sample_size_ratio": 0.67,
    "high_sample_size_ratio": 1.5,
    "ambiguous_palindromic_eaf_lower": 0.4,
    "ambiguous_palindromic_eaf_upper": 0.6,
    "genome_wide_significance": 5e-8,
    "suggestive_significance": 1e-5,
    "lead_variant_limit": 5,
}

FIELD_ALIASES = {
    "variant_id": ["variant_id", "variant", "snp", "rsid", "markername", "marker", "id"],
    "chromosome": ["chromosome", "chr", "chrom", "#chrom"],
    "position": ["position", "pos", "bp", "base_pair_location"],
    "effect_allele": ["effect_allele", "ea", "a1", "alt", "effectallele"],
    "other_allele": ["other_allele", "nea", "a2", "ref", "non_effect_allele"],
    "beta": ["beta", "effect", "beta_gc"],
    "odds_ratio": ["or", "odds_ratio"],
    "se": ["se", "stderr", "standard_error"],
    "p": ["p", "pval", "p_value", "pvalue"],
    "n": ["n", "n_total", "samplesize"],
    "eaf": ["eaf", "effect_allele_frequency", "frq", "freq1", "af"],
    "info": ["info", "imputation_info", "info_score"],
}

BASE_COLUMNS = [
    "variant_id",
    "chromosome",
    "position",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    "n",
    "eaf",
    "info",
    "qc_status",
    "include_in_downstream",
    "qc_flags",
]


def load_config(path: Path | None) -> dict:
    config = dict(DEFAULT_CONFIG)
    if path is None:
        return config
    loaded = json.loads(path.read_text(encoding="utf-8"))
    config.update(loaded)
    return config


def detect_delimiter(path: Path) -> str:
    header = path.read_text(encoding="utf-8").splitlines()[0]
    return "\t" if header.count("\t") >= header.count(",") else ","


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    delimiter = detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"No rows found in {path}.")
    return rows, fieldnames


def resolve_column_map(fieldnames: list[str]) -> dict[str, str]:
    normalized = {name.lower(): name for name in fieldnames}
    mapping: dict[str, str] = {}
    for canonical, aliases in FIELD_ALIASES.items():
        for alias in aliases:
            if alias.lower() in normalized:
                mapping[canonical] = normalized[alias.lower()]
                break
    required = ["chromosome", "position", "effect_allele", "other_allele", "p"]
    missing = [field for field in required if field not in mapping]
    if missing:
        raise KeyError(f"Missing required GWAS columns: {', '.join(missing)}")
    if "beta" not in mapping and "odds_ratio" not in mapping:
        raise KeyError("Missing effect-size column: provide BETA or OR.")
    return mapping


def parse_float(value: str, field: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric value for {field}: {value!r}") from exc


def parse_intish(value: str, field: str) -> int:
    try:
        return int(float(value))
    except ValueError as exc:
        raise ValueError(f"Invalid integer-like value for {field}: {value!r}") from exc


def is_palindromic(effect_allele: str, other_allele: str) -> bool:
    pair = {effect_allele.upper(), other_allele.upper()}
    return pair in [{"A", "T"}, {"C", "G"}]


def normalize_variant(row: dict[str, str], mapping: dict[str, str], config: dict) -> dict[str, object]:
    flags: list[str] = []

    chromosome = row[mapping["chromosome"]]
    position = parse_intish(row[mapping["position"]], "position")
    effect_allele = row[mapping["effect_allele"]].upper()
    other_allele = row[mapping["other_allele"]].upper()
    original_variant_id = row.get(mapping.get("variant_id", ""), "").strip()
    variant_id = original_variant_id or f"{chromosome}:{position}:{other_allele}:{effect_allele}"

    p_value = parse_float(row[mapping["p"]], "p")
    if not (0.0 < p_value <= 1.0):
        flags.append("invalid_p")

    beta = None
    if "beta" in mapping and row[mapping["beta"]]:
        beta = parse_float(row[mapping["beta"]], "beta")
    elif "odds_ratio" in mapping and row[mapping["odds_ratio"]]:
        odds_ratio = parse_float(row[mapping["odds_ratio"]], "odds_ratio")
        if odds_ratio <= 0.0:
            flags.append("invalid_odds_ratio")
        else:
            beta = math.log(odds_ratio)
    else:
        flags.append("missing_effect_size")

    se = None
    if "se" in mapping and row[mapping["se"]]:
        se = parse_float(row[mapping["se"]], "se")
        if se <= 0.0:
            flags.append("invalid_se")

    n = None
    if "n" in mapping and row[mapping["n"]]:
        n = parse_intish(row[mapping["n"]], "n")
        if n <= 0:
            flags.append("invalid_n")

    eaf = None
    if "eaf" in mapping and row[mapping["eaf"]]:
        eaf = parse_float(row[mapping["eaf"]], "eaf")
        if not (0.0 <= eaf <= 1.0):
            flags.append("invalid_eaf")

    info = None
    if "info" in mapping and row[mapping["info"]]:
        info = parse_float(row[mapping["info"]], "info")
        if not (0.0 <= info <= 1.0):
            flags.append("invalid_info")

    if info is not None and info < float(config["min_info"]):
        flags.append("low_info")
    if eaf is not None:
        maf = min(eaf, 1.0 - eaf)
        if maf < float(config["min_minor_allele_frequency"]):
            flags.append("low_maf")
        lower = float(config["ambiguous_palindromic_eaf_lower"])
        upper = float(config["ambiguous_palindromic_eaf_upper"])
        if is_palindromic(effect_allele, other_allele) and lower <= eaf <= upper:
            flags.append("ambiguous_palindromic")

    return {
        "variant_id": variant_id,
        "chromosome": chromosome,
        "position": position,
        "effect_allele": effect_allele,
        "other_allele": other_allele,
        "beta": beta,
        "se": se,
        "p": p_value,
        "n": n,
        "eaf": eaf,
        "info": info,
        "qc_flags": flags,
        "original_variant_id": original_variant_id,
    }


def apply_duplicate_flags(rows: list[dict[str, object]]) -> None:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["variant_id"]), []).append(row)
    for duplicates in grouped.values():
        if len(duplicates) == 1:
            continue
        ranked = sorted(duplicates, key=lambda item: (float(item["p"]), -abs(float(item["beta"] or 0.0))))
        for duplicate in ranked[1:]:
            duplicate["qc_flags"].append("duplicate_variant_id")


def apply_sample_size_flags(rows: list[dict[str, object]], config: dict) -> float | None:
    sample_sizes = [int(row["n"]) for row in rows if isinstance(row.get("n"), int)]
    if not sample_sizes:
        return None
    median_n = float(median(sample_sizes))
    low_cutoff = median_n * float(config["low_sample_size_ratio"])
    high_cutoff = median_n * float(config["high_sample_size_ratio"])
    for row in rows:
        n_value = row.get("n")
        if isinstance(n_value, int):
            if n_value < low_cutoff:
                row["qc_flags"].append("low_sample_size")
            elif n_value > high_cutoff:
                row["qc_flags"].append("high_sample_size")
    return median_n


def finalize_rows(rows: list[dict[str, object]]) -> None:
    exclusionary = {
        "invalid_p",
        "invalid_odds_ratio",
        "invalid_se",
        "invalid_n",
        "invalid_eaf",
        "invalid_info",
        "missing_effect_size",
        "low_info",
        "ambiguous_palindromic",
        "duplicate_variant_id",
    }
    for row in rows:
        unique_flags = sorted(set(str(flag) for flag in row["qc_flags"]))
        row["qc_flags"] = unique_flags
        include = not any(flag in exclusionary for flag in unique_flags)
        row["include_in_downstream"] = include
        row["qc_status"] = "pass" if include else "fail"


def summarize_rows(rows: list[dict[str, object]], config: dict, input_path: Path, out_tsv: Path, median_n: float | None) -> dict[str, object]:
    flag_counts = Counter(flag for row in rows for flag in row["qc_flags"])
    included = [row for row in rows if row["include_in_downstream"]]
    lead_variants = sorted(included, key=lambda item: (float(item["p"]), -abs(float(item["beta"] or 0.0))))[: int(config["lead_variant_limit"])]
    genome_wide = sum(1 for row in included if float(row["p"]) <= float(config["genome_wide_significance"]))
    suggestive = sum(1 for row in included if float(row["p"]) <= float(config["suggestive_significance"]))

    recommended_steps = [
        "Resolve failed rows before LDSC, meta-analysis refreshes, or PRS derivation.",
        "Clump retained lead variants with PLINK 2.0 before reporting independent loci.",
    ]
    recommended_resources = [
        {"resource_id": "plink2-docs", "reason": "Use for clumping, filtering, and downstream GWAS follow-up."},
        {"resource_id": "ldsc-repo", "reason": "Use after harmonization for SNP-heritability and genetic-correlation analyses."},
        {"resource_id": "fuma-docs", "reason": "Use for SNP2GENE and gene-set interpretation once lead loci are stable."},
    ]
    if any("ambiguous_palindromic" in row["qc_flags"] or "duplicate_variant_id" in row["qc_flags"] for row in rows):
        recommended_steps.insert(0, "Harmonize variants and headers with GWASLab before exporting to LDSC, FUMA, or GWAS-SSF.")
        recommended_resources.insert(
            0,
            {"resource_id": "gwaslab-docs", "reason": "Use for header normalization, harmonization, palindromic-SNP handling, and format conversion."},
        )
    if genome_wide:
        recommended_steps.append("Interpret the retained genome-wide significant loci with FUMA or a comparable locus-to-gene layer.")
    recommended_steps.append("If you need an exchange or submission format, export the cleaned results toward the GWAS Catalog summary-statistics standard.")
    recommended_resources.append(
        {
            "resource_id": "gwas-catalog-summary-statistics-docs",
            "reason": "Use for canonical summary-statistics field expectations and downstream interchange.",
        }
    )

    return {
        "input_path": str(input_path.resolve()),
        "flagged_tsv": str(out_tsv.resolve()),
        "row_count": len(rows),
        "qc_pass_count": len(included),
        "qc_fail_count": len(rows) - len(included),
        "median_sample_size": median_n,
        "flag_counts": dict(sorted(flag_counts.items())),
        "genome_wide_significant_count": genome_wide,
        "suggestive_hit_count": suggestive,
        "lead_variants": [
            {
                "variant_id": row["variant_id"],
                "chromosome": row["chromosome"],
                "position": row["position"],
                "p": row["p"],
                "beta": row["beta"],
                "qc_flags": row["qc_flags"],
            }
            for row in lead_variants
        ],
        "recommended_next_steps": recommended_steps,
        "recommended_resources": recommended_resources,
        "qc_protocol_anchor": "gwas-meta-analysis-qc-protocol",
    }


def write_flagged_tsv(rows: list[dict[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "variant_id": row["variant_id"],
                    "chromosome": row["chromosome"],
                    "position": row["position"],
                    "effect_allele": row["effect_allele"],
                    "other_allele": row["other_allele"],
                    "beta": "" if row["beta"] is None else row["beta"],
                    "se": "" if row["se"] is None else row["se"],
                    "p": row["p"],
                    "n": "" if row["n"] is None else row["n"],
                    "eaf": "" if row["eaf"] is None else row["eaf"],
                    "info": "" if row["info"] is None else row["info"],
                    "qc_status": row["qc_status"],
                    "include_in_downstream": str(bool(row["include_in_downstream"])).lower(),
                    "qc_flags": ",".join(row["qc_flags"]),
                }
            )


def write_json(payload: dict[str, object], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="GWAS summary-statistics TSV or CSV.")
    parser.add_argument("--config", type=Path, help="Optional JSON config with QC thresholds.")
    parser.add_argument("--out-tsv", type=Path, required=True, help="Flagged TSV output path.")
    parser.add_argument("--summary-out", type=Path, required=True, help="JSON summary output path.")
    args = parser.parse_args()

    config = load_config(args.config)
    raw_rows, fieldnames = read_rows(args.input)
    mapping = resolve_column_map(fieldnames)
    normalized = [normalize_variant(row, mapping, config) for row in raw_rows]
    apply_duplicate_flags(normalized)
    median_n = apply_sample_size_flags(normalized, config)
    finalize_rows(normalized)
    summary = summarize_rows(normalized, config, args.input, args.out_tsv, median_n)
    write_flagged_tsv(normalized, args.out_tsv)
    write_json(summary, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
