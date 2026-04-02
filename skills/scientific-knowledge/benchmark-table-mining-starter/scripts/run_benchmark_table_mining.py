#!/usr/bin/env python3
"""Parse markdown benchmark notes and summarize the strongest leaderboard table."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "examples" / "benchmark_note.md"
FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
SCORE_LIKE_EXACT = {
    "accuracy",
    "ari",
    "auc",
    "auprc",
    "auroc",
    "bleu",
    "f1",
    "macro f1",
    "mae",
    "map",
    "micro f1",
    "mrr",
    "nmi",
    "pearson",
    "recall",
    "rmse",
    "rouge",
    "score",
    "spearman",
}
SCORE_LIKE_TOKENS = (
    "acc",
    "accuracy",
    "ari",
    "auc",
    "auprc",
    "auroc",
    "bleu",
    "f1",
    "map",
    "mrr",
    "nmi",
    "pearson",
    "precision",
    "recall",
    "rmse",
    "rouge",
    "score",
    "spearman",
)
LOWER_IS_BETTER_TOKENS = ("cer", "error", "loss", "mae", "mse", "perplexity", "rmse", "wer")


def split_table_blocks(lines: list[str]) -> list[list[str]]:
    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("|"):
            current.append(stripped)
            continue
        if current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    return blocks


def parse_table_block(block: list[str]) -> list[dict[str, str]]:
    if len(block) < 3:
        return []
    headers = [cell.strip() for cell in block[0].strip("|").split("|")]
    rows: list[dict[str, str]] = []
    for line in block[2:]:
        values = [cell.strip() for cell in line.strip("|").split("|")]
        if len(values) != len(headers):
            continue
        rows.append(dict(zip(headers, values)))
    return rows


def parse_tables(lines: list[str]) -> list[list[dict[str, str]]]:
    tables = []
    for block in split_table_blocks(lines):
        rows = parse_table_block(block)
        if rows:
            tables.append(rows)
    return tables


def normalize_header(header: str) -> str:
    return " ".join(header.strip().lower().split())


def normalize_key(header: str) -> str:
    normalized = normalize_header(header)
    key = re.sub(r"[^a-z0-9]+", "_", normalized).strip("_")
    return key or "column"


def to_float(value: str) -> float | None:
    cleaned = value.strip().replace(",", "")
    if cleaned.endswith("%"):
        cleaned = cleaned[:-1]
    match = FLOAT_RE.search(cleaned)
    if match is None:
        return None
    try:
        return float(match.group())
    except ValueError:
        return None


def header_lookup(headers: list[str]) -> list[tuple[str, str]]:
    return [(normalize_header(header), header) for header in headers]


def choose_header(headers: list[str], exact: tuple[str, ...], contains: tuple[str, ...] = ()) -> str | None:
    normalized_headers = header_lookup(headers)
    for target in exact:
        for normalized, header in normalized_headers:
            if normalized == target:
                return header
    for target in contains:
        for normalized, header in normalized_headers:
            if target in normalized:
                return header
    return None


def choose_name_column(headers: list[str]) -> str | None:
    return choose_header(headers, ("method", "model", "approach"), ("method", "model", "approach", "system"))


def choose_dataset_column(headers: list[str]) -> str | None:
    return choose_header(headers, ("dataset", "benchmark", "corpus"), ("dataset", "benchmark", "corpus"))


def choose_task_column(headers: list[str]) -> str | None:
    return choose_header(headers, ("task", "subtask", "problem"), ("task", "subtask", "problem"))


def choose_metric_column(headers: list[str]) -> str | None:
    return choose_header(headers, ("metric", "measure"), ("metric", "measure"))


def is_score_like(header: str) -> bool:
    normalized = normalize_header(header)
    return normalized in SCORE_LIKE_EXACT or any(token in normalized for token in SCORE_LIKE_TOKENS)


def choose_score_column(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return None
    headers = list(rows[0])
    candidates: list[tuple[int, int, int, int, str]] = []
    for index, header in enumerate(headers):
        numeric_count = sum(to_float(row.get(header, "")) is not None for row in rows)
        if numeric_count == 0:
            continue
        candidates.append((1 if is_score_like(header) else 0, numeric_count, -index, len(header), header))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    return candidates[0][4]


def infer_sort_direction(score_header: str, metric_header: str | None, rows: list[dict[str, str]]) -> str:
    labels = [normalize_header(score_header)]
    if metric_header is not None:
        labels.extend(normalize_header(row.get(metric_header, "")) for row in rows if row.get(metric_header))
    joined = " ".join(labels)
    if any(token in joined for token in LOWER_IS_BETTER_TOKENS):
        return "ascending"
    return "descending"


def describe_table(index: int, rows: list[dict[str, str]]) -> dict[str, object]:
    headers = list(rows[0]) if rows else []
    score_header = choose_score_column(rows)
    score_numeric_count = 0
    if score_header is not None:
        score_numeric_count = sum(to_float(row.get(score_header, "")) is not None for row in rows)
    name_header = choose_name_column(headers)
    dataset_header = choose_dataset_column(headers)
    task_header = choose_task_column(headers)
    metric_header = choose_metric_column(headers)
    leaderboard_score = 0
    if score_header is not None:
        leaderboard_score += 100
        leaderboard_score += score_numeric_count * 2
        leaderboard_score += len(rows)
        if is_score_like(score_header):
            leaderboard_score += 30
    if name_header is not None:
        leaderboard_score += 20
    if dataset_header is not None:
        leaderboard_score += 10
    if task_header is not None:
        leaderboard_score += 10
    if metric_header is not None:
        leaderboard_score += 10
    return {
        "index": index,
        "headers": headers,
        "rows": rows,
        "score_header": score_header,
        "score_numeric_count": score_numeric_count,
        "name_header": name_header,
        "dataset_header": dataset_header,
        "task_header": task_header,
        "metric_header": metric_header,
        "leaderboard_score": leaderboard_score,
        "is_score_like": bool(score_header and is_score_like(score_header)),
    }


def select_leaderboard_table(tables: list[list[dict[str, str]]]) -> dict[str, object]:
    candidates = [describe_table(index, rows) for index, rows in enumerate(tables)]
    candidates = [candidate for candidate in candidates if candidate["score_header"] is not None]
    if not candidates:
        raise SystemExit("Could not identify a numeric score column in any markdown table.")
    candidates.sort(
        key=lambda candidate: (
            int(candidate["leaderboard_score"]),
            1 if candidate["is_score_like"] else 0,
            int(candidate["score_numeric_count"]),
            len(candidate["rows"]),
            -int(candidate["index"]),
        ),
        reverse=True,
    )
    return candidates[0]


def build_ranked_rows(table: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]], str]:
    rows = table["rows"]
    if not isinstance(rows, list):
        raise SystemExit("Selected table rows were not available.")
    score_header = table["score_header"]
    if not isinstance(score_header, str):
        raise SystemExit("Selected table is missing a score column.")
    name_header = table.get("name_header")
    dataset_header = table.get("dataset_header")
    task_header = table.get("task_header")
    metric_header = table.get("metric_header")
    sort_direction = infer_sort_direction(score_header, metric_header if isinstance(metric_header, str) else None, rows)

    ranked_rows: list[dict[str, object]] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        score_text = row.get(score_header, "")
        score_value = to_float(score_text)
        if score_value is None:
            continue
        normalized_row = {normalize_key(header): value for header, value in row.items()}
        metric_value = row.get(metric_header, score_header) if isinstance(metric_header, str) else score_header
        name_value = row.get(name_header, "") if isinstance(name_header, str) else ""
        ranked_row: dict[str, object] = {
            **normalized_row,
            "metric_name": score_header,
            "metric": metric_value,
            "metric_value": score_value,
            "score": score_text,
            "score_value": score_value,
        }
        if isinstance(task_header, str):
            ranked_row["task"] = row.get(task_header, "")
        if isinstance(dataset_header, str):
            ranked_row["dataset"] = row.get(dataset_header, "")
        if name_value:
            ranked_row["model"] = name_value
            if normalize_header(name_header) in {"approach", "method"}:
                ranked_row["method"] = name_value
        ranked_rows.append(ranked_row)

    ranked_rows.sort(
        key=lambda item: (
            float(item["metric_value"]) if sort_direction == "ascending" else -float(item["metric_value"]),
            str(item.get("method") or item.get("model") or ""),
        )
    )

    rows_payload: list[dict[str, object]] = []
    for rank, ranked_row in enumerate(ranked_rows, start=1):
        ranked_row["rank"] = rank
        canonical_row: dict[str, object] = {"rank": rank, "score": ranked_row["score"]}
        if ranked_row.get("task"):
            canonical_row["task"] = ranked_row["task"]
        if ranked_row.get("dataset"):
            canonical_row["dataset"] = ranked_row["dataset"]
        if ranked_row.get("metric"):
            canonical_row["metric"] = ranked_row["metric"]
        if ranked_row.get("model"):
            canonical_row["model"] = ranked_row["model"]
        elif ranked_row.get("method"):
            canonical_row["model"] = ranked_row["method"]
        rows_payload.append(canonical_row)
    return ranked_rows, rows_payload, sort_direction


def build_summary(input_path: Path) -> dict[str, object]:
    tables = parse_tables(input_path.read_text(encoding="utf-8").splitlines())
    if not tables:
        raise SystemExit("Expected at least one markdown table with header, separator, and data rows.")
    selected_table = select_leaderboard_table(tables)
    ranked_rows, rows_payload, sort_direction = build_ranked_rows(selected_table)
    score_header = selected_table["score_header"]
    headers = selected_table["headers"]
    if not isinstance(score_header, str) or not isinstance(headers, list):
        raise SystemExit("Selected table metadata was incomplete.")
    return {
        "input_path": str(input_path),
        "table_count": len(tables),
        "selected_table_index": selected_table["index"],
        "headers": headers,
        "columns": [normalize_key(header) for header in headers if isinstance(header, str)],
        "row_count": len(selected_table["rows"]),
        "metric_name": score_header,
        "score_column": score_header,
        "sort_direction": sort_direction,
        "best_method": ranked_rows[0] if ranked_rows else None,
        "best_row": ranked_rows[0] if ranked_rows else None,
        "ranked_rows": ranked_rows,
        "rows": rows_payload,
    }


def write_json(payload: dict[str, object], out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    if not args.input.exists():
        raise SystemExit(f"Input note not found: {args.input}")
    write_json(build_summary(args.input), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
