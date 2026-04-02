#!/usr/bin/env python3
"""Route a scientific task query to candidate skills from the local registry."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
REGISTRY = ROOT / "registry"
TOKEN_RE = re.compile(r"[a-z0-9]+")


def normalize_token(token: str) -> str:
    if token.endswith("ing") and len(token) > 5:
        token = token[:-3]
    elif token.endswith("ed") and len(token) > 4:
        token = token[:-2]
    if token.endswith("s") and len(token) > 4 and not token.endswith(("ss", "is", "ics")):
        token = token[:-1]
    return token


def tokenize(text: str) -> list[str]:
    return [normalize_token(token) for token in TOKEN_RE.findall(text.lower())]


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def candidate_alias_phrases(skill: dict, aliases: dict[str, list[str]]) -> list[str]:
    phrases: list[str] = []
    skill_blob = " ".join(
        [
            skill["slug"],
            skill["name"].lower(),
            skill["domain"],
            " ".join(skill.get("tags", [])),
            " ".join(skill.get("topic_path", [])),
        ]
    ).lower()
    for key, values in aliases.items():
        key_blob = key.replace("_", " ").lower()
        if key_blob in skill_blob or key in skill_blob:
            phrases.extend(values)
    return sorted(set(phrases))


def route_query(query: str, top_k: int = 5) -> dict:
    skills = load_jsonl(REGISTRY / "skills.jsonl")
    aliases = load_json(REGISTRY / "aliases.json")
    query_tokens = set(tokenize(query))
    ranked = []
    lowered_query = query.lower()

    for skill in skills:
        alias_phrases = candidate_alias_phrases(skill, aliases)
        matched_aliases = [phrase for phrase in alias_phrases if phrase.lower() in lowered_query]
        topic_tokens = set(tokenize(" ".join(skill.get("topic_path", []))))
        tag_tokens = set(tokenize(" ".join(skill.get("tags", []))))
        name_tokens = set(tokenize(skill["name"]))
        domain_tokens = set(tokenize(skill["domain"]))
        overlap = sorted(query_tokens & (topic_tokens | tag_tokens | name_tokens | domain_tokens))
        lexical_score = 0
        lexical_score += 60 * len(matched_aliases)
        lexical_score += 12 * len(query_tokens & tag_tokens)
        lexical_score += 10 * len(query_tokens & name_tokens)
        lexical_score += 8 * len(query_tokens & topic_tokens)
        lexical_score += 4 * len(query_tokens & domain_tokens)
        if lexical_score == 0:
            continue
        score = lexical_score
        if skill["status"] == "slurm_verified":
            score += 6
        elif skill["status"] == "sandbox_verified":
            score += 4
        rationale = []
        if matched_aliases:
            rationale.append("matched aliases: " + ", ".join(matched_aliases[:3]))
        if overlap:
            rationale.append("token overlap: " + ", ".join(overlap[:6]))
        rationale.append(f"status: {skill['status']}")
        ranked.append(
            {
                "skill_id": skill["skill_id"],
                "name": skill["name"],
                "slug": skill["slug"],
                "domain": skill["domain"],
                "status": skill["status"],
                "score": score,
                "matched_aliases": matched_aliases,
                "matched_tokens": overlap,
                "rationale": rationale,
                "path": skill["path"],
            }
        )

    ranked.sort(key=lambda item: (-item["score"], item["name"]))
    return {"query": query, "top_k": top_k, "matches": ranked[:top_k], "match_count": len(ranked[:top_k])}


def write_payload(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", required=True, help="Free-text scientific task query.")
    parser.add_argument("--top-k", type=int, default=5, help="Maximum number of matches to return.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    if args.top_k <= 0:
        raise SystemExit("top-k must be a positive integer.")
    payload = route_query(query=args.query, top_k=args.top_k)
    write_payload(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
