#!/usr/bin/env python3
"""Build a deterministic planning/execution graph over the local skill registry."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path
from typing import TypedDict


ROOT = Path(__file__).resolve().parents[4]
REGISTRY = ROOT / "registry"
ROUTER_SCRIPT = (
    ROOT
    / "skills"
    / "scientific-agents-and-automation"
    / "skill-registry-router-starter"
    / "scripts"
    / "route_skill_query.py"
)
TOKEN_RE = re.compile(r"[a-z0-9]+")


def load_langgraph():
    try:
        from langgraph.graph import END, START, StateGraph
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires the agents prefix at slurm/envs/agents.") from exc
    return StateGraph, START, END


def load_router_module():
    spec = importlib.util.spec_from_file_location("route_skill_query", ROUTER_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def skill_index() -> dict[str, dict]:
    return {skill["slug"]: skill for skill in load_jsonl(REGISTRY / "skills.jsonl")}


def tokenize(text: str) -> set[str]:
    return set(TOKEN_RE.findall(text.lower()))


def derive_step_queries(goal: str) -> list[str]:
    lowered = goal.lower()
    queries: list[str] = []
    if any(term in lowered for term in ["literature", "paper", "publication", "review"]):
        queries.append("literature search")
    if any(term in lowered for term in ["single-cell", "marker", "rank genes", "ranked genes", "scanpy"]):
        queries.append("single-cell marker ranking")
    if any(term in lowered for term in ["bayesian", "posterior", "uncertainty"]):
        queries.append("bayesian regression")
    if any(term in lowered for term in ["survival", "kaplan", "clinical outcome"]):
        queries.append("survival analysis")
    if any(term in lowered for term in ["pathway", "enrichment", "reactome"]):
        queries.append("reactome enrichment")
    if any(term in lowered for term in ["report", "dashboard", "interactive", "visual", "figure", "plot"]):
        queries.append("interactive html report" if "interactive" in lowered or "dashboard" in lowered or "report" in lowered else "publication figure")
    if not queries:
        queries.append(goal)
    deduped: list[str] = []
    for query in queries:
        if query not in deduped:
            deduped.append(query)
    return deduped


class AgentState(TypedDict):
    goal: str
    step_queries: list[str]
    routed_steps: list[dict]
    unique_skill_slugs: list[str]
    summary: str


def plan_node(state: AgentState) -> dict:
    return {"step_queries": derive_step_queries(state["goal"])}


def route_node(state: AgentState) -> dict:
    router = load_router_module()
    by_slug = skill_index()
    routed_steps = []
    unique_skill_slugs: list[str] = []
    for query in state["step_queries"]:
        result = router.route_query(query, top_k=2)
        if not result["matches"]:
            routed_steps.append(
                {
                    "step_query": query,
                    "selected_skill_slug": None,
                    "selected_skill_name": None,
                    "selected_skill_status": None,
                    "execution_preview": None,
                    "alternate_skill_slugs": [],
                    "rationale": ["no registry match"],
                }
            )
            continue
        primary = result["matches"][0]
        skill = by_slug[primary["slug"]]
        if primary["slug"] not in unique_skill_slugs:
            unique_skill_slugs.append(primary["slug"])
        routed_steps.append(
            {
                "step_query": query,
                "selected_skill_slug": primary["slug"],
                "selected_skill_name": primary["name"],
                "selected_skill_status": primary["status"],
                "execution_preview": (skill.get("test_commands") or [None])[0],
                "alternate_skill_slugs": [match["slug"] for match in result["matches"][1:]],
                "rationale": primary["rationale"],
            }
        )
    return {"routed_steps": routed_steps, "unique_skill_slugs": unique_skill_slugs}


def finalize_node(state: AgentState) -> dict:
    routed = state["routed_steps"]
    summary = (
        f"Planned {len(state['step_queries'])} step(s) from goal '{state['goal']}' and "
        f"selected {len(state['unique_skill_slugs'])} unique skill(s)."
    )
    return {"summary": summary}


def run_agent(goal: str) -> dict:
    StateGraph, START, END = load_langgraph()
    graph = StateGraph(AgentState)
    graph.add_node("plan", plan_node)
    graph.add_node("route", route_node)
    graph.add_node("finalize", finalize_node)
    graph.add_edge(START, "plan")
    graph.add_edge("plan", "route")
    graph.add_edge("route", "finalize")
    graph.add_edge("finalize", END)
    app = graph.compile()
    payload = app.invoke(
        {
            "goal": goal,
            "step_queries": [],
            "routed_steps": [],
            "unique_skill_slugs": [],
            "summary": "",
        }
    )
    matched = [step for step in payload["routed_steps"] if step["selected_skill_slug"]]
    recommended_execution_order = [step["selected_skill_slug"] for step in matched]
    return {
        "goal": goal,
        "plan_step_count": len(payload["step_queries"]),
        "matched_step_count": len(matched),
        "step_queries": payload["step_queries"],
        "routed_steps": payload["routed_steps"],
        "unique_skill_slugs": payload["unique_skill_slugs"],
        "recommended_execution_order": recommended_execution_order,
        "summary": payload["summary"],
    }


def write_json(payload: dict, out_path: Path | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path is None:
        print(text)
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--goal", required=True, help="Scientific task goal to plan.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()
    payload = run_agent(goal=args.goal)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
