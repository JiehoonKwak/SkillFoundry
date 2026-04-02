#!/usr/bin/env python3
"""Run a deterministic Inspect evaluation harness over toy scientific-agent cases."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from inspect_ai import Task, eval
from inspect_ai.dataset import Sample
from inspect_ai.model import ModelOutput
from inspect_ai.scorer import match
from inspect_ai.solver import Generate, TaskState, solver


DEFAULT_CASES = Path(__file__).resolve().parents[1] / "examples" / "toy_eval_cases.json"


def load_cases(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(f"Cases file not found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list) or not payload:
        raise ValueError("Cases file must contain a non-empty JSON list.")
    required = {"id", "prompt", "target", "candidate", "baseline"}
    for index, row in enumerate(payload):
        missing = sorted(required - row.keys())
        if missing:
            raise ValueError(f"Case {index} missing required keys: {', '.join(missing)}")
    return payload


@solver
def metadata_answer_solver(answer_key: str):
    async def solve(state: TaskState, generate: Generate) -> TaskState:
        answer = str(state.metadata[answer_key])
        state.output = ModelOutput.from_content(answer_key, answer)
        state.completed = True
        return state

    return solve


def build_task(cases: list[dict], answer_key: str) -> Task:
    dataset = [
        Sample(
            input=case["prompt"],
            target=case["target"],
            id=case["id"],
            metadata={
                "candidate": case["candidate"],
                "baseline": case["baseline"],
            },
        )
        for case in cases
    ]
    return Task(
        dataset=dataset,
        solver=metadata_answer_solver(answer_key),
        scorer=match(location="exact"),
        model=None,
        name=f"{answer_key}_scientific_agent_eval",
    )


def run_eval(task: Task, log_dir: Path) -> dict:
    log_dir.mkdir(parents=True, exist_ok=True)
    for path in log_dir.glob("*.json"):
        path.unlink()
    logs = eval(task, model=None, log_dir=str(log_dir), log_format="json")
    log = logs[0]
    accuracy = float(log.results.scores[0].metrics["accuracy"].value)
    samples = []
    for sample in log.samples or []:
        score = sample.scores["match"]
        target = sample.target.text if hasattr(sample.target, "text") else str(sample.target)
        samples.append(
            {
                "id": sample.id,
                "output": sample.output.completion,
                "score": score.value,
                "target": target,
            }
        )
    log_files = sorted(path.name for path in log_dir.glob("*.json"))
    return {
        "accuracy": round(accuracy, 6),
        "sample_count": len(samples),
        "log_dir": str(log_dir),
        "log_file_count": len(log_files),
        "samples": samples,
    }


def summarize(cases: list[dict], log_root: Path) -> dict:
    candidate = run_eval(build_task(cases, "candidate"), log_root / "candidate")
    baseline = run_eval(build_task(cases, "baseline"), log_root / "baseline")
    return {
        "case_count": len(cases),
        "candidate_accuracy": candidate["accuracy"],
        "baseline_accuracy": baseline["accuracy"],
        "accuracy_gain": round(candidate["accuracy"] - baseline["accuracy"], 6),
        "candidate": candidate,
        "baseline": baseline,
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
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES, help="JSON cases file.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--log-dir", type=Path, required=True, help="Directory for Inspect JSON logs.")
    args = parser.parse_args()

    try:
        cases = load_cases(args.cases)
        payload = summarize(cases, args.log_dir)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
