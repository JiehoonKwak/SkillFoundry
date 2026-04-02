#!/usr/bin/env python3
"""Optimize a small synthetic objective with Optuna and summarize the best trial."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import optuna


def objective(trial: optuna.Trial) -> float:
    x = trial.suggest_float("x", -2.0, 2.0)
    y = trial.suggest_float("y", -2.0, 2.0)
    return (x - 0.5) ** 2 + (y + 0.25) ** 2


def optimize(n_trials: int, seed: int) -> dict:
    if n_trials < 4:
        raise ValueError("n_trials must be at least 4 to exercise both startup and model-guided trials.")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=seed),
    )
    study.optimize(objective, n_trials=n_trials, show_progress_bar=False)
    ranked_trials = sorted(
        [
            {
                "number": trial.number,
                "value": round(float(trial.value), 6),
                "params": {key: round(float(value), 6) for key, value in trial.params.items()},
            }
            for trial in study.trials
            if trial.value is not None
        ],
        key=lambda item: item["value"],
    )
    best = ranked_trials[0]
    return {
        "study_direction": "minimize",
        "sampler": "TPESampler",
        "seed": int(seed),
        "trial_count": int(n_trials),
        "best_trial_number": int(best["number"]),
        "best_value": best["value"],
        "best_params": best["params"],
        "top_trials": ranked_trials[:5],
        "search_space": {
            "x": [-2.0, 2.0],
            "y": [-2.0, 2.0],
        },
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
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    parser.add_argument("--trials", type=int, default=32, help="Number of optimization trials.")
    parser.add_argument("--seed", type=int, default=7, help="Deterministic sampler seed.")
    args = parser.parse_args()

    try:
        payload = optimize(args.trials, args.seed)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
