#!/usr/bin/env python3
"""Integrate a deterministic Lotka-Volterra ODE system with SciPy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def lotka_volterra(_: float, state: np.ndarray) -> list[float]:
    prey, predator = state
    alpha = 1.1
    beta = 0.4
    delta = 0.1
    gamma = 0.4
    return [
        alpha * prey - beta * prey * predator,
        delta * prey * predator - gamma * predator,
    ]


def run_simulation(duration: float, points: int) -> dict[str, object]:
    t_eval = np.linspace(0.0, duration, points)
    solution = solve_ivp(
        lotka_volterra,
        t_span=(0.0, duration),
        y0=(10.0, 5.0),
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-9,
    )
    if not solution.success:  # pragma: no cover - guarded by solve_ivp runtime
        raise SystemExit(f"ODE solve failed: {solution.message}")

    prey = solution.y[0]
    predator = solution.y[1]
    prey_peak_idx = int(np.argmax(prey))
    predator_peak_idx = int(np.argmax(predator))
    return {
        "model": "lotka_volterra",
        "duration": round(float(duration), 6),
        "point_count": int(points),
        "initial_state": {"prey": 10.0, "predator": 5.0},
        "final_state": {
            "prey": round(float(prey[-1]), 6),
            "predator": round(float(predator[-1]), 6),
        },
        "prey_peak_time": round(float(solution.t[prey_peak_idx]), 6),
        "predator_peak_time": round(float(solution.t[predator_peak_idx]), 6),
        "prey_max": round(float(prey.max()), 6),
        "predator_max": round(float(predator.max()), 6),
        "coexistence_floor": round(float(min(prey.min(), predator.min())), 6),
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
    parser.add_argument("--duration", type=float, default=15.0, help="Simulation duration.")
    parser.add_argument("--points", type=int, default=151, help="Number of output points.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.duration <= 0:
        raise SystemExit("--duration must be positive.")
    if args.points < 20:
        raise SystemExit("--points must be at least 20.")

    payload = run_simulation(duration=args.duration, points=args.points)
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
