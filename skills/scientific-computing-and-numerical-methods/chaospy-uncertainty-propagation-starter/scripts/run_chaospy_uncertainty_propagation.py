#!/usr/bin/env python3
"""Run a deterministic Chaospy uncertainty-propagation toy workflow."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path


DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "examples" / "toy_parameters.json"


def load_chaospy():
    try:
        import chaospy as cp
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit("This script requires Chaospy in slurm/envs/numerics.") from exc
    return cp


def load_config(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "normal_mean",
        "normal_std",
        "uniform_low",
        "uniform_high",
        "polynomial_order",
        "quadrature_order",
    }
    missing = sorted(required - payload.keys())
    if missing:
        raise ValueError(f"Config missing required keys: {', '.join(missing)}")
    if payload["normal_std"] <= 0:
        raise ValueError("normal_std must be positive.")
    if payload["uniform_low"] >= payload["uniform_high"]:
        raise ValueError("uniform_low must be less than uniform_high.")
    if payload["polynomial_order"] <= 0 or payload["quadrature_order"] <= 0:
        raise ValueError("polynomial_order and quadrature_order must be positive integers.")
    return payload


def response_model(stiffness, damping):
    return stiffness**2 + 0.5 * damping


def summarize_uncertainty(config: dict) -> dict:
    cp = load_chaospy()
    distribution = cp.J(
        cp.Normal(config["normal_mean"], config["normal_std"]),
        cp.Uniform(config["uniform_low"], config["uniform_high"]),
    )
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*where' used without 'out'.*",
            category=UserWarning,
        )
        basis = cp.expansion.stieltjes(config["polynomial_order"], distribution)
        nodes, weights = cp.generate_quadrature(config["quadrature_order"], distribution, rule="gaussian")
        values = response_model(nodes[0], nodes[1])
        surrogate = cp.fit_quadrature(basis, nodes, weights, values)
        mean = float(cp.E(surrogate, distribution))
        std = float(cp.Std(surrogate, distribution))
        p05, p50, p95 = [float(value) for value in cp.Perc(surrogate, [5, 50, 95], distribution)]

    evaluations = []
    for stiffness, damping, response in zip(nodes[0].tolist(), nodes[1].tolist(), values.tolist(), strict=True):
        evaluations.append(
            {
                "stiffness": round(float(stiffness), 6),
                "damping": round(float(damping), 6),
                "response": round(float(response), 6),
            }
        )

    return {
        "model": "stiffness^2 + 0.5*damping",
        "distribution": {
            "stiffness": {
                "family": "Normal",
                "mean": round(float(config["normal_mean"]), 6),
                "std": round(float(config["normal_std"]), 6),
            },
            "damping": {
                "family": "Uniform",
                "low": round(float(config["uniform_low"]), 6),
                "high": round(float(config["uniform_high"]), 6),
            },
        },
        "polynomial_order": int(config["polynomial_order"]),
        "quadrature_order": int(config["quadrature_order"]),
        "node_count": int(len(weights)),
        "mean": round(mean, 6),
        "std": round(std, 6),
        "percentiles": {
            "p05": round(p05, 6),
            "p50": round(p50, 6),
            "p95": round(p95, 6),
        },
        "evaluations": evaluations,
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
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="JSON configuration path.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    try:
        config = load_config(args.config)
        payload = summarize_uncertainty(config)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
