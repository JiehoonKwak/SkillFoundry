#!/usr/bin/env python3
"""Estimate a deterministic toy average treatment effect with DoWhy."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path


def load_stack():
    try:
        import numpy as np
        import pandas as pd
        from dowhy import CausalModel
    except ImportError as exc:  # pragma: no cover - exercised via runtime tests
        raise SystemExit(
            "This script requires the causal prefix at slurm/envs/causal."
        ) from exc
    return np, pd, CausalModel


def build_toy_dataset(sample_count: int, seed: int):
    np, pd, _ = load_stack()
    rng = np.random.default_rng(seed)
    z = np.linspace(-1.0, 1.0, sample_count)
    treatment_noise = rng.normal(0.0, 0.08, sample_count)
    outcome_noise = rng.normal(0.0, 0.05, sample_count)
    treatment = (0.7 * z + treatment_noise > 0).astype(int)
    outcome = 1.9 * treatment + 0.8 * z + outcome_noise
    return pd.DataFrame({"z": z, "treatment": treatment, "outcome": outcome})


def estimate_average_treatment_effect(sample_count: int, seed: int, placebo_seed: int) -> dict:
    _, _, CausalModel = load_stack()
    data = build_toy_dataset(sample_count=sample_count, seed=seed)

    warnings.filterwarnings("ignore", category=FutureWarning)
    model = CausalModel(
        data=data,
        treatment="treatment",
        outcome="outcome",
        common_causes=["z"],
    )
    identified = model.identify_effect(proceed_when_unidentifiable=False)
    estimate = model.estimate_effect(identified, method_name="backdoor.linear_regression")
    refutation = model.refute_estimate(
        identified,
        estimate,
        method_name="placebo_treatment_refuter",
        placebo_type="permute",
        random_state=placebo_seed,
    )

    treatment_rate = float(data["treatment"].mean())
    confounder_treatment_correlation = float(data["z"].corr(data["treatment"]))
    observed_outcome_gap = float(
        data.loc[data["treatment"] == 1, "outcome"].mean()
        - data.loc[data["treatment"] == 0, "outcome"].mean()
    )

    return {
        "sample_count": int(sample_count),
        "seed": int(seed),
        "placebo_seed": int(placebo_seed),
        "estimand_type": "nonparametric_ate",
        "estimation_method": "backdoor.linear_regression",
        "refutation_method": "placebo_treatment_refuter",
        "treatment_rate": round(treatment_rate, 6),
        "confounder_treatment_correlation": round(confounder_treatment_correlation, 6),
        "observed_outcome_gap": round(observed_outcome_gap, 6),
        "average_treatment_effect": round(float(estimate.value), 6),
        "placebo_effect": round(float(refutation.new_effect), 6),
        "target_effect": 1.9,
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
    parser.add_argument("--sample-count", type=int, default=40, help="Number of synthetic samples.")
    parser.add_argument("--seed", type=int, default=7, help="Random seed for the toy dataset.")
    parser.add_argument("--placebo-seed", type=int, default=13, help="Random seed for the placebo refuter.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON output path.")
    args = parser.parse_args()

    if args.sample_count < 10:
        raise SystemExit("sample-count must be at least 10.")

    payload = estimate_average_treatment_effect(
        sample_count=args.sample_count,
        seed=args.seed,
        placebo_seed=args.placebo_seed,
    )
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
