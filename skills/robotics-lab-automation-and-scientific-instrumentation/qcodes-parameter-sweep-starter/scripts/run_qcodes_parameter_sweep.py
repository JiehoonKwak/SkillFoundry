#!/usr/bin/env python3
"""Run a deterministic local QCoDeS parameter sweep without hardware."""

from __future__ import annotations

import argparse
import csv
import io
import json
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from qcodes.dataset import Measurement, initialise_or_create_database_at, load_by_run_spec, load_or_create_experiment
from qcodes.parameters import ManualParameter


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SETPOINTS = ROOT / "examples" / "toy_setpoints.tsv"


def read_setpoints(path: Path) -> list[float]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "setpoint_v" not in (reader.fieldnames or []):
            raise ValueError("Expected a 'setpoint_v' column in the setpoints TSV.")
        values = [float(row["setpoint_v"]) for row in reader]
    if not values:
        raise ValueError("Setpoints TSV is empty.")
    return values


def response_model(value: float) -> float:
    return round(1.2 * value * value + 0.1, 6)


def build_summary(setpoints: list[float], db_out: Path, experiment_name: str, sample_name: str) -> dict[str, object]:
    if db_out.exists():
        db_out.unlink()
    db_out.parent.mkdir(parents=True, exist_ok=True)
    initialise_or_create_database_at(str(db_out))
    experiment = load_or_create_experiment(experiment_name, sample_name=sample_name)
    drive = ManualParameter("gate_voltage", unit="V", initial_value=setpoints[0])
    response = ManualParameter("signal", unit="nA", initial_value=response_model(setpoints[0]))
    measurement = Measurement(exp=experiment)
    measurement.register_parameter(drive)
    measurement.register_parameter(response, setpoints=[drive])

    buffer = io.StringIO()
    with redirect_stdout(buffer), redirect_stderr(buffer):
        with measurement.run() as datasaver:
            captured_points: list[dict[str, float]] = []
            for value in setpoints:
                signal = response_model(value)
                drive(value)
                response(signal)
                datasaver.add_result((drive, value), (response, signal))
                captured_points.append({"setpoint_v": round(value, 6), "response_na": signal})
            run_id = datasaver.run_id

    dataset = load_by_run_spec(captured_run_id=run_id)
    parameter_data = dataset.get_parameter_data()["signal"]
    responses = [round(float(value), 6) for value in parameter_data["signal"].tolist()]
    stored_setpoints = [round(float(value), 6) for value in parameter_data["gate_voltage"].tolist()]
    return {
        "db_path": str(db_out),
        "run_id": int(run_id),
        "experiment_name": experiment_name,
        "sample_name": sample_name,
        "point_count": len(stored_setpoints),
        "setpoint_min": min(stored_setpoints),
        "setpoint_max": max(stored_setpoints),
        "response_mean": round(sum(responses) / len(responses), 6),
        "response_max": max(responses),
        "captured_points": [
            {"setpoint_v": setpoint, "response_na": response_value}
            for setpoint, response_value in zip(stored_setpoints, responses)
        ],
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
    parser.add_argument("--setpoints", type=Path, default=DEFAULT_SETPOINTS, help="TSV file with a setpoint_v column.")
    parser.add_argument("--db-out", type=Path, required=True, help="SQLite database output path.")
    parser.add_argument("--experiment-name", default="toy_parameter_sweep", help="Experiment name.")
    parser.add_argument("--sample-name", default="synthetic_demo", help="Sample name.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    args = parser.parse_args()

    if not args.setpoints.exists():
        raise SystemExit(f"Setpoints TSV not found: {args.setpoints}")

    setpoints = read_setpoints(args.setpoints)
    payload = build_summary(setpoints, args.db_out, args.experiment_name, args.sample_name)
    write_json(payload, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
