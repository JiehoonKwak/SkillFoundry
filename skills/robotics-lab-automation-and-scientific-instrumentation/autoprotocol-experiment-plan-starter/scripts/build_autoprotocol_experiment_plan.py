#!/usr/bin/env python3
"""Build a deterministic Autoprotocol liquid-handling plan from a tiny transfer table."""

from __future__ import annotations

import argparse
import csv
import json
import warnings
from pathlib import Path


DEFAULT_TRANSFERS = Path(__file__).resolve().parents[1] / "examples" / "toy_transfers.tsv"


def load_transfers(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Transfer plan not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"source_well", "dest_well", "volume_ul"}
    if not rows:
        raise ValueError("Transfer plan is empty.")
    if missing := sorted(required - set(rows[0].keys())):
        raise ValueError(f"Transfer plan missing required columns: {', '.join(missing)}")
    for row in rows:
        volume = float(row["volume_ul"])
        if volume <= 0:
            raise ValueError("All transfer volumes must be positive.")
    return rows


def build_protocol(transfers: list[dict[str, str]], incubation_minutes: int, temperature: str):
    warnings.filterwarnings("ignore", message=r"pkg_resources is deprecated.*", category=UserWarning)
    from autoprotocol.protocol import Protocol

    protocol = Protocol()
    source = protocol.ref("source_plate", id=None, cont_type="96-pcr", storage="cold_4")
    destination = protocol.ref("assay_plate", id=None, cont_type="96-pcr", discard=True)
    for row in transfers:
        protocol.transfer(
            source.well(row["source_well"]),
            destination.well(row["dest_well"]),
            f"{float(row['volume_ul']):.1f}:microliter",
        )
    protocol.incubate(destination, temperature, f"{incubation_minutes}:minute")
    return protocol


def summarize_protocol(protocol_dict: dict, transfers: list[dict[str, str]], protocol_path: Path, incubation_minutes: int) -> dict:
    instructions = protocol_dict["instructions"]
    ops = [instruction["op"] for instruction in instructions]
    return {
        "ref_names": sorted(protocol_dict["refs"].keys()),
        "ref_count": len(protocol_dict["refs"]),
        "instruction_ops": ops,
        "instruction_count": len(instructions),
        "transfer_count": len(transfers),
        "source_wells": [row["source_well"] for row in transfers],
        "dest_wells": [row["dest_well"] for row in transfers],
        "incubation_minutes": int(incubation_minutes),
        "sealed_before_incubate": "seal" in ops,
        "protocol_json_path": str(protocol_path),
    }


def write_json(payload: dict, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transfers", type=Path, default=DEFAULT_TRANSFERS, help="TSV transfer plan.")
    parser.add_argument("--protocol-json", type=Path, required=True, help="Output Autoprotocol JSON path.")
    parser.add_argument("--summary-out", type=Path, default=None, help="Optional JSON summary output path.")
    parser.add_argument("--incubation-minutes", type=int, default=30, help="Incubation duration in minutes.")
    parser.add_argument("--temperature", default="warm_37", help="Autoprotocol incubation location token.")
    args = parser.parse_args()

    if args.incubation_minutes <= 0:
        raise SystemExit("incubation-minutes must be positive.")

    try:
        transfers = load_transfers(args.transfers)
        protocol = build_protocol(transfers, args.incubation_minutes, args.temperature)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    protocol_dict = protocol.as_dict()
    write_json(protocol_dict, args.protocol_json)
    summary = summarize_protocol(protocol_dict, transfers, args.protocol_json, args.incubation_minutes)
    if args.summary_out is None:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        write_json(summary, args.summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
