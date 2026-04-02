#!/usr/bin/env python3
"""Render and simulate a minimal Opentrons liquid-handling protocol."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from opentrons import simulate


PROTOCOL_TEMPLATE = """metadata = {{"protocolName": "Toy Transfer", "apiLevel": "{api_level}"}}

def run(protocol):
    plate = protocol.load_labware("{labware}", 1)
    tiprack = protocol.load_labware("{tiprack}", 2)
    pipette = protocol.load_instrument("{pipette}", mount="{mount}", tip_racks=[tiprack])
    pipette.transfer({volume}, plate["{source}"], plate["{destination}"], new_tip="once")
"""


def render_protocol(
    api_level: str,
    labware: str,
    tiprack: str,
    pipette: str,
    mount: str,
    source: str,
    destination: str,
    volume: float,
) -> str:
    return PROTOCOL_TEMPLATE.format(
        api_level=api_level,
        labware=labware,
        tiprack=tiprack,
        pipette=pipette,
        mount=mount,
        source=source,
        destination=destination,
        volume=volume,
    )


def simulate_protocol(protocol_path: Path) -> dict[str, object]:
    with protocol_path.open("r", encoding="utf-8") as handle:
        commands, _bundle = simulate.simulate(handle, file_name=str(protocol_path))
    command_text = [entry["payload"]["text"] for entry in commands]
    return {
        "command_count": len(command_text),
        "first_command": command_text[0] if command_text else None,
        "last_command": command_text[-1] if command_text else None,
        "commands": command_text,
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
    parser.add_argument("--api-level", default="2.15", help="Protocol API level.")
    parser.add_argument("--labware", default="corning_96_wellplate_360ul_flat", help="Source/destination labware.")
    parser.add_argument("--tiprack", default="opentrons_96_tiprack_300ul", help="Tiprack definition.")
    parser.add_argument("--pipette", default="p300_single_gen2", help="Instrument name.")
    parser.add_argument("--mount", default="left", choices=["left", "right"], help="Mount side.")
    parser.add_argument("--source", default="A1", help="Source well.")
    parser.add_argument("--destination", default="B1", help="Destination well.")
    parser.add_argument("--volume", type=float, default=50.0, help="Transfer volume in microliters.")
    parser.add_argument("--protocol-out", type=Path, required=True, help="Path to write the rendered protocol script.")
    parser.add_argument("--out", type=Path, default=None, help="Optional JSON summary path.")
    args = parser.parse_args()

    if args.volume <= 0:
        raise SystemExit("--volume must be positive.")

    protocol_text = render_protocol(
        api_level=args.api_level,
        labware=args.labware,
        tiprack=args.tiprack,
        pipette=args.pipette,
        mount=args.mount,
        source=args.source,
        destination=args.destination,
        volume=args.volume,
    )
    args.protocol_out.parent.mkdir(parents=True, exist_ok=True)
    args.protocol_out.write_text(protocol_text, encoding="utf-8")
    simulation = simulate_protocol(args.protocol_out)
    payload = {
        "api_level": args.api_level,
        "labware": args.labware,
        "tiprack": args.tiprack,
        "pipette": args.pipette,
        "mount": args.mount,
        "source": args.source,
        "destination": args.destination,
        "volume_ul": round(float(args.volume), 6),
        "protocol_path": str(args.protocol_out),
        **simulation,
    }
    write_json(payload, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
