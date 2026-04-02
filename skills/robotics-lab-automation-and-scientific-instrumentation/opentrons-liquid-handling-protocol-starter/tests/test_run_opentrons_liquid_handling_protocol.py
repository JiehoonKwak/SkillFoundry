from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
PYTHON = ROOT / "slurm" / "envs" / "automation" / "bin" / "python"
SCRIPT = (
    ROOT
    / "skills"
    / "robotics-lab-automation-and-scientific-instrumentation"
    / "opentrons-liquid-handling-protocol-starter"
    / "scripts"
    / "run_opentrons_liquid_handling_protocol.py"
)


class OpentronsLiquidHandlingProtocolTests(unittest.TestCase):
    def test_renders_and_simulates_one_transfer_protocol(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            protocol_path = Path(tmp_dir) / "protocol.py"
            out_path = Path(tmp_dir) / "summary.json"
            subprocess.run(
                [
                    str(PYTHON),
                    str(SCRIPT),
                    "--protocol-out",
                    str(protocol_path),
                    "--out",
                    str(out_path),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
                timeout=240,
            )
            payload = json.loads(out_path.read_text(encoding="utf-8"))
            protocol_text = protocol_path.read_text(encoding="utf-8")
            self.assertIn("def run(protocol):", protocol_text)
            self.assertIn("pipette.transfer", protocol_text)
            self.assertEqual(payload["command_count"], 5)
            self.assertIn("Transferring 50.0", payload["first_command"])
            self.assertIn("Dropping tip", payload["last_command"])

    def test_rejects_nonpositive_volume(self) -> None:
        completed = subprocess.run(
            [str(PYTHON), str(SCRIPT), "--protocol-out", "ignored.py", "--volume", "0"],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=240,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("--volume must be positive", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
