from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class RepositoryValidationTests(unittest.TestCase):
    def test_repository_validator_passes(self) -> None:
        completed = subprocess.run(
            ["python3", "scripts/validate_repository.py"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertIn("Repository validation passed", completed.stdout)


if __name__ == "__main__":
    unittest.main()
