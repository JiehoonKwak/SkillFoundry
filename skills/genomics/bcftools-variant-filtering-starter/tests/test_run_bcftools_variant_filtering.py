from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / "skills" / "genomics" / "bcftools-variant-filtering-starter" / "scripts" / "run_bcftools_variant_filtering.py"


class BcftoolsVariantFilteringTests(unittest.TestCase):
    def test_runtime_summary_matches_expected_variants(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "variant_filtering_summary.json"
            filtered_vcf = Path(tmp_dir) / "toy.filtered.vcf.gz"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    "--input",
                    "skills/genomics/bcftools-variant-filtering-starter/examples/toy_variants.vcf",
                    "--filtered-vcf-out",
                    str(filtered_vcf),
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
            self.assertEqual(payload["input_record_count"], 3)
            self.assertEqual(payload["passing_record_count"], 2)
            self.assertEqual(payload["kept_ids"], ["varA", "varC"])
            self.assertTrue(filtered_vcf.exists())
            self.assertTrue(Path(f"{filtered_vcf}.tbi").exists())


if __name__ == "__main__":
    unittest.main()
