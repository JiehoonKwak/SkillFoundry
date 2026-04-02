from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SKILL_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = SKILL_ROOT / "scripts" / "fetch_uniprot_sequence_feature_summary.py"


def load_module():
    spec = importlib.util.spec_from_file_location("uniprot_sequence_feature_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


MODULE = load_module()


class UniProtSequenceFeatureSummaryTests(unittest.TestCase):
    def test_build_summary_collapses_features_deterministically(self) -> None:
        payload = {
            "primaryAccession": "PTEST",
            "uniProtkbId": "PTEST_HUMAN",
            "annotationScore": 5.0,
            "proteinDescription": {"recommendedName": {"fullName": {"value": "Example protein"}}},
            "genes": [{"geneName": {"value": "GENE1"}, "synonyms": [{"value": "ALT1"}]}],
            "organism": {"scientificName": "Homo sapiens", "commonName": "Human", "taxonId": 9606},
            "sequence": {"length": 42},
            "features": [
                {
                    "type": "Natural variant",
                    "location": {"start": {"value": 17, "modifier": "EXACT"}, "end": {"value": 17, "modifier": "EXACT"}},
                    "description": "Example variant",
                    "featureId": "VAR_000001",
                    "alternativeSequence": {"originalSequence": "A", "alternativeSequences": ["V"]},
                },
                {
                    "type": "Chain",
                    "location": {"start": {"value": 1, "modifier": "EXACT"}, "end": {"value": 42, "modifier": "EXACT"}},
                    "description": "Example chain",
                    "featureId": "PRO_000001",
                },
                {
                    "type": "Motif",
                    "location": {"start": {"value": 5, "modifier": "EXACT"}, "end": {"value": 10, "modifier": "EXACT"}},
                    "description": "Activation motif",
                },
                {
                    "type": "Region",
                    "location": {"start": {"value": 2, "modifier": "EXACT"}, "end": {"value": 20, "modifier": "EXACT"}},
                    "description": "Binding region",
                },
                {
                    "type": "Modified residue",
                    "location": {"start": {"value": 8, "modifier": "EXACT"}, "end": {"value": 8, "modifier": "EXACT"}},
                    "description": "Phosphoserine",
                },
            ],
        }

        summary = MODULE.build_summary(payload, feature_limit=4)

        self.assertEqual(summary["accession"], "PTEST")
        self.assertEqual(summary["recommended_name"], "Example protein")
        self.assertEqual(summary["gene_names"], ["GENE1", "ALT1"])
        self.assertEqual(summary["sequence_length"], 42)
        self.assertEqual(summary["feature_count"], 5)
        self.assertEqual(summary["feature_type_counts"]["Natural variant"], 1)
        self.assertEqual(
            [feature["type"] for feature in summary["representative_features"]],
            ["Chain", "Region", "Motif", "Modified residue"],
        )
        self.assertEqual(summary["representative_features"][0]["span"], "1-42")

    def test_normalize_feature_includes_variation(self) -> None:
        feature = {
            "type": "Natural variant",
            "location": {"start": {"value": 9, "modifier": "EXACT"}, "end": {"value": 9, "modifier": "EXACT"}},
            "description": "Example variant",
            "featureId": "VAR_123",
            "alternativeSequence": {"originalSequence": "Q", "alternativeSequences": ["H"]},
        }

        normalized = MODULE.normalize_feature(feature)

        self.assertEqual(normalized["type"], "Natural variant")
        self.assertEqual(normalized["span"], "9")
        self.assertEqual(normalized["variation"], "Q->H")

    def test_live_fetch_writes_json_for_p04637(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            out_path = Path(tmp_dir) / "p04637.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--accession",
                    "P04637",
                    "--out",
                    str(out_path),
                ],
                check=True,
                capture_output=True,
                text=True,
                timeout=120,
            )

            payload = json.loads(out_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["accession"], "P04637")
            self.assertEqual(payload["entry_id"], "P53_HUMAN")
            self.assertEqual(payload["organism_scientific_name"], "Homo sapiens")
            self.assertIn("TP53", payload["gene_names"])
            self.assertEqual(payload["sequence_length"], 393)
            self.assertGreaterEqual(payload["feature_count"], 1000)
            self.assertGreaterEqual(payload["feature_type_counts"]["Natural variant"], 1000)
            self.assertEqual(payload["representative_features"][0]["type"], "Chain")
            self.assertEqual(payload["representative_features"][0]["span"], "1-393")

    def test_rejects_invalid_accession_format(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--accession",
                "P0_NOT_REAL",
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Invalid UniProt accession format", completed.stderr + completed.stdout)


if __name__ == "__main__":
    unittest.main()
