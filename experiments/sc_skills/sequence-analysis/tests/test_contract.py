from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_toy_input_match_raw_contract(self) -> None:
        metadata = json.loads((SKILL_DIR / 'metadata.yaml').read_text(encoding='utf-8'))
        example = json.loads((SKILL_DIR / 'examples' / 'toy_input.json').read_text(encoding='utf-8'))

        self.assertNotIn('deliverables', example)
        self.assertEqual(
            [item['path'] for item in metadata['deliverables']],
            ['sequence_summary.tsv', 'feature_annotations.tsv', 'primer_candidates.tsv'],
        )
        self.assertEqual(example['expected_invariants']['query_count'], len(example['queries']))
        self.assertEqual(example['expected_invariants']['lookup_candidates_per_query'], len(example['references']))

        query_ids = [item['query_id'] for item in example['queries']]
        self.assertEqual(len(query_ids), len(set(query_ids)))

        reference_ids = [item['reference_id'] for item in example['references']]
        self.assertEqual(len(reference_ids), len(set(reference_ids)))

        for reference in example['references']:
            sequence_length = len(reference['sequence'])
            for feature in reference['features']:
                self.assertGreaterEqual(feature['start'], 1)
                self.assertLessEqual(feature['end'], sequence_length)
                self.assertLessEqual(feature['start'], feature['end'])
            for variant in reference['variants']:
                self.assertGreaterEqual(variant['position'], 1)
                self.assertLessEqual(variant['position'], sequence_length)
                self.assertNotEqual(variant['ref_base'], variant['alt_base'])


if __name__ == '__main__':
    unittest.main()
