from __future__ import annotations

import json
import unittest
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parents[1]


class SkillContractMetadataTests(unittest.TestCase):
    def test_metadata_and_raw_toy_input_stay_aligned(self) -> None:
        metadata = json.loads((SKILL_DIR / 'metadata.yaml').read_text(encoding='utf-8'))
        example = json.loads((SKILL_DIR / 'examples' / 'toy_input.json').read_text(encoding='utf-8'))
        skill_text = (SKILL_DIR / 'SKILL.md').read_text(encoding='utf-8')

        self.assertEqual(
            sorted(item['path'] for item in metadata['deliverables']),
            ['niche_labels.tsv', 'niche_markers.tsv', 'tissue_niche_report.md'],
        )
        self.assertEqual(
            metadata['starter_qc_files'],
            [
                'qc_spatial_graph.tsv',
                'qc_neighborhood_profiles.tsv',
                'qc_niche_scores.tsv',
                'run_summary.json',
            ],
        )
        self.assertNotIn('deliverables', example)
        self.assertEqual(example['expected_invariants']['cell_count'], len(example['cells']))
        self.assertEqual(example['expected_invariants']['gene_count'], len(example['genes']))
        self.assertEqual(example['expected_invariants']['cell_type_count'], len(example['cell_types']))
        self.assertEqual(example['expected_invariants']['niche_count'], len(example['candidate_niches']))

        for niche_name in example['candidate_niches']:
            self.assertEqual(
                sorted(example['niche_prototypes'][niche_name]['composition'].keys()),
                sorted(example['cell_types']),
            )
        for cell in example['cells']:
            self.assertEqual(sorted(cell['expression'].keys()), sorted(example['genes']))
            self.assertIn(cell['cell_type'], example['cell_types'])

        self.assertTrue((SKILL_DIR / 'scripts' / 'check_runtime.py').exists())
        self.assertTrue((SKILL_DIR / 'scripts' / 'plan_niche_run.py').exists())
        self.assertEqual(
            [item['path'] for item in metadata['helper_scripts']],
            [
                'scripts/check_runtime.py',
                'scripts/plan_niche_run.py',
                'scripts/run_tissue_niche_annotation.py',
            ],
        )
        self.assertIn('## Real-run directory contract', skill_text)
        self.assertIn('## Step 0. Check the runtime before touching the data', skill_text)
        self.assertIn('## Step 2. Scaffold the niche run plan before graph construction', skill_text)
        self.assertIn('## Missing-tool fallback rules', skill_text)
        self.assertIn('annotation_table.tsv', skill_text)
        self.assertIn('conda', skill_text)


if __name__ == '__main__':
    unittest.main()
