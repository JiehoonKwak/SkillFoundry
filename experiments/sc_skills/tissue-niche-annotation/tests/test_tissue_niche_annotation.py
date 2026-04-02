from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / 'scripts' / 'run_exercise.py'
CHECK_RUNTIME = SKILL_DIR / 'scripts' / 'check_runtime.py'
PLAN_RUN = SKILL_DIR / 'scripts' / 'plan_niche_run.py'


class SkillLocalToyRunTests(unittest.TestCase):
    def test_runtime_check_reports_missing_modules_and_fallbacks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / 'runtime_check.json'
            completed = subprocess.run(
                ['python3', str(CHECK_RUNTIME), '--out', str(out_path)],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(out_path.read_text(encoding='utf-8'))
            self.assertIn('missing_modules', payload)
            self.assertIn('bootstrap_commands', payload)
            self.assertIn('fallback_guidance', payload)
            self.assertGreaterEqual(len(payload['bootstrap_commands']), 2)
            self.assertIn('conda_available', payload)

    def test_run_planner_emits_expected_outputs_and_stop_conditions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / 'niche_run_plan.json'
            completed = subprocess.run(
                [
                    'python3',
                    str(PLAN_RUN),
                    '--adata',
                    'query.h5ad',
                    '--cell-type-column',
                    'cell_type',
                    '--x-key',
                    'x_coord',
                    '--y-key',
                    'y_coord',
                    '--sample-key',
                    'sample_id',
                    '--graph-mode',
                    'knn',
                    '--k',
                    '12',
                    '--out',
                    str(out_path),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            payload = json.loads(out_path.read_text(encoding='utf-8'))
            self.assertEqual(payload['graph']['mode'], 'knn')
            self.assertEqual(payload['graph']['k'], 12)
            self.assertIn('03_spatial_graph.tsv', payload['expected_outputs'])
            self.assertIn('missing spatial coordinates', payload['stop_conditions'])

    def test_toy_run_computes_graph_derived_niches(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            completed = subprocess.run(
                ['python3', str(SCRIPT), '--outdir', tmpdir],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            outdir = Path(tmpdir)
            labels = pd.read_csv(outdir / 'niche_labels.tsv', sep='\t').set_index('cell_id')
            self.assertEqual(labels.shape[0], 9)
            self.assertEqual(
                labels['niche_label'].value_counts().to_dict(),
                {
                    'epithelial_core': 3,
                    'immune_aggregate': 4,
                    'immune_epithelial_interface': 2,
                },
            )
            self.assertEqual(labels.loc['cell_004', 'base_best_niche'], 'epithelial_core')
            self.assertEqual(labels.loc['cell_004', 'niche_label'], 'immune_epithelial_interface')
            self.assertTrue(bool(labels.loc['cell_004', 'graph_shifted']))
            self.assertGreater(float(labels.loc['cell_004', 'boundary_score']), 0.75)
            self.assertLess(float(labels.loc['cell_001', 'boundary_score']), 0.35)
            self.assertEqual(labels.loc['cell_006', 'niche_label'], 'immune_aggregate')
            self.assertGreater(float(labels.loc['cell_006', 'niche_confidence']), 0.7)
            self.assertEqual(
                labels.loc[labels['cell_type'] == 'Epithelial', 'niche_label'].nunique(),
                2,
            )

            markers = pd.read_csv(outdir / 'niche_markers.tsv', sep='\t').set_index('niche_label')
            self.assertIn('EPCAM', markers.loc['epithelial_core', 'marker_context'])
            self.assertIn('CXCL9', markers.loc['immune_aggregate', 'marker_context'])
            self.assertIn('CXCL10', markers.loc['immune_epithelial_interface', 'marker_context'])
            self.assertGreater(float(markers.loc['epithelial_core', 'support_score']), 0.55)
            self.assertGreater(float(markers.loc['immune_aggregate', 'support_score']), 0.55)

            graph = pd.read_csv(outdir / 'qc_spatial_graph.tsv', sep='\t')
            row_sums = graph.groupby('source_cell')['weight'].sum().to_numpy(dtype=float)
            self.assertTrue(np.allclose(row_sums, 1.0, atol=1e-6))
            self.assertEqual(set(graph['relation']), {'neighbor', 'self'})

            profiles = pd.read_csv(outdir / 'qc_neighborhood_profiles.tsv', sep='\t').set_index('cell_id')
            local_columns = [column for column in profiles.columns if column.startswith('local_') and column.endswith('_fraction')]
            smoothed_columns = [column for column in profiles.columns if column.startswith('smoothed_') and column.endswith('_fraction')]
            self.assertTrue(np.allclose(profiles[local_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6))
            self.assertTrue(np.allclose(profiles[smoothed_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6))
            self.assertGreater(float(profiles.loc['cell_005', 'boundary_score']), 0.75)
            self.assertEqual(profiles.loc['cell_006', 'dominant_neighbor_type'], 'T cell')

            scores = pd.read_csv(outdir / 'qc_niche_scores.tsv', sep='\t').set_index('cell_id')
            probability_columns = [column for column in scores.columns if column.endswith('_probability')]
            self.assertTrue(np.allclose(scores[probability_columns].sum(axis=1).to_numpy(dtype=float), 1.0, atol=1e-6))
            self.assertEqual(scores.loc['cell_004', 'base_best_niche'], 'epithelial_core')
            self.assertEqual(scores.loc['cell_004', 'assigned_niche'], 'immune_epithelial_interface')
            self.assertTrue(bool(scores.loc['cell_004', 'graph_shifted']))
            self.assertGreater(float(scores.loc['cell_004', 'confidence_margin']), 0.1)

            summary = json.loads((outdir / 'run_summary.json').read_text(encoding='utf-8'))
            self.assertEqual(
                summary['method_steps'],
                [
                    'spatial_graph_construction',
                    'neighborhood_composition',
                    'graph_smoothed_niche_scoring',
                    'marker_summarization',
                ],
            )
            self.assertIn('cell_004', summary['shifted_cells'])
            self.assertIn('cell_005', summary['boundary_cells'])
            self.assertEqual(summary['niche_counts']['immune_epithelial_interface'], 2)
            self.assertIn('run_summary.json', summary['written_files'])


if __name__ == '__main__':
    unittest.main()
