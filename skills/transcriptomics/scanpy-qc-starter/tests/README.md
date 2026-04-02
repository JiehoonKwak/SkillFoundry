# Tests

- Preflight smoke test: `python3 skills/transcriptomics/scanpy-qc-starter/scripts/preflight_counts.py --input skills/transcriptomics/scanpy-qc-starter/examples/toy_counts.tsv`
- Runtime QC smoke test: `slurm/envs/scanpy/bin/python skills/transcriptomics/scanpy-qc-starter/scripts/run_scanpy_qc.py --input skills/transcriptomics/scanpy-qc-starter/examples/toy_counts.tsv --summary-out scratch/scanpy/summary.json --h5ad-out scratch/scanpy/toy_counts.h5ad`
