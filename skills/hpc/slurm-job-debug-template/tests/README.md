# Tests

- Renderer smoke test: `python3 skills/hpc/slurm-job-debug-template/scripts/render_sbatch.py --command "echo hello" --job-name smoke`
- Cluster smoke test: `python3 skills/hpc/slurm-job-debug-template/scripts/submit_smoke_job.py --partition cpu --job-name smoke --sleep 1`
