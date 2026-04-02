#!/bin/bash
#SBATCH --job-name=smoke
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --error=slurm/logs/%x-%j.err
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -euo pipefail
echo "hostname: $(hostname)"
echo "started: $(date --iso-8601=seconds)"
echo hello
