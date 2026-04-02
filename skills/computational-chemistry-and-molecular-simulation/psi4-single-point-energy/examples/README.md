# Example

Run the bundled water single-point calculation:

```bash
slurm/envs/psi4/bin/python \
  skills/computational-chemistry-and-molecular-simulation/psi4-single-point-energy/scripts/run_psi4_single_point.py \
  --out skills/computational-chemistry-and-molecular-simulation/psi4-single-point-energy/assets/water_hf_sto3g_summary.json
```

Override the method and basis for exploration:

```bash
slurm/envs/psi4/bin/python \
  skills/computational-chemistry-and-molecular-simulation/psi4-single-point-energy/scripts/run_psi4_single_point.py \
  --method hf \
  --basis sto-3g
```
