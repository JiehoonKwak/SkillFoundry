# Example Data

`toy_measurements.tsv` is a tiny deterministic dataset used to render the starter figure.

```bash
slurm/envs/statistics/bin/python skills/visualization-and-reporting/matplotlib-publication-plot-starter/scripts/render_publication_plot.py \
  --input skills/visualization-and-reporting/matplotlib-publication-plot-starter/examples/toy_measurements.tsv \
  --png-out scratch/matplotlib/publication_plot.png \
  --summary-out scratch/matplotlib/publication_plot_summary.json
```
