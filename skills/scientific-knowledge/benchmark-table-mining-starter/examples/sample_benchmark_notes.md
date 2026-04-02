# Example Benchmark Section

The paper first summarizes dataset scale, then reports a leaderboard snapshot.

| Dataset | Cells | Species |
| --- | --- | --- |
| PBMC 68k | 68000 | Human |
| Tabula Sapiens | 480000 | Human |
| Norman 2019 | 105000 | Human |

The leaderboard below is the table that should be mined for benchmark rows.

| Task | Dataset | Metric | Model | Score |
| --- | --- | --- | --- | --- |
| Cell type annotation | PBMC 68k | Macro F1 | scBERT | 0.921 |
| Cell type annotation | Tabula Sapiens | Macro F1 | Geneformer | 0.908 |
| Perturbation response | Norman 2019 | AUROC | GEARS | 0.873 |
