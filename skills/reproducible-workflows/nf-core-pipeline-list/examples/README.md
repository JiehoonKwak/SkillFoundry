# Example Invocation

Use the wrapper script instead of calling `nf-core pipelines list --json` directly, because the CLI currently emits malformed multiline JSON in this environment.

```bash
python3 skills/reproducible-workflows/nf-core-pipeline-list/scripts/list_nfcore_pipelines.py \
  --sort pulled \
  --limit 5
```
