# References

## Canonical model and documentation anchors

- `esm-github`: official ESM repository with ESM-1/ESM-2 lineage and usage notes
- `esm-transformers-docs`: Hugging Face Transformers documentation for loading ESM checkpoints through a maintained Python API
- `esm-foundation-model-paper`: Rives et al. 2021, foundational protein language model paper showing structure and function signals in unsupervised sequence embeddings
- `prottrans-github`: official ProtTrans repository for ProtBert, ProtT5, and related transfer-learning workflows
- `prottrans-paper`: Elnaggar et al. 2022, large-scale transfer-learning paper for protein sequences
- `prott5-model-card`: maintained model-card surface for `Rostlab/prot_t5_xl_half_uniref50-enc`
- `deepfri-github`: canonical downstream function-prediction repository for sequence-to-function analysis

## Practical implementation notes

- ESM-2 and ProtT5 are the two primary real-model backends exposed by this skill.
- ESM-family tokenization uses raw amino-acid strings.
- ProtT5-style tokenization usually expects amino acids separated by spaces, with uncommon residues normalized to `X`.
- The default repo smoke path uses the bundled deterministic `mock` backend so tests remain sandbox-safe and do not depend on `torch`, `transformers`, or model downloads.

## Suggested escalation path

1. Run the bundled toy example with the `mock` backend to validate the contract.
2. Swap in a real FASTA plus labels from a curated benchmark or in-house dataset.
3. Move to `--backend transformers` with ESM-2 or ProtT5 in a dedicated environment.
4. Use the emitted embedding table as input to downstream DeepFRI, nearest-neighbor transfer, or supervised evaluation.

## Slurm note

For large checkpoints, run the same CLI under a batch job such as [`slurm/jobs/protein-language-model-function-analysis-smoke.sbatch`](../../../slurm/jobs/protein-language-model-function-analysis-smoke.sbatch) and replace the `mock` config with a real `transformers` backend plus the appropriate environment activation.
