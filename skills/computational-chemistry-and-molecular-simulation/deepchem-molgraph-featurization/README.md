# DeepChem MolGraph Featurization

This starter skill turns a small TSV of `molecule_id` and `smiles` rows into a structured DeepChem molecular graph summary.

## What It Does

- Featurizes each SMILES string with `DeepChem`'s `MolGraphConvFeaturizer`
- Records graph size information for each molecule
- Preserves basic provenance such as the input file path and DeepChem version

## Typical Use

Run the repo-managed `chemtools` prefix:

```bash
slurm/envs/chemtools/bin/python skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/scripts/featurize_molecules.py \
  --input skills/computational-chemistry-and-molecular-simulation/deepchem-molgraph-featurization/examples/molecules.tsv \
  --out scratch/deepchem/featurization.json
```

## Benchmark Notes

The layer 2 benchmark compares this maintained wrapper against an ad hoc inline DeepChem script that computes the same graphs but omits wrapper-level provenance and feature-dimension metadata.
