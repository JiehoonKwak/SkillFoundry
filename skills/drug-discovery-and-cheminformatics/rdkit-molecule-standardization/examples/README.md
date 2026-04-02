# Examples

- `sodium_acetate.smiles` contains the deterministic salt-form example used by the committed asset.
- Example command:
  `slurm/envs/chem-tools/bin/python skills/drug-discovery-and-cheminformatics/rdkit-molecule-standardization/scripts/standardize_rdkit_molecule.py --smiles-file skills/drug-discovery-and-cheminformatics/rdkit-molecule-standardization/examples/sodium_acetate.smiles --name sodium-acetate`
- Intended outcome after salt stripping and uncharging: acetic acid, `CC(=O)O`
