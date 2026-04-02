# Task Inventory

This inventory maps every Spatial Agent reference under [`ref/skill/`](/home/shuaikes/projects/agent/SciSkillUniverse/ref/skill) to an experiment-only portable skill target under [`experiments/sc_skills/`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills).

| Task slug | Title | Family | Deliverables |
| --- | --- | --- | --- |
| annotation | Spatial Cell Annotation and Niche Labeling | spatial_foundation | annotation_table.tsv, marker_evidence.md, niche_summary.md |
| tissue-niche-annotation | Tissue Niche Annotation | spatial_structure | niche_labels.tsv, niche_markers.tsv, tissue_niche_report.md |
| cell-cell-communication | Cell-Cell Communication Analysis | communication | communication_results.tsv, priority_pairs.tsv, interpretation_report.md |
| cell-deconvolution | Tangram-Based Cell Deconvolution | tangram_extensions | mapping_matrix.h5ad, cell_assignment.tsv, deconvolution_report.md |
| cellphonedb-analysis | CellPhoneDB Ligand-Receptor Analysis | communication | cellphonedb_significant_means.tsv, interaction_ranked.tsv, cellphonedb_report.md |
| database-query | Bioinformatics Database Query | bioinformatics_utilities | query_results.tsv, source_provenance.md, resolution_notes.md |
| gene-imputation | Spatial Gene Imputation | tangram_extensions | imputed_expression.h5ad, heldout_metrics.tsv, imputation_report.md |
| liana-analysis | LIANA Communication Analysis | communication | liana_rankings.tsv, resource_overlap.tsv, liana_report.md |
| ligand-receptor-discovery | Ligand-Receptor Discovery | communication | candidate_pairs.tsv, evidence_table.tsv, discovery_summary.md |
| mapping-validation | Spatial Mapping Validation | tangram_extensions | validation_metrics.tsv, holdout_summary.md, mapping_validation_report.md |
| multimodal-integration | Multimodal Integration | multimodal_modeling | integrated_latent.h5mu, modality_qc.tsv, integration_report.md |
| panel-design | Spatial Panel Design | bioinformatics_utilities | panel_candidates.tsv, panel_rationale.md, platform_notes.md |
| sequence-analysis | Sequence Analysis and Primer Design | bioinformatics_utilities | sequence_summary.tsv, feature_annotations.tsv, primer_candidates.tsv |
| spatial-deconvolution | Spatial Deconvolution | spatial_deconvolution | cell_type_abundance.h5ad, model_comparison.tsv, deconvolution_summary.md |
| spatial-domain-detection | Spatial Domain Detection | spatial_structure | domain_labels.tsv, domain_markers.tsv, domain_detection_report.md |
| spatial-mapping | Spatial Mapping | spatial_foundation | mapped_labels.tsv, mapping_scores.tsv, spatial_mapping_report.md |
| squidpy-analysis | Squidpy Spatial Statistics | spatial_structure | spatial_stats.tsv, neighbor_enrichment.tsv, squidpy_report.md |
| trajectory-inference | Trajectory and Velocity Inference | trajectory_modeling | trajectory_coordinates.tsv, fate_probabilities.tsv, trajectory_report.md |

- The machine-readable source of truth for this inventory is [`batch_design_manifest.json`](/home/shuaikes/projects/agent/SciSkillUniverse/experiments/sc_skills/batch_design_manifest.json).
- Every experiment package is expected to carry docs plus `metadata.yaml`, `scripts/`, `examples/`, `tests/`, and `assets/` after hardening.

