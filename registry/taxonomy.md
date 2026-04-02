# SciSkillUniverse Taxonomy

## Scientific knowledge access and method extraction
- Literature search
- Paper triage and ranking
- Review-paper mining
- Benchmark-paper mining
- Method-section extraction
- Figure / table extraction
- Citation chaining
- Reproducibility cue extraction
- Dataset / code / package link extraction
- Protocol and workflow extraction
- Scientific summarization
- Claim / assumption / limitation extraction

## Data acquisition and dataset handling
- Public dataset discovery
- Metadata harmonization
- Format conversion
- Compression / decompression
- Chunking / sharding
- Data validation
- Data provenance tracking
- Synthetic toy dataset generation for tests

## Genomics
- FASTQ / BAM / CRAM basics
- Read QC and trimming
- Alignment and mapping
- Quantification
- Variant calling
- Variant filtering
- Annotation and effect prediction
- Somatic pipelines
- Germline pipelines
- GWAS
- Fine mapping
- Polygenic risk scoring
- Comparative genomics
- Phylogenomics
- Long-read genomics
- Metagenomics

## Transcriptomics
- Bulk RNA-seq QC and normalization
- Differential expression
- Isoform / transcript-level analysis
- Single-cell RNA-seq preprocessing
- Single-cell integration / batch correction
- Cell type annotation
- Trajectory inference
- RNA velocity
- Pseudobulk analysis
- Perturb-seq
- Spatial transcriptomics
- Multi-sample / atlas workflows

## Epigenomics and chromatin
- ATAC-seq
- ChIP-seq
- CUT&RUN / CUT&Tag
- Methylation analysis
- Peak calling
- Motif analysis
- Footprinting
- Chromatin interaction / Hi-C
- Multiome integration

## Proteomics and protein biology
- Proteomics dataset discovery
- Protein accession and metadata lookup
- MS proteomics preprocessing
- Protein identification / quantification
- Differential proteomics
- PTM analysis
- Protein families and domains
- Sequence feature annotation
- Protein embeddings
- Protein structure cross-links
- Interface analysis
- Sequence-to-function modeling

## Metabolomics and other omics
- Metabolomics preprocessing
- Feature annotation
- Differential metabolomics
- Lipidomics
- Microbiome analysis
- Multi-omics integration

## Structural biology and molecular modeling
- Structure parsing and cleanup
- PDB / mmCIF utilities
- Visualization
- Fold comparison
- Binding site analysis
- Sequence/structure alignment
- Structure-to-sequence cross-linking
- Protein complex metadata

## Systems biology and network science
- Gene set enrichment
- Pathway analysis
- Reactome event lookup
- Reactome identifier enrichment
- Pathway traversal and hierarchy walks
- Gene-set tooling from Bioconductor
- Regulatory network inference
- Protein-protein interaction analysis
- Graph construction
- Network propagation
- Causal network analysis

## Imaging and phenotype analysis
- Microscopy pipelines
- Segmentation
- Cell tracking
- Pathology / histology workflows
- Feature extraction
- Representation learning for scientific images
- Multi-modal image-omics integration

## Drug discovery and cheminformatics
- Molecule standardization
- Compound identifier lookup
- Descriptor computation
- Fingerprints and similarity search
- Scaffold analysis
- QSAR / property prediction
- Virtual screening
- ADMET modeling
- Molecular generation
- Benchmarks and evaluation

## Computational chemistry and molecular simulation
- Small-molecule conformer generation
- Force-field assignment
- Molecular dynamics setup
- Molecular dynamics execution
- Trajectory analysis
- Free-energy workflows
- Quantum chemistry setup
- Single-point energy calculations
- Geometry optimization
- QM/MM workflows

## Materials science and engineering
- Crystal structure parsing
- Composition featurization
- Materials-property prediction
- Surrogate models for simulation
- Materials benchmark datasets
- Generative design for materials
- Phase / stability analysis
- Laboratory-to-simulation data linkage

## Earth, climate, and geospatial science
- Raster / vector ingestion
- Remote sensing preprocessing
- Time-series environmental data handling
- Climate reanalysis access
- Geospatial feature engineering
- Physics-informed forecasting
- Spatial interpolation and uncertainty
- Scientific map and dashboard generation

## Clinical / biomedical data science
- EHR preprocessing
- Survival analysis
- Time-series clinical modeling
- Cohort extraction
- Phenotyping
- Missing-data handling
- Fairness / bias analysis
- Privacy-preserving analysis

## Statistical and machine learning foundations for science
- Experimental design
- Statistical testing
- Bayesian workflows
- Causal inference
- Dimensionality reduction
- Clustering
- Representation learning
- Foundation models for biology / chemistry / science
- Multimodal fusion
- Uncertainty estimation
- Active learning
- Bayesian optimization
- Simulation-based inference

## Scientific agents and automation
- Literature-to-hypothesis agents
- Tool-using analysis agents
- Planning and execution agents
- Workflow generation agents
- Code-generation agents for scientific tasks
- Evaluation harnesses for scientific agents
- Long-horizon agent memory and provenance
- Agent orchestration over skills and registries

## Reproducible workflows and workflow engines
- Snakemake
- Nextflow / nf-core
- CWL command-line tools
- CWL workflows
- WDL tasks
- WDL workflows
- Workflow interoperability and translation
- Docker / Apptainer / Singularity
- Conda / Mamba / Pixi / uv
- Environment locking
- Reproducible notebooks
- CI for scientific pipelines

## HPC, Slurm, and scaling
- Interactive debug workflows
- Batch jobs
- Job arrays
- GPU jobs
- Multi-node jobs
- CPU / memory tuning
- I/O-aware workflows
- Scratch-space usage
- Cluster-safe installs
- Monitoring and accounting

## Visualization and reporting
- Publication plots
- Dashboards
- Interactive reports
- Notebook-to-report conversion
- Figure recreation from papers
- Summary pages and catalogs
- Skill browser / mindmap generation

## Meta-maintenance
- Taxonomy refinement
- Resource deduplication
- Skill deduplication
- Freshness audits
- Broken-link audits
- Regression testing
- Deprecation / migration
- Documentation quality improvement

## Neuroscience and neuroimaging
- Neuroimaging I/O and formats
- fMRI preprocessing and denoising
- EEG / MEG preprocessing
- Connectomics and graph analysis
- Spike sorting and electrophysiology
- Neural decoding and encoding models
- Brain atlas integration
- Multimodal neuroimaging fusion
- Computational psychiatry and clinical neuroinformatics

## Physics and astronomy
- PDE / CFD simulation workflows
- Operator learning and surrogate physics models
- Inverse problems and scientific reconstruction
- Differentiable simulation
- High-energy / detector data analysis
- Cosmology and sky-survey analysis
- Telescope image preprocessing
- Spectroscopy and signal extraction
- Symbolic discovery of physical laws
- Scientific time-series anomaly detection

## Ecology, evolution, and biodiversity
- Species distribution modeling
- Biodiversity dataset discovery
- Population dynamics and ecological forecasting
- Phylogenetic comparative workflows
- Conservation genomics and eDNA
- Wildlife sensing and bioacoustics
- Ecosystem remote sensing
- Climate-ecology integration

## Agriculture, food, and plant science
- Plant phenotyping
- Crop genomics and breeding
- Precision agriculture sensing
- Soil / rhizosphere / microbiome analysis
- Yield forecasting
- Disease and stress detection
- Agronomic experiment design
- Food chemistry and quality modeling

## Robotics, lab automation, and scientific instrumentation
- Liquid-handling protocol generation
- Robotic experiment planning
- Active-learning experiment loops
- Instrument control and scheduling
- Automated microscopy acquisition
- ELN / LIMS / sample-tracking integration
- Laboratory robotics safety and monitoring
- Closed-loop optimization for experiments

## Scientific computing and numerical methods
- ODE / SDE simulation workflows
- PDE discretization and solvers
- Sparse / iterative linear algebra
- Automatic differentiation for scientific models
- Differentiable programming for simulators
- Inverse problems and data assimilation
- Uncertainty-aware simulation
- Numerical benchmarking and verification
