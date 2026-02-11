# Changelog

All notable changes to the Scanpy Nextflow Pipeline will be documented in this file.

## [1.0.0] - 2026-02-11

### Added
- Initial release of the Scanpy single-cell RNA-seq Nextflow pipeline
- Complete DSL2 Nextflow workflow for scRNA-seq analysis
- Modular architecture with 8 independent process modules:
  - LOAD_AND_QC: Data loading and quality control
  - NORMALIZE_AND_HVG: Normalization and highly variable gene selection
  - DIMENSIONALITY_REDUCTION: PCA and UMAP
  - CLUSTERING: Leiden clustering algorithm
  - MARKER_GENES: Differential expression and marker identification
  - ANNOTATE_CELL_TYPES: Cell type annotation
  - VISUALIZATION: Publication-quality plot generation
  - GENERATE_REPORT: HTML/Markdown summary reports
- Comprehensive configuration file with customizable parameters
- Support for multiple execution profiles (local, Docker, Singularity, SLURM)
- Automatic resource management and error handling
- Resume capability for interrupted runs
- Complete documentation:
  - README.md: Full pipeline documentation
  - QUICKSTART.md: Quick start guide
  - WORKFLOW.md: Detailed workflow description
  - PIPELINE_DIAGRAM.txt: Visual workflow diagram
  - CHANGELOG.md: Version history
- Example execution script (run_example.sh)
- Git ignore configuration

### Features
- Quality control with customizable thresholds
- Normalization to 10,000 counts per cell
- Highly variable gene selection (default 2000 genes)
- PCA with 50 components
- UMAP dimensionality reduction
- Leiden clustering with adjustable resolution
- Wilcoxon rank-sum test for marker genes
- Cell type annotation based on canonical markers
- Multiple visualization outputs:
  - UMAP plots (clusters, cell types, QC metrics)
  - Marker gene heatmaps
  - Dotplots for canonical markers
  - QC violin plots
  - PCA variance plots
- Automated HTML and Markdown report generation
- Execution tracking and logging

### Default Parameters
- min_genes: 200
- max_genes: 2500
- max_mt_percent: 5
- n_top_genes: 2000
- n_pcs: 50
- n_neighbors: 10
- n_pcs_umap: 40
- leiden_resolution: 1.0
- annotate_celltypes: true

### Resource Defaults
- max_cpus: 4
- max_memory: 16 GB
- max_time: 4 hours
- Retry strategy: 2 retries per process

### Supported Input Formats
- 10X Genomics format (matrix.mtx, genes.tsv, barcodes.tsv)

### Supported Output Formats
- HDF5 AnnData (.h5ad)
- CSV tables
- PNG visualizations
- HTML reports
- Markdown documents

### Dependencies
- Nextflow >= 23.04.0
- Python >= 3.9
- Scanpy >= 1.9.0
- igraph and leidenalg for clustering
- matplotlib, pandas, openpyxl for analysis and visualization

## Future Enhancements (Planned)

### [1.1.0] - Planned
- Support for additional input formats (Seurat, Loom, CSV)
- Batch correction module (Harmony, BBKNN, Scanorama)
- Cell cycle scoring
- Trajectory inference (PAGA, RNA velocity)
- Differential expression between conditions
- Gene set enrichment analysis
- Interactive HTML reports with embedded plots
- Multi-sample integration workflow

### [1.2.0] - Planned
- GPU acceleration support
- Automatic optimal parameter selection
- Cell type annotation using reference databases (CellTypist, scType)
- Support for spatial transcriptomics data
- Advanced QC metrics (doublet detection with Scrublet)
- Pseudobulk aggregation
- Advanced visualization (interactive UMAP with Plotly)

### [1.3.0] - Planned
- Integration with public databases (CellxGene, Human Cell Atlas)
- Automated benchmarking and quality assessment
- Cloud execution support (AWS Batch, Google Cloud, Azure)
- Nextflow Tower integration
- Container definitions (Dockerfile, Singularity recipe)
- Comprehensive test suite

## Version Numbering

This project follows [Semantic Versioning](https://semver.org/):
- MAJOR version: Incompatible API changes
- MINOR version: New functionality (backwards compatible)
- PATCH version: Bug fixes (backwards compatible)

## Contributing

To contribute to this pipeline:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Update this CHANGELOG
5. Submit a pull request

## License

MIT License - See LICENSE file for details
