# Scanpy Single-Cell RNA-seq Nextflow Pipeline

A comprehensive Nextflow pipeline for single-cell RNA-seq analysis using Scanpy.

## Features

- Fully automated workflow from raw 10X data to annotated cell types
- Quality control with customizable thresholds
- Normalization and highly variable gene selection
- PCA and UMAP dimensionality reduction
- Leiden clustering
- Marker gene identification
- Cell type annotation
- Comprehensive visualizations
- HTML and Markdown reports

## Environment Management

The pipeline supports **4 execution modes** to match your needs:

1. **Standard** (default) - Uses your current Python environment
   ```bash
   nextflow run main.nf --input data --outdir results
   ```

2. **Conda** (recommended) - Nextflow manages conda environment automatically
   ```bash
   nextflow run main.nf -profile conda --input data --outdir results
   ```

3. **Docker** - Maximum reproducibility with containers
   ```bash
   nextflow run main.nf -profile docker --input data --outdir results
   ```

4. **Singularity** - For HPC clusters
   ```bash
   nextflow run main.nf -profile singularity --input data --outdir results
   ```

📖 **See [ENVIRONMENT.md](ENVIRONMENT.md) for complete setup instructions** including:
- Detailed installation steps for each mode
- Conda environment.yml
- Dockerfile and build scripts
- HPC cluster configurations
- Troubleshooting guide

## Requirements

### Software
- Nextflow >= 23.04.0
- Python >= 3.9
- Scanpy >= 1.9.0

### Python packages
```bash
pip install scanpy matplotlib pandas openpyxl
pip install igraph leidenalg  # For clustering
```

## Quick Start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
./nextflow -version
```

### 2. Run the pipeline

```bash
# Basic usage with default parameters
nextflow run main.nf --input /path/to/10x/data --outdir results

# With custom parameters
nextflow run main.nf \\
    --input /path/to/10x/data \\
    --outdir results \\
    --min_genes 200 \\
    --max_genes 2500 \\
    --max_mt_percent 5 \\
    --leiden_resolution 1.0 \\
    --annotate_celltypes true
```

### 3. Resume a failed run

```bash
nextflow run main.nf --input /path/to/10x/data --outdir results -resume
```

## Input Data

The pipeline expects 10X Genomics format data:
```
input_directory/
├── matrix.mtx      # Expression matrix
├── genes.tsv       # Gene information
└── barcodes.tsv    # Cell barcodes
```

## Parameters

### Input/Output
- `--input`: Path to 10X directory (required)
- `--outdir`: Output directory (default: `./results`)

### QC Parameters
- `--min_genes`: Minimum genes per cell (default: 200)
- `--min_cells`: Minimum cells per gene (default: 3)
- `--max_genes`: Maximum genes per cell (default: 2500)
- `--max_mt_percent`: Maximum mitochondrial percentage (default: 5)

### Analysis Parameters
- `--n_top_genes`: Number of highly variable genes (default: 2000)
- `--n_pcs`: Number of principal components (default: 50)
- `--n_neighbors`: Number of neighbors for graph (default: 10)
- `--n_pcs_umap`: PCs used for UMAP (default: 40)
- `--leiden_resolution`: Leiden clustering resolution (default: 1.0)
- `--annotate_celltypes`: Enable cell type annotation (default: true)

### Resource Parameters
- `--max_cpus`: Maximum CPUs (default: 4)
- `--max_memory`: Maximum memory (default: '16.GB')
- `--max_time`: Maximum time (default: '4.h')

## Output Structure

```
results/
├── 01_qc/
│   ├── adata_qc.h5ad
│   ├── qc_metrics.csv
│   └── qc_violin.png
├── 02_normalize/
│   ├── adata_normalized.h5ad
│   └── hvg_plot.png
├── 03_dimred/
│   ├── adata_dimred.h5ad
│   └── pca_variance.png
├── 04_clustering/
│   ├── adata_clustered.h5ad
│   └── cluster_sizes.csv
├── 05_markers/
│   ├── adata_with_markers.h5ad
│   ├── marker_genes.csv
│   └── marker_genes_heatmap.png
├── 06_annotation/
│   ├── adata_annotated.h5ad
│   └── cell_type_counts.csv
├── 07_visualization/
│   ├── umap_clusters.png
│   ├── umap_cell_types.png
│   ├── umap_qc_metrics.png
│   └── marker_dotplot.png
├── analysis_summary.html
├── analysis_summary.md
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Execution Profiles

### Local (default)
```bash
nextflow run main.nf --input data --outdir results
```

### Docker
```bash
nextflow run main.nf --input data --outdir results -profile docker
```

### Singularity
```bash
nextflow run main.nf --input data --outdir results -profile singularity
```

### SLURM
```bash
nextflow run main.nf --input data --outdir results -profile slurm
```

## Customizing Cell Type Annotation

To customize cell type annotations for your specific dataset, edit the `modules/annotate.nf` file:

```python
cell_type_map = {
    '0': 'Your Cell Type 1',
    '1': 'Your Cell Type 2',
    # ... add your annotations
}
```

Base your annotations on the marker genes identified in `05_markers/marker_genes.csv`.

## Advanced Usage

### Running specific steps only

You can modify the `main.nf` to skip certain steps by commenting them out.

### Adjusting cluster resolution

Higher resolution = more clusters:
```bash
nextflow run main.nf --input data --leiden_resolution 1.5
```

Lower resolution = fewer clusters:
```bash
nextflow run main.nf --input data --leiden_resolution 0.5
```

### Parallel execution on HPC

Create a custom profile in `nextflow.config`:

```groovy
profiles {
    hpc {
        process.executor = 'slurm'
        process.queue = 'normal'
        process.cpus = 8
        process.memory = '32.GB'
        process.time = '8.h'
    }
}
```

Run with: `nextflow run main.nf -profile hpc`

## Troubleshooting

### Out of memory errors
Increase memory: `--max_memory '32.GB'`

### Clustering fails
Install igraph: `pip install igraph leidenalg`

### Pipeline hangs
Check logs in `work/` directory

### Resume from checkpoint
Always use `-resume` flag to continue from where it stopped

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology.
- **Scanpy**: Wolf, F.A., et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biology.

## Support

For issues and questions:
- Check the [Scanpy documentation](https://scanpy.readthedocs.io/)
- Check the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

## License

MIT License
