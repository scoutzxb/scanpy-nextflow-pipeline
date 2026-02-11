# Scanpy Pipeline Workflow

## Pipeline Architecture

```
Input (10X Data)
    ↓
┌───────────────────────────────────────────────┐
│  LOAD_AND_QC                                  │
│  - Load 10X data (matrix.mtx, genes, barcodes)│
│  - Calculate QC metrics                       │
│  - Filter cells and genes                     │
│  Output: adata_qc.h5ad                        │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  NORMALIZE_AND_HVG                            │
│  - Normalize to 10,000 counts                 │
│  - Log transform                              │
│  - Identify highly variable genes             │
│  - Regress out confounders                    │
│  - Scale data                                 │
│  Output: adata_normalized.h5ad                │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  DIMENSIONALITY_REDUCTION                     │
│  - PCA (50 components)                        │
│  - Compute neighbor graph                     │
│  - UMAP embedding                             │
│  Output: adata_dimred.h5ad                    │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  CLUSTERING                                   │
│  - Leiden clustering                          │
│  - Generate cluster assignments               │
│  Output: adata_clustered.h5ad                 │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  MARKER_GENES                                 │
│  - Wilcoxon rank-sum test                     │
│  - Identify top markers per cluster           │
│  Output: adata_with_markers.h5ad              │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  ANNOTATE_CELL_TYPES (optional)               │
│  - Map clusters to cell types                 │
│  - Based on marker gene expression            │
│  Output: adata_annotated.h5ad                 │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  VISUALIZATION                                │
│  - Generate UMAP plots                        │
│  - Create dotplots                            │
│  - QC metric overlays                         │
│  Output: Multiple PNG files                   │
└───────────────────────────────────────────────┘
    ↓
┌───────────────────────────────────────────────┐
│  GENERATE_REPORT                              │
│  - Compile all results                        │
│  - Create HTML/Markdown summary               │
│  Output: analysis_summary.html                │
└───────────────────────────────────────────────┘
```

## Module Details

### 1. LOAD_AND_QC
**Purpose**: Initial data loading and quality control

**Inputs**:
- 10X directory path
- QC thresholds (min_genes, max_genes, max_mt_percent)

**Processing**:
1. Load matrix.mtx, genes.tsv, barcodes.tsv
2. Calculate QC metrics (genes per cell, counts, MT%)
3. Filter low-quality cells
4. Filter lowly expressed genes

**Outputs**:
- `adata_qc.h5ad`: Filtered AnnData object
- `qc_metrics.csv`: Before/after statistics
- `qc_violin.png`: QC metric distributions

**Key Parameters**:
- `min_genes`: Minimum genes per cell (default: 200)
- `min_cells`: Minimum cells per gene (default: 3)
- `max_genes`: Maximum genes per cell (default: 2500)
- `max_mt_percent`: Maximum mitochondrial % (default: 5)

---

### 2. NORMALIZE_AND_HVG
**Purpose**: Normalization and feature selection

**Processing**:
1. Normalize total counts to 10,000 per cell
2. Log1p transform
3. Identify highly variable genes (HVGs)
4. Regress out total counts and MT%
5. Scale data (max value = 10)

**Outputs**:
- `adata_normalized.h5ad`: Normalized and scaled data
- `hvg_plot.png`: Highly variable gene plot

**Key Parameters**:
- `n_top_genes`: Number of HVGs to select (default: 2000)

---

### 3. DIMENSIONALITY_REDUCTION
**Purpose**: Reduce dimensions for visualization and clustering

**Processing**:
1. Principal Component Analysis (PCA)
2. Construct k-nearest neighbor graph
3. Compute UMAP embedding

**Outputs**:
- `adata_dimred.h5ad`: Data with PCA and UMAP
- `pca_variance.png`: Variance explained plot

**Key Parameters**:
- `n_pcs`: Number of PCs to compute (default: 50)
- `n_neighbors`: Neighbors for graph (default: 10)
- `n_pcs_umap`: PCs used for UMAP (default: 40)

---

### 4. CLUSTERING
**Purpose**: Group similar cells into clusters

**Processing**:
1. Leiden community detection algorithm
2. Assigns cluster labels to each cell

**Outputs**:
- `adata_clustered.h5ad`: Data with cluster assignments
- `cluster_sizes.csv`: Number of cells per cluster

**Key Parameters**:
- `leiden_resolution`: Clustering resolution (default: 1.0)
  - Higher values = more clusters
  - Lower values = fewer clusters

---

### 5. MARKER_GENES
**Purpose**: Identify genes that define each cluster

**Processing**:
1. Wilcoxon rank-sum test for each cluster vs rest
2. Rank genes by fold change and p-value
3. Extract top markers per cluster

**Outputs**:
- `adata_with_markers.h5ad`: Data with marker results
- `marker_genes.csv`: Top markers per cluster
- `marker_genes_heatmap.png`: Heatmap visualization

**Statistics Generated**:
- Log2 fold change
- P-values
- Adjusted p-values (FDR)
- Gene names

---

### 6. ANNOTATE_CELL_TYPES
**Purpose**: Assign biological cell type labels

**Processing**:
1. Map cluster IDs to cell type names
2. Based on known marker gene expression
3. Can be customized for specific tissues

**Outputs**:
- `adata_annotated.h5ad`: Data with cell type labels
- `cell_type_counts.csv`: Cells per cell type

**Customization**:
Edit `modules/annotate.nf` to define your own cell type mappings:
```python
cell_type_map = {
    '0': 'Your Cell Type Name',
    '1': 'Another Cell Type',
    # etc.
}
```

---

### 7. VISUALIZATION
**Purpose**: Generate publication-quality plots

**Plots Generated**:
- `umap_clusters.png`: UMAP colored by cluster
- `umap_cell_types.png`: UMAP colored by cell type
- `umap_qc_metrics.png`: UMAP colored by QC metrics
- `marker_dotplot.png`: Canonical marker expression

**Customization**:
- Adjust figure sizes in module code
- Add additional marker genes
- Change color schemes

---

### 8. GENERATE_REPORT
**Purpose**: Summarize entire analysis

**Report Contents**:
- Dataset statistics
- QC metrics (before/after)
- Cluster/cell type distributions
- Top marker genes per cluster
- Analysis parameters used

**Outputs**:
- `analysis_summary.html`: Interactive HTML report
- `analysis_summary.md`: Markdown version

---

## Pipeline Execution Flow

### Sequential Dependencies
Each step depends on the previous:
```
LOAD_AND_QC → NORMALIZE_AND_HVG → DIMENSIONALITY_REDUCTION → 
CLUSTERING → MARKER_GENES → ANNOTATE_CELL_TYPES → 
VISUALIZATION → GENERATE_REPORT
```

### Resource Allocation
- QC, normalization: 2 CPUs, 8 GB RAM
- Dimensionality reduction: 4 CPUs, 12 GB RAM
- Clustering: 4 CPUs, 12 GB RAM
- Marker genes: 4 CPUs, 12 GB RAM
- Annotation: 2 CPUs, 8 GB RAM
- Visualization: 2 CPUs, 8 GB RAM

### Error Handling
- Each process retries up to 2 times on failure
- Work directory preserves intermediate files
- Use `-resume` to continue from checkpoint

---

## Data Flow

### Input Format
```
input_directory/
├── matrix.mtx       # Sparse expression matrix
├── genes.tsv        # Gene annotations (2 columns)
└── barcodes.tsv     # Cell barcodes (1 column)
```

### Intermediate Files
All `.h5ad` files contain:
- Expression matrix
- Cell metadata (observations)
- Gene metadata (variables)
- Embeddings (PCA, UMAP)
- Cluster assignments
- QC metrics

### Final Output
Complete processed AnnData object with:
- Raw counts (optional)
- Normalized, log-transformed counts
- Scaled data
- PCA coordinates
- UMAP coordinates
- Cluster labels
- Cell type annotations
- Marker gene results
- QC metrics

---

## Extending the Pipeline

### Adding a New Module

1. Create `modules/your_module.nf`:
```groovy
process YOUR_MODULE {
    tag "Description"
    publishDir "${params.outdir}/XX_your_step", mode: 'copy'
    
    input:
    path(input_file)
    
    output:
    path("output_file"), emit: result
    
    script:
    """
    #!/usr/bin/env python3
    # Your Python code here
    """
}
```

2. Import in `main.nf`:
```groovy
include { YOUR_MODULE } from './modules/your_module.nf'
```

3. Add to workflow:
```groovy
YOUR_MODULE(previous_output)
```

### Common Extensions

**Trajectory Analysis**:
- Add PAGA or RNA velocity modules
- Insert after clustering

**Differential Expression**:
- Compare conditions or cell types
- Add after marker gene identification

**Integration**:
- Combine multiple samples/batches
- Add before or after normalization

**Cell Cycle Scoring**:
- Score cells for cell cycle phase
- Add after normalization

**Gene Set Enrichment**:
- Test pathways in marker genes
- Add after marker identification

---

## Best Practices

1. **Always use `-resume`** when re-running
2. **Check QC metrics** before proceeding
3. **Validate cluster resolution** by testing multiple values
4. **Customize cell type annotations** based on your data
5. **Review marker genes** before final annotation
6. **Keep raw data separate** from pipeline outputs
7. **Version control** your custom configurations

---

## Performance Optimization

### For Large Datasets (>50K cells)
```bash
nextflow run main.nf \
    --max_cpus 16 \
    --max_memory '64.GB' \
    --input data \
    --outdir results
```

### For Quick Testing
```bash
# Reduce resolution and features
nextflow run main.nf \
    --n_top_genes 1000 \
    --n_pcs 30 \
    --input data \
    --outdir test_results
```

### For HPC Clusters
Use SLURM profile with appropriate queue settings in `nextflow.config`.
