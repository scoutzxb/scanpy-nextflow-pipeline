process GENERATE_REPORT {
    tag "Generate summary report"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(adata)
    path(qc_metrics)
    path(markers)
    path(plots)
    
    output:
    path("analysis_summary.html"), emit: report
    path("analysis_summary.md"), emit: markdown
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    from datetime import datetime
    
    # Load data
    adata = sc.read_h5ad('${adata}')
    qc_df = pd.read_csv('${qc_metrics}')
    markers_df = pd.read_csv('${markers}')
    
    # Get stats
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    n_clusters = adata.obs['leiden'].nunique()
    
    # Create markdown report
    md_content = f'''# Single-Cell RNA-seq Analysis Report

## Analysis Summary

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

### Dataset Statistics
- **Total cells after QC:** {n_cells:,}
- **Total genes after QC:** {n_genes:,}
- **Number of clusters:** {n_clusters}

### QC Metrics
'''
    
    for _, row in qc_df.iterrows():
        md_content += f"- **{row['metric']}:** {row['value']:,}\\n"
    
    md_content += '''
### Analysis Pipeline

1. **Data Loading & QC**: Filtered low-quality cells and genes
2. **Normalization**: Total count normalization to 10,000 counts per cell
3. **Feature Selection**: Identified highly variable genes
4. **Dimensionality Reduction**: PCA and UMAP
5. **Clustering**: Leiden algorithm
6. **Marker Gene Discovery**: Wilcoxon rank-sum test
'''
    
    if 'cell_type' in adata.obs.columns:
        md_content += '7. **Cell Type Annotation**: Identified cell types based on marker genes\\n\\n'
        md_content += '### Cell Type Distribution\\n\\n'
        
        cell_type_counts = adata.obs['cell_type'].value_counts()
        md_content += '| Cell Type | Cell Count | Percentage |\\n'
        md_content += '|-----------|-----------|-----------|\\n'
        for ct, count in cell_type_counts.items():
            pct = (count / n_cells) * 100
            md_content += f'| {ct} | {count:,} | {pct:.1f}% |\\n'
    else:
        md_content += '\\n### Cluster Distribution\\n\\n'
        cluster_counts = adata.obs['leiden'].value_counts().sort_index()
        md_content += '| Cluster | Cell Count |\\n'
        md_content += '|---------|-----------|\\n'
        for cluster, count in cluster_counts.items():
            md_content += f'| {cluster} | {count:,} |\\n'
    
    md_content += '''
### Top Marker Genes per Cluster

'''
    
    for cluster in markers_df['cluster'].unique():
        cluster_markers = markers_df[markers_df['cluster'] == cluster].head(5)
        md_content += f"**Cluster {cluster}:** "
        md_content += ", ".join(cluster_markers['names'].tolist()) + "\\n\\n"
    
    md_content += '''
## Visualizations

See the `07_visualization` directory for all plots.

## Files Generated

- **adata_annotated.h5ad**: Final annotated AnnData object
- **marker_genes.csv**: Complete list of marker genes
- **cluster_sizes.csv**: Cell counts per cluster
- **cell_type_counts.csv**: Cell counts per cell type (if annotated)

## Next Steps

1. Validate cell type annotations with biological expertise
2. Perform differential expression analysis between conditions
3. Investigate specific cell populations of interest
4. Integrate with additional datasets
5. Perform trajectory inference if applicable
'''
    
    # Write markdown
    with open('analysis_summary.md', 'w') as f:
        f.write(md_content)
    
    # Convert to HTML (basic)
    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <title>scRNA-seq Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        table {{ border-collapse: collapse; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        code {{ background-color: #f4f4f4; padding: 2px 5px; border-radius: 3px; }}
    </style>
</head>
<body>
    <pre>{md_content}</pre>
</body>
</html>'''
    
    with open('analysis_summary.html', 'w') as f:
        f.write(html_content)
    
    print("Report generation complete")
    """
}
