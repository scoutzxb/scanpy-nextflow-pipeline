process VISUALIZATION {
    tag "Generate visualizations"
    publishDir "${params.outdir}/07_visualization", mode: 'copy'
    
    input:
    path(adata)
    val(annotated)
    
    output:
    path("*.png"), emit: plots
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white')
    
    # Load data
    adata = sc.read_h5ad('${adata}')
    
    # UMAP by cluster
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='leiden', ax=ax, title='UMAP - Clusters', 
               show=False, legend_loc='right margin')
    plt.tight_layout()
    plt.savefig('umap_clusters.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # UMAP by cell type (if annotated)
    if ${annotated} and 'cell_type' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.umap(adata, color='cell_type', ax=ax, title='UMAP - Cell Types',
                   show=False, legend_loc='right margin')
        plt.tight_layout()
        plt.savefig('umap_cell_types.png', dpi=150, bbox_inches='tight')
        plt.close()
    
    # UMAP by QC metrics
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[0], 
               title='Genes per Cell', show=False)
    sc.pl.umap(adata, color='total_counts', ax=axes[1],
               title='Total Counts', show=False)
    plt.tight_layout()
    plt.savefig('umap_qc_metrics.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Marker gene dotplot
    if ${annotated} and 'cell_type' in adata.obs.columns:
        marker_genes = {
            'T cells': ['CD3D', 'CD3E'],
            'NK cells': ['NKG7', 'GZMA'],
            'B cells': ['CD79A', 'MS4A1'],
            'Monocytes': ['CD14', 'LYZ', 'FCGR3A'],
            'DC': ['FCER1A', 'CST3']
        }
        
        # Get available markers
        available_markers = {}
        for cell_type, genes in marker_genes.items():
            available = [g for g in genes if g in adata.var_names]
            if available:
                available_markers[cell_type] = available
        
        if available_markers:
            sc.pl.dotplot(adata, var_names=available_markers, 
                         groupby='cell_type', dendrogram=True,
                         show=False, figsize=(14, 6))
            plt.savefig('marker_dotplot.png', dpi=150, bbox_inches='tight')
            plt.close()
    
    print("Visualization complete")
    """
}
