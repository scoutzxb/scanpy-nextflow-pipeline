process MARKER_GENES {
    tag "Marker gene identification"
    publishDir "${params.outdir}/05_markers", mode: 'copy'
    
    input:
    path(adata_clustered)
    
    output:
    path("adata_with_markers.h5ad"), emit: adata
    path("marker_genes.csv"), emit: markers
    path("marker_genes_heatmap.png"), emit: heatmap
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=150, facecolor='white')
    
    # Load data
    adata = sc.read_h5ad('${adata_clustered}')
    
    # Find marker genes
    print("Identifying marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    # Extract top markers
    n_genes = 10
    markers_list = []
    
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster, key='rank_genes_groups')
        cluster_markers['cluster'] = cluster
        cluster_markers = cluster_markers.head(n_genes)
        markers_list.append(cluster_markers)
    
    markers_df = pd.concat(markers_list)
    markers_df.to_csv('marker_genes.csv', index=False)
    
    # Plot heatmap of top markers
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby='leiden', 
                                     show_gene_labels=True, show=False, 
                                     figsize=(12, 8))
    plt.savefig('marker_genes_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save
    adata.write('adata_with_markers.h5ad')
    print("Marker gene identification complete")
    """
}
