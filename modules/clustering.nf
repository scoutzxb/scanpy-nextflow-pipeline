process CLUSTERING {
    tag "Leiden clustering"
    publishDir "${params.outdir}/04_clustering", mode: 'copy'
    
    input:
    path(adata_dimred)
    val(resolution)
    
    output:
    path("adata_clustered.h5ad"), emit: adata
    path("cluster_sizes.csv"), emit: cluster_sizes
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    
    sc.settings.verbosity = 3
    
    # Load data
    adata = sc.read_h5ad('${adata_dimred}')
    
    # Leiden clustering
    print("Performing Leiden clustering...")
    sc.tl.leiden(adata, resolution=${resolution})
    
    # Save cluster sizes
    cluster_df = adata.obs['leiden'].value_counts().reset_index()
    cluster_df.columns = ['cluster', 'n_cells']
    cluster_df = cluster_df.sort_values('cluster')
    cluster_df.to_csv('cluster_sizes.csv', index=False)
    
    print(f"Identified {adata.obs['leiden'].nunique()} clusters")
    print(cluster_df)
    
    # Save
    adata.write('adata_clustered.h5ad')
    """
}
