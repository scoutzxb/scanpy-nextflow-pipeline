process DIMENSIONALITY_REDUCTION {
    tag "PCA and UMAP"
    publishDir "${params.outdir}/03_dimred", mode: 'copy'
    
    input:
    path(adata_normalized)
    val(n_pcs)
    val(n_neighbors)
    val(n_pcs_umap)
    
    output:
    path("adata_dimred.h5ad"), emit: adata
    path("pca_variance.png"), emit: pca_plot
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    sc.settings.verbosity = 3
    
    # Load data
    adata = sc.read_h5ad('${adata_normalized}')
    
    # PCA
    print("Computing PCA...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=${n_pcs})
    
    # Plot variance explained
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    plt.savefig('pca_variance.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Compute neighborhood graph
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=${n_neighbors}, n_pcs=${n_pcs_umap})
    
    # UMAP
    print("Computing UMAP...")
    sc.tl.umap(adata)
    
    # Save
    adata.write('adata_dimred.h5ad')
    print("Dimensionality reduction complete")
    """
}
