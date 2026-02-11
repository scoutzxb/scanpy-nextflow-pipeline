process ANNOTATE_CELL_TYPES {
    tag "Cell type annotation"
    publishDir "${params.outdir}/06_annotation", mode: 'copy'
    
    input:
    path(adata_with_markers)
    
    output:
    path("adata_annotated.h5ad"), emit: adata
    path("cell_type_counts.csv"), emit: cell_type_counts
    
    script:
    """
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    
    sc.settings.verbosity = 3
    
    # Load data
    adata = sc.read_h5ad('${adata_with_markers}')
    
    # Define cell type annotations based on marker genes
    # This is a basic annotation - should be customized for your data
    cell_type_map = {
        '0': 'Naive CD4+ T cells',
        '1': 'CD14+ Monocytes',
        '2': 'B cells',
        '3': 'CD8+ T cells',
        '4': 'NK cells',
        '5': 'FCGR3A+ Monocytes',
        '6': 'Dendritic cells',
        '7': 'Platelets'
    }
    
    # Map clusters to cell types
    adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type_map)
    
    # For unmapped clusters, keep cluster number
    adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Cluster_' + adata.obs['leiden'])
    
    # Save cell type counts
    cell_type_df = adata.obs['cell_type'].value_counts().reset_index()
    cell_type_df.columns = ['cell_type', 'n_cells']
    cell_type_df.to_csv('cell_type_counts.csv', index=False)
    
    print("Cell type distribution:")
    print(cell_type_df)
    
    # Save
    adata.write('adata_annotated.h5ad')
    print("Cell type annotation complete")
    """
}
