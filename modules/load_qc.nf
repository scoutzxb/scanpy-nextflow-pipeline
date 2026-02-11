process LOAD_AND_QC {
    tag "QC filtering"
    publishDir "${params.outdir}/01_qc", mode: 'copy'
    
    input:
    path(input_dir)
    val(min_genes)
    val(min_cells)
    val(max_genes)
    val(max_mt_percent)
    
    output:
    path("adata_qc.h5ad"), emit: adata
    path("qc_metrics.csv"), emit: qc_metrics
    path("qc_violin.png"), emit: qc_plot
    
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
    print("Loading 10X data...")
    adata = sc.read_10x_mtx('${input_dir}', var_names='gene_symbols', cache=True)
    
    # Initial stats
    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, 
                                log1p=False, inplace=True)
    
    # Filter
    print("Filtering cells...")
    sc.pp.filter_cells(adata, min_genes=${min_genes})
    sc.pp.filter_genes(adata, min_cells=${min_cells})
    
    adata = adata[adata.obs.n_genes_by_counts < ${max_genes}, :]
    adata = adata[adata.obs.pct_counts_mt < ${max_mt_percent}, :]
    
    # Save QC metrics
    qc_df = pd.DataFrame({
        'metric': ['initial_cells', 'initial_genes', 'filtered_cells', 'filtered_genes'],
        'value': [n_cells_initial, n_genes_initial, adata.n_obs, adata.n_vars]
    })
    qc_df.to_csv('qc_metrics.csv', index=False)
    
    # QC violin plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, ax=axes, show=False)
    plt.tight_layout()
    plt.savefig('qc_violin.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save
    adata.write('adata_qc.h5ad')
    print(f"QC complete: {adata.n_obs} cells, {adata.n_vars} genes")
    """
}
