#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Single-cell RNA-seq Analysis Pipeline using Scanpy
========================================================================================
    Author: Xiaobin Zheng
    Description: Automated pipeline for scRNA-seq QC, normalization, clustering, 
                 and cell type identification
========================================================================================
*/

// Print pipeline header
def printHeader() {
    log.info """\
    ===================================================
     S C A N P Y   P I P E L I N E
    ===================================================
     Input directory : ${params.input}
     Output directory: ${params.outdir}
     Min genes/cell  : ${params.min_genes}
     Max genes/cell  : ${params.max_genes}
     Max MT%         : ${params.max_mt_percent}
     Leiden res.     : ${params.leiden_resolution}
     Annotate cells  : ${params.annotate_celltypes}
    ===================================================
     """.stripIndent()
}

// Validate parameters
def validateParams() {
    if (!params.input) {
        error "Please provide an input directory with --input"
    }
    if (!file(params.input).exists()) {
        error "Input directory does not exist: ${params.input}"
    }
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { LOAD_AND_QC } from './modules/load_qc.nf'
include { NORMALIZE_AND_HVG } from './modules/normalize.nf'
include { DIMENSIONALITY_REDUCTION } from './modules/dimensionality_reduction.nf'
include { CLUSTERING } from './modules/clustering.nf'
include { MARKER_GENES } from './modules/marker_genes.nf'
include { ANNOTATE_CELL_TYPES } from './modules/annotate.nf'
include { VISUALIZATION } from './modules/visualization.nf'
include { GENERATE_REPORT } from './modules/report.nf'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    printHeader()
    validateParams()
    
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Step 1: Load data and perform QC
    LOAD_AND_QC(
        input_ch,
        params.min_genes,
        params.min_cells,
        params.max_genes,
        params.max_mt_percent
    )
    
    // Step 2: Normalize and identify highly variable genes
    NORMALIZE_AND_HVG(
        LOAD_AND_QC.out.adata,
        params.n_top_genes
    )
    
    // Step 3: Dimensionality reduction (PCA and UMAP)
    DIMENSIONALITY_REDUCTION(
        NORMALIZE_AND_HVG.out.adata,
        params.n_pcs,
        params.n_neighbors,
        params.n_pcs_umap
    )
    
    // Step 4: Clustering
    CLUSTERING(
        DIMENSIONALITY_REDUCTION.out.adata,
        params.leiden_resolution
    )
    
    // Step 5: Identify marker genes
    MARKER_GENES(
        CLUSTERING.out.adata
    )
    
    // Step 6: Annotate cell types (if enabled)
    if (params.annotate_celltypes) {
        ANNOTATE_CELL_TYPES(
            MARKER_GENES.out.adata
        )
        final_adata = ANNOTATE_CELL_TYPES.out.adata
    } else {
        final_adata = MARKER_GENES.out.adata
    }
    
    // Step 7: Generate visualizations
    VISUALIZATION(
        final_adata,
        params.annotate_celltypes
    )
    
    // Step 8: Generate summary report
    GENERATE_REPORT(
        final_adata,
        LOAD_AND_QC.out.qc_metrics,
        MARKER_GENES.out.markers,
        VISUALIZATION.out.plots
    )
}

workflow.onComplete {
    log.info """\
    ===================================================
     Pipeline completed at: ${workflow.complete}
     Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Execution duration: ${workflow.duration}
     Output directory: ${params.outdir}
    ===================================================
    """.stripIndent()
}
