#!/usr/bin/env nextflow

/*
================================================================================
Nextflow DSL2 scRNA-seq Pipeline
================================================================================

Modular pipeline for processing 10x scRNA-seq data:
- Quantification with simpleaf/alevin-fry
- Per-sample QC with Scanpy
- Batch-corrected integration with BBKNN

================================================================================
*/

nextflow.enable.dsl = 2

// Include modules
include { PARSE_SAMPLESHEET } from './modules/local/parse_samplesheet'
include { SIMPLEAF_INDEX } from './modules/local/simpleaf_index'
include { SIMPLEAF_QUANT } from './modules/local/simpleaf_quant'
include { SCANPY_QC } from './modules/local/scanpy_qc'
include { SCANPY_BBKNN } from './modules/local/scanpy_bbknn'


// Main workflow
workflow {
    // Validate inputs
    if (!params.samplesheet) {
        error "Error: --samplesheet parameter is required"
    }

    samplesheet_file = file(params.samplesheet)
    if (!samplesheet_file.exists()) {
        error "Error: Samplesheet file not found: ${params.samplesheet}"
    }

    // Parse samplesheet
    ch_parsed = PARSE_SAMPLESHEET(samplesheet_file)

    // Convert parsed TSV to channel tuples
    ch_samples = ch_parsed.samples
        .splitCsv(header: false, sep: '\t')
        .map { row ->
            [row[0], file(row[1]), file(row[2]), row[3], row[4]]
        }

    // Conditional: Build index
    if (params.run_mode in ['all', 'index_only']) {
        if (!params.fasta) {
            error "Error: --fasta parameter is required for indexing"
        }
        if (!params.annotation) {
            error "Error: --annotation parameter is required for indexing"
        }

        fasta_file = file(params.fasta)
        gtf_file = file(params.annotation)

        if (!fasta_file.exists()) {
            error "Error: FASTA file not found: ${params.fasta}"
        }
        if (!gtf_file.exists()) {
            error "Error: GTF file not found: ${params.annotation}"
        }

        ch_index = SIMPLEAF_INDEX(fasta_file, gtf_file)
    }

    // Conditional: Quantification
    if (params.run_mode in ['all', 'quant_only']) {
        ch_quants = SIMPLEAF_QUANT(ch_index.index, ch_samples)
    } else if (params.run_mode in ['qc_only', 'integration_only']) {
        // For QC/integration-only modes, assume quants exist in work directory
        ch_quants = ch_samples.map { sample_id, fastq_r1, fastq_r2, chemistry, meta_map ->
            [sample_id, file("${sample_id}_quants.h5ad"), meta_map]
        }
    }

    // Conditional: Per-sample QC
    if (params.run_mode in ['all', 'qc_only']) {
        qc_config = file(params.qc_config)
        if (!qc_config.exists()) {
            error "Error: QC config file not found: ${params.qc_config}"
        }

        ch_qc = SCANPY_QC(ch_quants, qc_config)
    } else if (params.run_mode in ['integration_only']) {
        // For integration-only, assume QC files exist
        ch_qc = ch_samples.map { sample_id, fastq_r1, fastq_r2, chemistry, meta_map ->
            [sample_id, file("${sample_id}_qc.h5ad"), meta_map]
        }
    }

    // Conditional: Integration
    if (params.run_mode in ['all', 'integration_only']) {
        bbknn_config = file(params.bbknn_config)
        if (!bbknn_config.exists()) {
            error "Error: BBKNN config file not found: ${params.bbknn_config}"
        }

        // Collect all filtered h5ad files
        ch_merged = ch_qc.filtered_adata
            .map { sample_id, h5ad, meta_map -> h5ad }
            .collect()

        ch_integrated = SCANPY_BBKNN(ch_merged, bbknn_config)
    }
}


workflow.onComplete {
    log.info """
    ================================================================================
    Pipeline execution summary
    ================================================================================
    Command line: ${workflow.commandLine}
    Working directory: ${workflow.workDir}
    Output directory: ${params.outdir}
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    ================================================================================
    """
}


workflow.onError {
    log.error """
    ================================================================================
    Pipeline execution error
    ================================================================================
    Error message: ${workflow.errorMessage}
    ================================================================================
    """
}
