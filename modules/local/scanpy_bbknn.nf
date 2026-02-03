process SCANPY_BBKNN {
    tag "integration"
    label 'process_high'
    stageInMode 'copy'
    publishDir "${params.outdir}/integration/bbknn", mode: 'copy'

    input:
    path qc_files
    path bbknn_config

    output:
    path "merged_bbknn.h5ad", emit: merged_adata
    path "*.png", emit: plots

    when:
    params.run_mode in ['all', 'integration_only']

    script:
    """
    python3 ${projectDir}/bin/scanpy_bbknn_integration.py \\
        --input-files ${qc_files.join(' ')} \\
        --output merged_bbknn.h5ad \\
        --config ${bbknn_config} \\
        --plots-dir .
    """
}
