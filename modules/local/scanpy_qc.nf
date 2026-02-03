process SCANPY_QC {
    tag "${sample_id}"
    label 'process_high'
    stageInMode 'copy'
    publishDir "${params.outdir}/scanpy/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(quants_h5ad), val(meta_map)
    path qc_config

    output:
    tuple val(sample_id), path("${sample_id}_qc.h5ad"), val(meta_map), emit: filtered_adata
    path "*.png", emit: plots
    path "${sample_id}_qc_stats.txt", emit: stats

    when:
    params.run_mode in ['all', 'qc_only']

    script:
    meta_json = groovy.json.JsonOutput.toJson(meta_map)
    """
    python3 ${projectDir}/bin/scanpy_qc.py \\
        --input ${quants_h5ad} \\
        --output ${sample_id}_qc.h5ad \\
        --sample-id ${sample_id} \\
        --config ${qc_config} \\
        --metadata '${meta_json}' \\
        --stats-file ${sample_id}_qc_stats.txt \\
        --plots-dir .
    """
}
