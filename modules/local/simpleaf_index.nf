process SIMPLEAF_INDEX {
    tag "reference index"
    label 'process_high'
    stageInMode 'copy'

    input:
    path fasta
    path gtf

    output:
    path "index_dir", emit: index

    when:
    params.run_mode in ['all', 'index_only']

    script:
    """
    simpleaf index \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        --output index_dir \\
        --threads ${task.cpus}

    if [ ! -d "index_dir" ]; then
        echo "Error: simpleaf index failed to create output directory"
        exit 1
    fi
    """
}
