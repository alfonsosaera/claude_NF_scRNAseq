process SIMPLEAF_INDEX {
    tag "reference index"
    label 'process_high'

    input:
    path fasta
    path gtf

    output:
    path "index/", emit: index

    when:
    params.run_mode in ['all', 'index_only']

    script:
    """
    # Export required environment variables
    export ALEVIN_FRY_HOME=.

    # Set maximum number of file descriptors for temp files
    ulimit -n 2048

    # Prepare simpleaf
    simpleaf set-paths

    # Run index
    simpleaf index \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        --output index \\
        --threads ${task.cpus}

    if [ ! -d "index/index" ]; then
        echo "Error: simpleaf index failed to create index directory"
        exit 1
    fi
    """
}
