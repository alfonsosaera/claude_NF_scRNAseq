process SIMPLEAF_QUANT {
    tag "${sample_id}"
    label 'process_high'

    input:
    path index_dir
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(chemistry), val(meta_map)

    output:
    tuple val(sample_id), path("${sample_id}_quants.h5ad"), val(meta_map), emit: quants

    when:
    params.run_mode in ['all', 'quant_only']

    script:
    """
    # Export required environment variables
    export ALEVIN_FRY_HOME=.

    # Set maximum number of file descriptors for temp files
    ulimit -n 2048

    # Prepare simpleaf
    simpleaf set-paths

    # Create output directory
    mkdir -p quant_output

    # Run simpleaf quant
    simpleaf quant \\
        --index ${index_dir}/index \\
        --reads1 ${fastq_r1} \\
        --reads2 ${fastq_r2} \\
        --chemistry ${chemistry} \\
        --t2g-map ${index_dir}/ref/t2g_3col.tsv \\
        --threads ${task.cpus} \\
        --resolution cr-like \\
        --unfiltered-pl \\
        --anndata-out \\
        --output quant_output

    # Check if quants.h5ad was created
    if [ ! -f "quant_output/af_quant/alevin/quants.h5ad" ]; then
        echo "Error: simpleaf quant failed to create quants.h5ad"
        exit 1
    fi

    # Copy and rename output
    cp quant_output/af_quant/alevin/quants.h5ad ${sample_id}_quants.h5ad

    # Verify output
    if [ ! -f "${sample_id}_quants.h5ad" ]; then
        echo "Error: Failed to copy quants file"
        exit 1
    fi
    """
}
