# Nextflow DSL2 scRNA-seq Pipeline

> This pipeline was built as an exercise to test Claude Code on building NextFlow pipelines as<br>
> show in this [blog post](https://lovednacodeblog.com/post/2026-02-10-claude-nf-simpleaf/). The single cell analysis code was AI generated and I did not review it.<br>
> Feel free to modify and use this pipeline if you wish.

A modular Nextflow pipeline for processing 10x scRNA-seq data with:
- **Quantification**: simpleaf/alevin-fry for rapid quantification
- **Per-sample QC**: Scanpy-based quality control, filtering, and normalization
- **Batch Integration**: BBKNN for batch-corrected integration

## Features

- ✅ Modular DSL2 design with independent processes
- ✅ Flexible run modes: `all`, `index_only`, `quant_only`, `qc_only`, `integration_only`
- ✅ Metadata propagation through entire pipeline
- ✅ Doublet detection with scDblFinder (optional)
- ✅ Comprehensive QC visualizations
- ✅ Batch-corrected UMAP integration with BBKNN
- ✅ Support for multiple executors (local, SLURM)
- ✅ Configurable parameters via JSON config files

## Quick Start

### 1. Create a Samplesheet

Create a tab-separated file `samplesheet.tsv`:
```
sample_id	fastq_dir	chemistry	batch	donor	condition
sample1	/path/to/fastq/sample1	10xv3	batch1	donor1	control
sample2	/path/to/fastq/sample2	10xv3	batch1	donor2	treated
sample3	/path/to/fastq/sample3	10xv3	batch2	donor1	control
```

**Requirements:**
- `sample_id`: Unique identifier for each sample
- `fastq_dir`: Directory containing paired R1/R2 FASTQ files (must contain `*_R1*.fastq.gz` and `*_R2*.fastq.gz`)
- `chemistry`: 10x chemistry version (10xv2, 10xv3, etc.)
- Additional columns (batch, donor, condition) are optional metadata propagated through the pipeline

### 2. Run the Pipeline

```bash
# Full pipeline: indexing → quantification → QC → integration
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -profile docker

# Or run individual stages
nextflow run main.nf --run_mode index_only --fasta reference.fa.gz --annotation reference.gtf.gz
nextflow run main.nf --run_mode quant_only --samplesheet samplesheet.tsv
nextflow run main.nf --run_mode qc_only --samplesheet samplesheet.tsv
nextflow run main.nf --run_mode integration_only --samplesheet samplesheet.tsv
```

### 3. Outputs

- **Per-sample QC**: `results/scanpy/qc/{sample_id}/{sample_id}_qc.h5ad` + plots
- **Integration**: `results/integration/bbknn/merged_bbknn.h5ad` + UMAP plots

## Full Documentation

### Installation

Requirements:
- Nextflow >= 23.0.0
- Docker or Singularity
- Python 3.8+

### Pipeline Parameters

**Key parameters:**
- `--samplesheet`: Path to input samplesheet (TSV)
- `--fasta`: Reference genome FASTA file
- `--annotation`: Reference annotation GTF file
- `--outdir`: Output directory (default: `./results`)
- `--run_mode`: Pipeline stage to run (default: `all`)

**QC Parameters:**
- `--min_genes_per_cell`: Minimum genes per cell (default: 200)
- `--max_genes_per_cell`: Maximum genes per cell (default: null)
- `--min_counts_per_cell`: Minimum counts per cell (default: 500)
- `--max_pct_mt`: Maximum mitochondrial percentage (default: 20)
- `--run_doublets`: Run doublet detection (default: true)
- `--n_hvgs`: Number of highly variable genes (default: 3000)

**Integration Parameters:**
- `--bbknn_neighbors_within_batch`: BBKNN k-nearest neighbors within batch (default: 3)
- `--leiden_resolution`: Leiden clustering resolution (default: 1.0)

**Container Images:**
- `--simpleaf_image`: simpleaf Docker image (default: `quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0`)
- `--scanpy_image`: Scanpy Docker image (default: `quay.io/biocontainers/scanpy:1.9.6--pyhdfd78af_0`)

### Configuration

Modify parameters via:

1. **Command line:**
```bash
nextflow run main.nf --min_genes_per_cell 150 --max_pct_mt 25
```

2. **Config files:**
   - `assets/qc_config.json`: QC parameters and thresholds
   - `assets/bbknn_config.json`: Integration parameters
   - `conf/test.config`: Test profile with minimal resources

3. **Nextflow profiles:**
```bash
-profile local      # Local executor
-profile slurm      # SLURM cluster
-profile docker     # Enable Docker containers
```

### Output Directory Structure

```
results/
├── reports/
│   ├── execution_trace.txt          # Execution metadata
│   └── execution_timeline.html       # Timeline visualization
├── scanpy/
│   └── qc/
│       └── {sample_id}/
│           ├── {sample_id}_qc.h5ad                    # Filtered AnnData
│           ├── {sample_id}_prefilter_*.png            # Pre-filter plots
│           ├── {sample_id}_postfilter_*.png           # Post-filter plots
│           └── {sample_id}_qc_stats.txt               # QC statistics
└── integration/
    └── bbknn/
        ├── merged_bbknn.h5ad                          # Integrated AnnData
        ├── umap_by_batch.png
        ├── umap_by_leiden.png
        ├── umap_by_sample_id.png
        └── umap_by_*.png                              # Metadata plots
```

### Per-Sample QC (`scanpy_qc.py`)

Performs:
1. Mitochondrial gene flagging
2. QC metrics calculation (n_genes, total_counts, pct_mt)
3. Cell filtering (genes/counts per cell, mitochondrial percentage)
4. Gene filtering (minimum cells per gene)
5. Doublet detection (optional, via scDblFinder)
6. Normalization (total count → log-transform)
7. Highly variable gene selection (Seurat v3)
8. PCA decomposition

**Outputs:**
- `{sample_id}_qc.h5ad`: Normalized, filtered, dimensionality-reduced AnnData
- QC plots: Histograms, violin plots, scatter plots

### Integration (`scanpy_bbknn_integration.py`)

Performs on merged data:
1. Merge all per-sample AnnData objects
2. Recompute HVGs on merged data
3. Scale and PCA
4. BBKNN batch-aware neighbor detection
5. UMAP embedding
6. Leiden clustering

**Outputs:**
- `merged_bbknn.h5ad`: Integrated AnnData with batch correction
- UMAP plots: Colored by batch, leiden, sample_id, and metadata

### Run Modes

- `all`: Complete pipeline (index → quant → QC → integration)
- `index_only`: Build reference index only
- `quant_only`: Run quantification only
- `qc_only`: Run per-sample QC only
- `integration_only`: Run integration only

### Executors

**Local:**
```bash
-profile local
```

**SLURM:**
```bash
-profile slurm
```

Customize executor settings in `nextflow.config` `executor` block.

### Test Mode

Test with example samplesheet and minimal resources:
```bash
nextflow run main.nf -profile test
```

Update `examples/samplesheet.tsv` with valid test data paths before running.

## Troubleshooting

### Container image not found
Pre-pull Docker images:
```bash
docker pull quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0
docker pull quay.io/biocontainers/scanpy:1.9.6--pyhdfd78af_0
```

### bbknn import error
Install bbknn in your Scanpy container or provide a custom image:
```bash
pip install bbknn
```

### scDblFinder not available
Doublet detection is optional and skipped if rpy2/R are unavailable. To enable:
```
# In custom Dockerfile:
RUN apt-get install -y r-base
RUN R -e "install.packages('BiocManager'); BiocManager::install('scDblFinder')"
```

### Out of memory during integration
For large datasets (>20 samples), increase memory allocation:
```bash
nextflow run main.nf --max_memory 64.GB ...
```

### Validate samplesheet

```bash
python bin/check_samplesheet.py samplesheet.tsv --check-files
```

## Project Structure

```
.
├── main.nf                              # Main workflow
├── nextflow.config                      # Pipeline configuration
├── modules/
│   └── local/
│       ├── parse_samplesheet.nf        # Samplesheet parser
│       ├── simpleaf_index.nf           # Reference indexing
│       ├── simpleaf_quant.nf           # Quantification
│       ├── scanpy_qc.nf                # QC module
│       └── scanpy_bbknn.nf             # Integration module
├── bin/
│   ├── scanpy_qc.py                    # QC processing script
│   ├── scanpy_bbknn_integration.py     # Integration script
│   └── check_samplesheet.py            # Samplesheet validator
├── conf/
│   ├── base.config                     # Base process configs
│   └── test.config                     # Test profile
├── assets/
│   ├── qc_config.json                  # QC default parameters
│   └── bbknn_config.json               # Integration parameters
└── examples/
    └── samplesheet.tsv                 # Example samplesheet
```

## Pipeline Architecture

### Module: PARSE_SAMPLESHEET
- **Input:** TSV samplesheet
- **Output:** Channels with (sample_id, fastq_r1, fastq_r2, chemistry, metadata)
- **Function:** Validates format, locates FASTQ files, parses metadata

### Module: SIMPLEAF_INDEX
- **Input:** FASTA, GTF
- **Output:** Indexed reference directory
- **Function:** Builds simpleaf-compatible index (runs once, broadcast to quantification)

### Module: SIMPLEAF_QUANT
- **Input:** Index, per-sample FASTQ
- **Output:** Quantification AnnData (quants.h5ad)
- **Function:** Runs simpleaf quantification with 10x chemistry detection

### Module: SCANPY_QC
- **Input:** Per-sample quants.h5ad
- **Output:** Filtered, normalized AnnData + QC plots
- **Function:** QC filtering, normalization, HVG selection, optional doublet detection

### Module: SCANPY_BBKNN
- **Input:** Collection of per-sample filtered AnnData
- **Output:** Integrated merged AnnData + UMAP plots
- **Function:** Merges samples, applies BBKNN batch correction, Leiden clustering

## Key Features

### Metadata Propagation
Columns in the samplesheet (beyond required: sample_id, fastq_dir, chemistry) are automatically propagated as metadata through:
1. Parsing → stored in channel tuples
2. QC processing → added to AnnData.obs
3. Integration → available in final merged object

### Doublet Detection
If `run_doublets=true` in QC config:
- Uses scDblFinder R package via rpy2 bridge
- Adds `doublet_score` and `doublet_class` to AnnData.obs
- Gracefully skips if dependencies unavailable

### Batch Correction
BBKNN automatically:
- Respects batch identity (from `batch` column or sample_id)
- Ensures k-nearest neighbor graph is balanced across batches
- Prevents batch artifacts while preserving true biological variation

## Performance

- **Quantification**: Scales with FASTQ size; CPU/memory depends on reference size
- **QC**: Runs in parallel per sample; ~5-30 min per 10k cells depending on compute
- **Integration**: Requires loading all samples into memory; ~64GB for 20 × 10k cell samples

## Citation

If you use this pipeline in your research, please cite:

- **Nextflow**: Di Tommaso et al. (2017). Nature Biotechnology
- **simpleaf**: Uman et al. (2023). bioRxiv
- **Scanpy**: Wolf et al. (2018). Genome Biology
- **BBKNN**: Polanski et al. (2020). Nature Communications

## License

MIT
