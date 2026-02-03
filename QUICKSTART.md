# Quick Start Guide

## 5-Minute Setup

### Prerequisites
- Nextflow installed: `nextflow -version`
- Docker/Singularity available: `docker --version` or `singularity --version`
- FASTQ files in paired R1/R2 format
- Reference genome (FASTA) and annotation (GTF)

### Step 1: Create Samplesheet

Create `samplesheet.tsv`:
```tsv
sample_id	fastq_dir	chemistry	batch	donor	condition
sample1	/data/fastq/sample1	10xv3	batch1	donor1	control
sample2	/data/fastq/sample2	10xv3	batch1	donor2	treated
sample3	/data/fastq/sample3	10xv3	batch2	donor1	control
```

**Important:**
- Each row is one sample
- `fastq_dir` must contain `*_R1*.fastq.gz` and `*_R2*.fastq.gz` files
- Additional columns (batch, donor, condition) are optional metadata

Validate samplesheet:
```bash
python bin/check_samplesheet.py samplesheet.tsv --check-files
```

### Step 2: Prepare Reference Files

Ensure you have:
- `reference.fa.gz`: Reference genome FASTA (can be gzipped)
- `reference.gtf.gz`: Gene annotation GTF (can be gzipped)

### Step 3: Run Full Pipeline

```bash
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  --outdir results \
  -profile docker
```

**For SLURM clusters:**
```bash
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -profile slurm
```

### Step 4: Monitor Progress

In another terminal:
```bash
tail -f .nextflow.log
```

### Step 5: Check Results

After completion:
```bash
# Per-sample QC
ls -la results/scanpy/qc/sample1/
# → sample1_qc.h5ad (filtered counts)
# → QC plots (PNG images)

# Integrated analysis
ls -la results/integration/bbknn/
# → merged_bbknn.h5ad (batch-corrected counts)
# → UMAP plots
```

## Common Commands

### Run Individual Stages Only

```bash
# Stage 1: Build index
nextflow run main.nf \
  --run_mode index_only \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz

# Stage 2: Quantification only
nextflow run main.nf \
  --run_mode quant_only \
  --samplesheet samplesheet.tsv

# Stage 3: QC only
nextflow run main.nf \
  --run_mode qc_only \
  --samplesheet samplesheet.tsv

# Stage 4: Integration only
nextflow run main.nf \
  --run_mode integration_only \
  --samplesheet samplesheet.tsv
```

### Customize QC Parameters

```bash
# Stricter filtering
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --min_genes_per_cell 300 \
  --max_pct_mt 15 \
  --n_hvgs 2000

# Disable doublet detection (faster)
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --run_doublets false
```

### Customize Integration Parameters

```bash
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --bbknn_neighbors_within_batch 5 \
  --leiden_resolution 0.8
```

### Resume Failed Run

```bash
# Nextflow automatically caches completed steps
# Rerun from where it failed:
nextflow run main.nf -resume
```

### Generate Reports

```bash
# HTML timeline
open results/reports/execution_timeline.html

# Trace file
cat results/reports/execution_trace.txt
```

## Output Files

After running, you'll have:

### Per-Sample QC
```
results/scanpy/qc/{sample_id}/
├── {sample_id}_qc.h5ad                   # Main output: filtered counts
├── {sample_id}_prefilter_histograms.png  # Before filtering
├── {sample_id}_postfilter_histograms.png # After filtering
├── {sample_id}_pca_variance.png          # PCA cumulative variance
└── {sample_id}_qc_stats.txt              # Statistics summary
```

**Open filtered counts in Python:**
```python
import scanpy as sc
adata = sc.read_h5ad('results/scanpy/qc/sample1/sample1_qc.h5ad')
print(adata)  # ~5k-20k cells × ~20k genes
adata.obs.head()  # Cell metadata including batch, QC metrics
adata.var.head()  # Gene metadata including highly_variable flags
```

### Integrated Analysis
```
results/integration/bbknn/
├── merged_bbknn.h5ad                     # Main output: batch-corrected counts
├── umap_by_batch.png                     # UMAP colored by batch
├── umap_by_leiden.png                    # UMAP with cluster assignments
├── umap_by_sample_id.png                 # UMAP colored by sample
└── umap_by_*.png                         # Additional metadata plots
```

**Open integrated counts in Python:**
```python
import scanpy as sc
adata = sc.read_h5ad('results/integration/bbknn/merged_bbknn.h5ad')
print(adata)  # All samples merged
adata.obs[['batch', 'leiden']].head()  # Batch and cluster assignments
# .obsm['X_umap'] - UMAP coordinates
# .obsm['X_pca'] - PCA coordinates
# .obs['pct_counts_mt'] - QC metric
```

## Troubleshooting

### Issue: "fastq files not found"
**Check:** Ensure your fastq_dir contains properly named files:
```bash
ls /your/data/fastq/sample1/*fastq.gz
# Should show: sample1_R1_*.fastq.gz and sample1_R2_*.fastq.gz
```

### Issue: "Out of memory"
**Solution:** Increase memory and run fewer parallel tasks:
```bash
nextflow run main.nf \
  --max_memory 64.GB \
  --max_cpus 8
```

### Issue: "simpleaf not found in container"
**Solution:** Pre-pull the Docker image:
```bash
docker pull quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0
```

### Issue: "bbknn import error"
**Solution:** The default Scanpy container might not have bbknn. Use:
```bash
nextflow run main.nf \
  --scanpy_image quay.io/biocontainers/scanpy:1.9.6--pyhdfd78af_0
```

Or build custom image (see INSTALLATION.md).

## Next Steps

After running the pipeline:

1. **Explore QC results** in your favorite analysis tool:
   ```python
   import scanpy as sc
   import matplotlib.pyplot as plt

   adata = sc.read_h5ad('results/integration/bbknn/merged_bbknn.h5ad')
   sc.pl.umap(adata, color='leiden', show=True)
   ```

2. **Run downstream analysis:**
   - Differential expression: `sc.tl.rank_genes_groups()`
   - Cell type annotation: Manual or using marker genes
   - Trajectory inference: `scvelo` or `palantir`
   - Protein activity: `decoupler`, `progeny`

3. **Customize further:**
   - Edit QC parameters for stricter filtering
   - Adjust BBKNN neighbors for different batch integration strength
   - Change Leiden resolution for different clustering granularity

## Key Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `--min_genes_per_cell` | 200 | Minimum genes expressed per cell |
| `--max_pct_mt` | 20 | Max mitochondrial % per cell |
| `--min_cells_per_gene` | 3 | Min cells expressing each gene |
| `--n_hvgs` | 3000 | Number of highly variable genes |
| `--bbknn_neighbors_within_batch` | 3 | BBKNN neighborhood size within batch |
| `--leiden_resolution` | 1.0 | Leiden clustering resolution |
| `--run_doublets` | true | Run scDblFinder doublet detection |

## Getting Help

```bash
# Show all parameters
nextflow run main.nf --help

# View configuration
nextflow config main.nf

# Debug mode
nextflow run main.nf -debug
```

For detailed documentation, see:
- README.md: Full pipeline documentation
- INSTALLATION.md: Installation and environment setup
- nextflow.config: All configurable parameters
