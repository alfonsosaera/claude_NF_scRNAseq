# Installation Guide

## Prerequisites

### System Requirements
- Linux/macOS (Unix-like environment)
- 8+ CPU cores (16+ recommended)
- 32+ GB RAM (64+ GB for large datasets)
- Docker or Singularity for containerization

### Software Requirements
- **Nextflow** >= 23.0.0
- **Java** (required by Nextflow)
- **Docker** >= 19.0 OR **Singularity** >= 3.0
- **Python** 3.8+ (for validation scripts)

## Installation Steps

### 1. Install Java

```bash
# Ubuntu/Debian
sudo apt-get install openjdk-11-jdk

# macOS
brew install openjdk@11

# Verify installation
java -version
```

### 2. Install Nextflow

```bash
# Install Nextflow (requires Java)
curl -s https://get.nextflow.io | bash

# Make executable and move to PATH
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

Or use conda:
```bash
conda create -n nextflow -c bioconda nextflow
conda activate nextflow
```

### 3. Install Docker (for local execution)

```bash
# Ubuntu/Debian
sudo apt-get install docker.io

# Start Docker daemon
sudo systemctl start docker
sudo systemctl enable docker

# Add user to docker group (optional, to avoid sudo)
sudo usermod -aG docker $USER
newgrp docker

# Verify installation
docker --version
docker run hello-world
```

Or use macOS:
```bash
# Install Docker Desktop from https://www.docker.com/products/docker-desktop
docker --version
```

### 4. Install Singularity (alternative to Docker)

```bash
# Ubuntu/Debian
sudo apt-get install singularity-container

# Verify installation
singularity --version
```

Or from source:
```bash
# See: https://sylabs.io/guides/3.0/user-guide/installation.html
```

### 5. Clone/Download Pipeline

```bash
# Clone from GitHub (if available)
git clone https://github.com/your-org/claude_NF_scRNAseq
cd claude_NF_scRNAseq

# Or download as ZIP
wget https://github.com/your-org/claude_NF_scRNAseq/archive/refs/heads/main.zip
unzip main.zip
cd claude_NF_scRNAseq-main
```

### 6. Verify Installation

```bash
# Check Nextflow can find Docker images
nextflow info

# Verify pipeline configuration
nextflow config main.nf | head -20

# Test samplesheet validator
python bin/check_samplesheet.py examples/samplesheet.tsv
```

## Optional: Build Custom Scanpy Container

The pipeline uses the biocontainers Scanpy image by default. For optimal performance with all features (especially doublet detection), build a custom image:

### Option A: Using Provided Dockerfile

```bash
# Build custom image
docker build -f Dockerfile.scanpy -t myorg/scanpy:custom .

# Or with Singularity
singularity build scanpy_custom.sif Dockerfile.scanpy

# Use in pipeline
nextflow run main.nf --scanpy_image myorg/scanpy:custom ...
```

### Option B: Manual Installation

If you prefer to install dependencies in a running container:

```bash
# Start Scanpy container
docker run -it quay.io/biocontainers/scanpy:1.9.6--pyhdfd78af_0 /bin/bash

# Inside container, install additional packages
pip install bbknn rpy2
apt-get update && apt-get install -y r-base
R -e "install.packages('BiocManager'); BiocManager::install('scDblFinder')"

# Exit and commit container
# docker commit <container_id> myorg/scanpy:custom
```

## Testing Installation

### Minimal Test

Run the help command:
```bash
nextflow run main.nf --help
```

Expected output: Pipeline parameters and usage information

### Quick Validation Test

Test with example samplesheet (requires no real data):

```bash
# This will fail gracefully at the quantification step due to missing FASTQ files
# But validates the pipeline structure
nextflow run main.nf \
  --samplesheet examples/samplesheet.tsv \
  --run_mode qc_only \
  -preview
```

### Full Test with Test Data

Prepare small test datasets:

```bash
# Create test data directory
mkdir -p test_data/fastq/{sample1,sample2}

# Generate small test FASTQ files (example using art_illumina or similar)
# For testing without real sequences:
touch test_data/fastq/sample1/sample1_R1.fastq.gz
touch test_data/fastq/sample1/sample1_R2.fastq.gz
touch test_data/fastq/sample2/sample2_R1.fastq.gz
touch test_data/fastq/sample2/sample2_R2.fastq.gz

# Update test samplesheet
cat > test_samplesheet.tsv << 'EOF'
sample_id	fastq_dir	chemistry
sample1	$(pwd)/test_data/fastq/sample1	10xv3
sample2	$(pwd)/test_data/fastq/sample2	10xv3
EOF

# Run test profile
nextflow run main.nf \
  --samplesheet test_samplesheet.tsv \
  -profile test \
  -preview
```

## Configuration for Different Environments

### Local Machine (Docker)

```bash
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -profile docker
```

### Local Machine (Singularity)

```bash
# Create nextflow config
cat > local.config << 'EOF'
singularity.enabled = true
docker.enabled = false
process.container = 'quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0'
EOF

nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -c local.config
```

### SLURM Cluster

```bash
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -profile slurm
```

### HPC with Module System

```bash
# Load required modules (example)
module load java/11
module load nextflow
module load singularity

# Run with Singularity
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta reference.fa.gz \
  --annotation reference.gtf.gz \
  -profile slurm \
  -c hpc_config.nf
```

### Cloud Environments

#### AWS

```bash
# Requires AWS CLI configured
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta s3://my-bucket/reference.fa.gz \
  --annotation s3://my-bucket/reference.gtf.gz \
  -profile aws
```

#### Google Cloud

```bash
# Requires Google Cloud SDK configured
nextflow run main.nf \
  --samplesheet samplesheet.tsv \
  --fasta gs://my-bucket/reference.fa.gz \
  --annotation gs://my-bucket/reference.gtf.gz \
  -profile gcloud
```

## Troubleshooting

### Issue: Docker permission denied

**Solution:**
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Log out and log back in, or:
newgrp docker
```

### Issue: Nextflow command not found

**Solution:**
```bash
# Ensure Nextflow is in PATH
export PATH=$PATH:/path/to/nextflow
# Or move to system path:
sudo mv /path/to/nextflow /usr/local/bin/
```

### Issue: Image pull rate limiting

**Solution:**
Pre-pull Docker images before running:
```bash
docker pull quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0
docker pull quay.io/biocontainers/scanpy:1.9.6--pyhdfd78af_0
```

### Issue: Out of disk space

**Solution:**
```bash
# Check disk usage
df -h

# Clean Nextflow cache
nextflow clean -f -k

# Remove work directory
rm -rf work/
```

### Issue: Python version incompatibility

**Solution:**
```bash
# Check Python version
python --version  # Should be >= 3.8

# Use specific Python version if available
python3.10 bin/check_samplesheet.py samplesheet.tsv
```

## Environment Customization

### Create Nextflow Config for Your Setup

```bash
cat > myenv.config << 'EOF'
// Custom environment configuration

process {
    executor = 'slurm'
    cpus = 8
    memory = 16.GB
    time = 4.h

    withName: SIMPLEAF_QUANT {
        cpus = 16
        memory = 32.GB
    }
}

executor {
    queueSize = 50
    submitRateLimit = '10 sec'
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '/scratch/singularity'
}

params {
    outdir = '/scratch/results'
}
EOF

# Use custom config
nextflow run main.nf -c myenv.config ...
```

### Load Module Environment

```bash
# For systems with module system
module load java/11
module load nextflow/23.10.0
module load singularity/3.8
module load python/3.10

# Then run pipeline
nextflow run main.nf ...
```

## Updating Nextflow

```bash
# Update to latest version
nextflow self-update

# Update to specific version
nextflow self-update 24.10.0
```

## Getting Help

For issues or questions:

1. Check the logs:
```bash
tail -f .nextflow.log
```

2. Run with verbose output:
```bash
nextflow run main.nf ... -debug
```

3. Check Nextflow documentation:
```bash
nextflow help
nextflow run --help
```

4. For tool-specific issues, consult:
   - simpleaf: https://github.com/COMBINE-lab/simpleaf
   - Scanpy: https://scanpy.readthedocs.io
   - Nextflow: https://www.nextflow.io/docs
