process PARSE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_single'

    input:
    path samplesheet

    output:
    tuple val(sample_id), path(fastq_r1), path(fastq_r2), val(chemistry), val(meta_map), emit: samples

    script:
    """
    python3 << 'EOF'
import sys
import csv
import json
import re
from pathlib import Path

# Read and validate samplesheet
samplesheet_file = "${samplesheet}"

try:
    with open(samplesheet_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        # Check required columns
        if not reader.fieldnames:
            print("Error: Empty samplesheet", file=sys.stderr)
            sys.exit(1)

        required_cols = {'sample_id', 'fastq_dir', 'chemistry'}
        missing = required_cols - set(reader.fieldnames)

        if missing:
            print(f"Error: Missing required columns: {missing}", file=sys.stderr)
            sys.exit(1)

        # Process each sample
        for row in reader:
            sample_id = row['sample_id'].strip()
            fastq_dir = Path(row['fastq_dir'].strip())
            chemistry = row['chemistry'].strip()

            # Validate inputs
            if not sample_id:
                print("Error: sample_id cannot be empty", file=sys.stderr)
                sys.exit(1)

            if not fastq_dir.exists():
                print(f"Error: fastq_dir does not exist: {fastq_dir}", file=sys.stderr)
                sys.exit(1)

            if not chemistry:
                print("Error: chemistry cannot be empty", file=sys.stderr)
                sys.exit(1)

            # Find R1 and R2 files
            r1_pattern = re.compile(r'.*_R1(_001)?\\.fastq\\.gz\$')
            r2_pattern = re.compile(r'.*_R2(_001)?\\.fastq\\.gz\$')

            r1_files = sorted([f for f in fastq_dir.glob('*.fastq.gz') if r1_pattern.match(f.name)])
            r2_files = sorted([f for f in fastq_dir.glob('*.fastq.gz') if r2_pattern.match(f.name)])

            if not r1_files or not r2_files:
                print(f"Error: Could not find paired FASTQ files in {fastq_dir}", file=sys.stderr)
                sys.exit(1)

            if len(r1_files) != len(r2_files):
                print(f"Error: Mismatched R1/R2 file counts in {fastq_dir}", file=sys.stderr)
                sys.exit(1)

            # Concatenate multiple R1/R2 files if needed
            fastq_r1 = r1_files[0]
            fastq_r2 = r2_files[0]

            # Build metadata map from extra columns
            meta_map = {}
            for col in reader.fieldnames:
                if col not in required_cols:
                    value = row[col].strip()
                    if value:
                        meta_map[col] = value

            # Set default batch if not provided
            if 'batch' not in meta_map:
                meta_map['batch'] = sample_id

            # Output in Nextflow format
            print(f"{sample_id}\\t{fastq_r1}\\t{fastq_r2}\\t{chemistry}\\t{json.dumps(meta_map)}")

except Exception as e:
    print(f"Error: {str(e)}", file=sys.stderr)
    sys.exit(1)

EOF
    """
}
