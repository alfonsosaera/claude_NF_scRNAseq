#!/usr/bin/env python3
"""
Validate scRNA-seq samplesheet format.
Checks required columns and file existence.
"""

import sys
import csv
import re
from pathlib import Path
from typing import List, Tuple


def check_samplesheet(samplesheet_file: str, check_files: bool = False) -> Tuple[bool, List[str]]:
    """
    Validate samplesheet format.

    Args:
        samplesheet_file: Path to TSV samplesheet
        check_files: Whether to verify FASTQ files exist

    Returns:
        (is_valid, error_messages)
    """
    errors = []

    try:
        with open(samplesheet_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')

            if not reader.fieldnames:
                errors.append("Samplesheet is empty")
                return False, errors

            # Check required columns
            required_cols = {'sample_id', 'fastq_dir', 'chemistry'}
            missing = required_cols - set(reader.fieldnames)

            if missing:
                errors.append(f"Missing required columns: {missing}")
                return False, errors

            # Validate each row
            for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is 1)
                sample_id = row['sample_id'].strip()
                fastq_dir = row['fastq_dir'].strip()
                chemistry = row['chemistry'].strip()

                # Validate sample_id
                if not sample_id:
                    errors.append(f"Row {row_num}: sample_id is empty")
                elif not re.match(r'^[a-zA-Z0-9_\-\.]+$', sample_id):
                    errors.append(f"Row {row_num}: Invalid sample_id format: {sample_id}")

                # Validate chemistry
                if not chemistry:
                    errors.append(f"Row {row_num}: chemistry is empty")

                # Validate fastq_dir
                if not fastq_dir:
                    errors.append(f"Row {row_num}: fastq_dir is empty")
                else:
                    fastq_path = Path(fastq_dir)
                    if not fastq_path.exists():
                        errors.append(f"Row {row_num}: fastq_dir does not exist: {fastq_dir}")
                    elif check_files and fastq_path.is_dir():
                        # Check for FASTQ files
                        r1_files = list(fastq_path.glob('*_R1*.fastq.gz'))
                        r2_files = list(fastq_path.glob('*_R2*.fastq.gz'))

                        if not r1_files or not r2_files:
                            errors.append(
                                f"Row {row_num}: No R1/R2 paired FASTQ files found in {fastq_dir}"
                            )
                        elif len(r1_files) != len(r2_files):
                            errors.append(
                                f"Row {row_num}: Mismatched R1/R2 file counts in {fastq_dir}"
                            )

    except csv.Error as e:
        errors.append(f"CSV parsing error: {e}")
        return False, errors

    except Exception as e:
        errors.append(f"Unexpected error: {e}")
        return False, errors

    return len(errors) == 0, errors


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Validate scRNA-seq samplesheet'
    )
    parser.add_argument('samplesheet', help='Samplesheet TSV file')
    parser.add_argument(
        '--check-files',
        action='store_true',
        help='Verify FASTQ files exist'
    )

    args = parser.parse_args()

    # Validate
    is_valid, errors = check_samplesheet(args.samplesheet, check_files=args.check_files)

    if is_valid:
        print(f"✓ Samplesheet '{args.samplesheet}' is valid")
        sys.exit(0)
    else:
        print(f"✗ Samplesheet '{args.samplesheet}' has errors:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)


if __name__ == '__main__':
    main()
