#!/usr/bin/env python3

"""
FastFilter: Threaded Paired-End FASTQ Filter for STAR Alignment

Author: Gil Poiares-Oliveira
PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon

--------------------------------------------------------------------------------
DESCRIPTION
--------------------------------------------------------------------------------
FastFilter is a high-performance, multi-threaded Python tool for filtering
paired-end FASTQ sequences before RNA-seq analysis with STAR. It checks
each read pair against biological and technical quality thresholds to ensure
only high-quality reads are passed downstream.

--------------------------------------------------------------------------------
FEATURES
--------------------------------------------------------------------------------
- Filters paired-end FASTQ files (.fastq or .fastq.gz) by:
    1. Minimum sequence length
    2. Average Phred quality score
    3. Maximum homopolymer run length
    4. Ambiguous nucleotides (N or .)
- Produces compressed output files (.fastq.gz) compatible with STAR.
- Multi-threaded with per-file stacked progress bars using tqdm.
- Dry-run mode for testing parameters without writing outputs.
- Continuous TSV summary updates per sample.
- Uses logging instead of print statements for better tracking.

--------------------------------------------------------------------------------
INPUT
--------------------------------------------------------------------------------
- Directory containing paired FASTQ files.
- Files should follow naming convention:
    - Read 1: *_R1*.fastq or *.fastq.gz
    - Read 2: *_R2*.fastq or *.fastq.gz
- Files must be synchronized (same number of reads and matching IDs).

--------------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------------
- Filtered FASTQ files (compressed):
    - <original_name>_FILTERED.fastq.gz
- Summary file (TSV):
    - filtering_summary.tsv
    - Columns:
        - Sample: sample name
        - Input_Reads: total reads in original files
        - Passed_Pairs: read pairs passing all filters
        - Passed_R1: R1 reads passing filters
        - Passed_R2: R2 reads passing filters
        - Percent_Pairs_Passed: proportion of paired reads passed

--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
Basic filtering:

    python fastfilter.py -i /path/to/sequences -o /path/to/output

Multi-threaded execution (4 samples in parallel):

    python fastfilter.py -i /path/to/sequences -o /path/to/output -j 4

Custom filtering thresholds:

    python fastfilter.py -i /path/to/sequences -o /path/to/output \
                         -l 30 -s 35 -p 20

Dry-run mode (no files written):

    python fastfilter.py -i /path/to/sequences -d

--------------------------------------------------------------------------------
PARAMETERS
--------------------------------------------------------------------------------
- -i / --sequences-dir  : Input directory with paired FASTQ files (required)
- -o / --output-dir     : Output directory for filtered files (default: ./fastfilter)
- -l / --minlen         : Minimum sequence length (default: 25)
- -s / --min-score      : Minimum average Phred quality score (default: 30)
- -p / --homopolymerlen : Maximum homopolymer length (default: 25)
- -j / --cpus           : Number of threads to run in parallel (default: 1)
- -d / --dryrun         : Dry-run mode, does not write output files

--------------------------------------------------------------------------------
EXAMPLES
--------------------------------------------------------------------------------
1. Filter sequences with default thresholds:

    python fastfilter.py -i ./raw_data -o ./filtered_data

2. Use custom thresholds and 8 threads:

    python fastfilter.py -i ./raw_data -o ./filtered_data -l 30 -s 35 -p 20 -j 8

3. Test filtering without writing output:

    python fastfilter.py -i ./raw_data -d

--------------------------------------------------------------------------------
NOTES
--------------------------------------------------------------------------------
- Ensure R1/R2 files are synchronized before running.
- Logs are written to stdout and show sample start/end and filtering stats.
- The tool streams input files to minimize memory usage.
- Recommended for pre-processing before STAR alignment.
"""
