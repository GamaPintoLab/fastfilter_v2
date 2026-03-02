# FastFilter

**FastFilter** is a high-performance, threaded paired-end FASTQ filter designed for preprocessing RNA-seq data before alignment with STAR. It supports standard and gzipped FASTQ files, applies customizable quality and sequence filters, and produces STAR-compatible gzipped output with detailed logging and summary statistics.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Command-Line Options](#command-line-options)
- [Filtering Criteria](#filtering-criteria)
- [Output](#output)
- [Multi-threading & Progress](#multi-threading--progress)
- [Examples](#examples)
- [Logging & Summary](#logging--summary)
- [License](#license)

---

## Features

- Filters paired-end FASTQ files by:
  - Minimum sequence length
  - Mean base quality
  - Homopolymer runs
  - Ambiguous bases (`N` or `.`)
- Supports `.fastq` and `.fastq.gz` inputs
- Produces `.fastq.gz` outputs compatible with STAR
- Multi-threaded processing with stacked progress bars per sample
- Dry-run mode for testing
- Continuous TSV summary updates
- Logging-based progress tracking instead of print statements

---

## Installation

FastFilter is implemented in Python 3 and requires the following dependencies:

- Python >= 3.7
- [Biopython](https://biopython.org/)
- [tqdm](https://tqdm.github.io/)

You can install dependencies via pip:

pip install biopython tqdm

Clone the repository:

git clone https://github.com/yourusername/FastFilter.git
cd FastFilter

Make the script executable:

chmod +x fastfilter.py

---

## Usage

Run FastFilter from the command line:

./fastfilter.py -i /path/to/fastq_files -o /path/to/output_dir -j 4

**Arguments:**

- `-i / --sequences-dir` (required): Directory containing paired-end FASTQ files
- `-o / --output-dir` (optional): Directory to write filtered FASTQ files (default: `../fastfilter`)
- `-j / --cpus` (optional): Number of threads for parallel processing (default: 1)
- `-d / --dryrun` (optional): Perform a dry run without writing output

---

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| -l, --minlen | Minimum sequence length to retain | 25 |
| -s, --min-score | Minimum mean Phred quality score | 30 |
| -p, --homopolymerlen | Maximum allowed homopolymer run | 25 |
| -i, --sequences-dir | Directory containing input FASTQ files | Required |
| -o, --output-dir | Directory for filtered FASTQ files | `../fastfilter` |
| -j, --cpus | Number of threads to use | 1 |
| -d, --dryrun | Perform filtering without writing files | False |

---

## Filtering Criteria

Each read in a paired-end FASTQ file is evaluated against the following criteria:

1. **Minimum Length:** Reads shorter than the specified length are discarded.
2. **Mean Quality Score:** Reads with mean Phred score below the threshold are discarded.
3. **Homopolymer Runs:** Reads containing homopolymer runs longer than the threshold are discarded.
4. **Ambiguous Bases:** Reads containing `N` or `.` are discarded.

A pair is retained only if both reads pass all filters.

---

## Output

Filtered reads are written as gzipped FASTQ files with `_FILTERED` suffix:

- `<sample>_R1_FILTERED.fastq.gz`
- `<sample>_R2_FILTERED.fastq.gz`

Additionally, a TSV summary is generated:

`filtering_summary.tsv`

Columns include:

- Sample
- Input Reads
- Passed Pairs
- Passed R1
- Passed R2
- Percent Pairs Passed

---

## Multi-threading & Progress

FastFilter uses Python's `concurrent.futures.ThreadPoolExecutor` to process multiple samples in parallel. Each sample displays a dedicated progress bar using `tqdm` for real-time monitoring.

---

## Examples

Filter a directory with 4 threads:

./fastfilter.py -i /data/fastq_samples -o /data/fastq_filtered -j 4

Perform a dry run without writing output:

./fastfilter.py -i /data/fastq_samples -d

Filter using custom thresholds:

./fastfilter.py -i /data/fastq_samples -l 50 -s 35 -p 20 -j 8

---

## Logging & Summary

All progress, warnings, and errors are logged using Python's `logging` module with timestamps. At the end of the run, a detailed TSV summary is written, including read counts before and after filtering.

Example log entry:

2026-03-02 14:22:05 | INFO | Processing sample: Sample_01  
2026-03-02 14:25:13 | INFO | Finished sample: Sample_01 | 950000/1000000 pairs passed  

---

## License

FastFilter is released under the MIT License.

---

**Author:** Lucas Monteiro  
**PI:** Margarida Gama-Carvalho  
**Lab:** RNA Systems Biology Lab, BioISI, University of Lisbon
