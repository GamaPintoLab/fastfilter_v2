# FastFilter2

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

**FastFilter2** is a high-performance, production-ready Python tool for filtering paired-end FASTQ files. Designed for bioinformatics pipelines, it provides flexible, reliable, and fast filtering of sequencing data with built-in support for multi-threading, compression, and detailed summaries.

This tool is ideal for pre-processing RNA-seq, DNA-seq, or other high-throughput sequencing datasets prior to alignment, assembly, or downstream analysis.

---

## Key Features

- **Biological quality filters**:
  - Minimum read length
  - Maximum allowed ambiguous bases (Ns)
  - Homopolymer detection
  - Minimum average Phred quality score
- **Paired-end safe**: Ensures that reads are filtered in pairs, maintaining synchronization between R1 and R2.
- **High-performance I/O**: Writes uncompressed FASTQ first for speed, then compresses output automatically with `pigz` using multiple threads.
- **Batch processing**: Efficient batch writing to reduce I/O overhead.
- **Progress tracking**: Real-time progress bars via `tqdm` for monitoring large datasets.
- **Summary output**: Generates CSV reports with total reads, passing reads, and pass rates.

---

## Installation

Clone the repository and install dependencies:

```sh
git clone https://github.com/GamaPintoLab/fastfilter2.git
cd fastfilter2
pip install -r requirements.txt
```

**Dependencies**:

- Python 3.9 or higher
- Biopython
- tqdm
- pigz (for multi-threaded compression)

---

## Usage

Run the tool from the command line:

```sh
fastfilter2 -i /path/to/input_fastq_dir -o /path/to/output_dir -j 4
```

### Command-line Options

- `-i, --seq-dir` : Input directory containing paired-end FASTQ files (required)
- `-o, --output-dir` : Directory to write filtered outputs (defaults to `<input_dir>/fastfilter`)
- `-j, --cpus` : Number of threads for parallel processing and compression (default: 1)
- `-l, --minlen` : Minimum sequence length (default: 25)
- `-p, --homopolymerlen` : Maximum allowed homopolymer length (default: 25)
- `-s, --min-score` : Minimum average Phred quality score (default: 30)
- `--dryrun` : Run without writing outputs (for testing)

### Example

```sh
fastfilter2 -i samples/fastq -o results/filtered -j 8 -l 50 -s 20 --dryrun
```

This example processes paired-end FASTQ files in `samples/fastq` using 8 CPU threads, filters reads shorter than 50 bases or with average quality below 20, and performs a dry run without writing files.

---

## How It Works

1. **Input parsing**: Reads paired-end FASTQ files and validates file pairs.
2. **Filtering**: Applies multiple biological filters to each read:
   - Removes reads with ambiguous bases (N or .)
   - Filters out reads with homopolymers above a given threshold
   - Filters based on minimum length and average Phred score
3. **Batch writing**: Writes passing reads in batches to reduce I/O overhead.
4. **Compression**: Automatically compresses output FASTQ files with `pigz` for speed and storage efficiency.
5. **Reporting**: Produces a summary CSV with per-file statistics including total reads, passing reads, and pass rates.

---

## Output

Filtered paired-end files are named:

`<sample_name>_R1_FILTERED.fastq.gz`  
`<sample_name>_R2_FILTERED.fastq.gz`  

Summary CSV:

`fastfilter_summary.csv` containing:

- file: sample name
- total_reads: number of reads in input
- good_reads: reads passing filters
- pass_rate_pct: percentage of reads passing filters

---

## Performance

- Multi-threaded filtering and compression using `multiprocessing` and `pigz`.
- Efficient memory usage via batch processing of reads.
- Suitable for very large FASTQ datasets (tens to hundreds of millions of reads).

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Acknowledgements

- Built on top of the original FastFilter concept  
- Biopython community for sequence handling tools  
- `tqdm` for progress visualization  
- `pigz` for high-speed parallel compression

---

**Author:** Lucas Monteiro  
**PI:** Margarida Gama-Carvalho  
**Lab:** RNA Systems Biology Lab, BioISI, University of Lisbon
