# FastFilter2

**FastFilter2** is a high-performance, threaded FASTQ filter designed for paired-end sequencing data, optimized for downstream analysis with **STAR**. It efficiently filters sequences based on length, quality, homopolymer runs, and ambiguous nucleotides while supporting both `.fastq` and `.fastq.gz` formats. The script produces compressed outputs compatible with STAR and generates a continuous TSV summary of filtering results.

---

## **Table of Contents**

- [Features](#features)  
- [Installation](#installation)  
- [Usage](#usage)  
- [Command-Line Arguments](#command-line-arguments)  
- [Workflow](#workflow)  
- [Output Files](#output-files)  
- [Example](#example)  
- [License](#license)  
- [Acknowledgements](#acknowledgements)  

---

## **Features**

- Multi-threaded filtering using Python’s `ThreadPoolExecutor`.  
- Filters sequences by:  
  - Minimum length  
  - Average Phred quality score  
  - Homopolymer runs  
  - Presence of ambiguous nucleotides (`N`) or dots (`.`)  
- Supports both `.fastq` and `.fastq.gz` inputs.  
- Produces gzipped `.fastq.gz` outputs compatible with STAR.  
- Stacked progress bars per sample using `tqdm`.  
- Dry-run mode for testing without generating output files.  
- Continuous TSV summary updates during processing.  
- Thread-safe logging for professional progress monitoring.  

---

## **Installation**

1. Clone the repository:

```sh
git clone https://github.com/yourusername/FastFilter2.git  
cd FastFilter2  
```

2. Install dependencies (recommended in a virtual environment):

```sh
pip install biopython tqdm  
```

> Python 3.8 or higher is required.  

---

## **Usage**

Run FastFilter2 from the command line:

python FastFilter2.py -i /path/to/sequences -o /path/to/output -j 4  

- `-i` specifies the directory containing paired-end FASTQ files.  
- `-o` specifies the output directory. If not provided, a `fastfilter` folder will be created next to the input directory.  
- `-j` specifies the number of CPU threads to use.  

Dry-run example (does not write output files):

```sh
python FastFilter2.py -i /path/to/sequences -d  
```

---

## **Command-Line Arguments**

- `-i, --sequences-dir` **(required)**: Path to the input FASTQ/FASTQ.gz folder.  
- `-o, --output-dir`: Path for the output directory (default: `<input_dir>/fastfilter`).  
- `-l, --minlen`: Minimum sequence length to retain (default: 25).  
- `-s, --min-score`: Minimum average Phred quality score (default: 30).  
- `-p, --homopolymerlen`: Maximum allowed homopolymer run length (default: 25).  
- `-j, --cpus`: Number of threads to use for filtering (default: 1).  
- `-d, --dryrun`: Execute a dry run without generating output files.  

---

## **Workflow**

1. **File Discovery:** Detects paired-end FASTQ files (`*_R1*.fastq*` and `*_R2*.fastq*`).  
2. **Read Filtering:** Each R1/R2 pair is filtered in a separate thread:  
   - Compute sequence length and average Phred score.  
   - Detect maximum homopolymer length.  
   - Check for ambiguous nucleotides (`N`) or dots (`.`).  
   - Write passing reads to gzipped FASTQ files.  
3. **Progress Monitoring:** Each thread has a dedicated `tqdm` progress bar.  
4. **Summary Generation:** Continuously updates a TSV summary with per-sample statistics:  
   - Input reads  
   - Passed read pairs  
   - Reads passing individually in R1 and R2  
   - Percentage of read pairs passed  
5. **Completion Logging:** Reports total execution time.  

---

## **Output Files**

All output is stored in the specified output directory.

- **Filtered FASTQ Files:**  
  - `<sample>_R1_FILTERED.fastq.gz`  
  - `<sample>_R2_FILTERED.fastq.gz`  

- **Summary File:**  
  - `filtering_summary.tsv` – Tab-separated summary of all samples, including:  
    - Sample name  
    - Input reads  
    - Passed pairs  
    - Passed R1 and R2 reads  
    - Percent of pairs passing filters  

---

## **Example**

Filter paired-end FASTQ files with 4 threads:

```sh
python FastFilter2.py -i /data/project/fastq -o /data/project/fastfilter -j 4  
```

Dry-run test:

```sh
python FastFilter2.py -i /data/project/fastq -d  
```

Example snippet from `filtering_summary.tsv`:

Sample    Input_Reads    Passed_Pairs    Passed_R1    Passed_R2    Percent_Pairs_Passed  
sample1    1000000        950000          970000       960000       95.00  

---

## **License**

This project is released under the MIT License. See LICENSE file for details.  

---

## **Acknowledgements**

- **Author:** Lucas Monteiro  
- **PI:** Margarida Gama-Carvalho  
- **Lab:** RNA Systems Biology Lab, BioISI, University of Lisbon  
- Inspired by previous FastFilter implementations for RNA-Seq preprocessing.  

---

---

## License

FastFilter is released under the MIT License.

---

**Author:** Lucas Monteiro  
**PI:** Margarida Gama-Carvalho  
**Lab:** RNA Systems Biology Lab, BioISI, University of Lisbon
