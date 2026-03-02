# FastFilter2 Upgrade Guide

This document outlines the changes, improvements, and migration notes for users moving from **FastFilter v1** to **FastFilter2** (v2.0).

---

## **Version Overview**

| Feature | FastFilter v1 | FastFilter2 |
|---------|---------------|-------------|
| Input file support | `.fastq` | `.fastq` & `.fastq.gz` |
| Processing | Single-threaded or multiprocessing | Multi-threaded via `multiprocessing.Pool` |
| Paired-end support | Yes, basic R1/R2 matching | Yes, strict ID matching with synchronized filtering and progress bars |
| Output | Uncompressed FASTQ | gzipped FASTQ (`.fastq.gz`) compatible with STAR |
| Logging | Print statements | Thread-safe logging with timestamps and progress feedback |
| Summary | CSV & TXT reports | Continuous CSV/TSV summary with total reads, passed reads, and percent passed per sample |
| Dry-run | Yes | Yes, fully supported; no output written |
| Progress monitoring | None | Stacked progress bars per sample with `tqdm` |
| Batch writing | No | Yes, faster I/O with configurable batch sizes |
| Homopolymer detection | Yes, basic | Yes, optimized for speed and accuracy |
| Quality filtering | Average Phred score | Average Phred score with configurable thresholds |
| Error handling | Minimal | Robust: R1/R2 mismatches terminate with clear logs |

---

## **Key Improvements in FastFilter2**

1. **True Multi-threading for Large Datasets**  
   - Each paired-end sample is processed in a separate worker process.  
   - Reduces runtime significantly on multi-core systems.  
   - Visual stacked progress bars show live progress for all samples.

2. **Gzipped Input and Output**  
   - Supports `.fastq.gz` input files directly.  
   - Outputs are automatically compressed with `pigz` for storage efficiency and STAR compatibility.

3. **Improved Filtering Logic**  
   - Minimum sequence length, homopolymers, ambiguous bases, and Phred quality scores enforced per read.  
   - Reads fail if any criterion is not met, maintaining data integrity for downstream analysis.

4. **Robust Paired-End Validation**  
   - Ensures R1 and R2 reads remain synchronized.  
   - Detects and stops processing if R1/R2 files mismatch in read counts or IDs.

5. **Continuous Summary Reports**  
   - Generates `fastfilter_summary.csv` or TSV with total reads, reads passing each filter, and pass rates.  
   - Includes per-sample statistics for better QC tracking.

6. **Dry-Run Mode**  
   - Allows testing the pipeline without writing files.  
   - Useful for validating parameter choices or debugging workflows.

7. **Configurable Parameters via CLI**  
   - Minimum sequence length (`-l`)  
   - Minimum Phred score (`-s`)  
   - Maximum homopolymer length (`-p`)  
   - Number of CPU threads (`-j`)  
   - Input/output directories (`-i`, `-o`)

8. **Better Memory & I/O Management**  
   - Batch writing reduces file I/O overhead.  
   - Handles large datasets efficiently without overloading RAM.

---

## **Migration Notes**

- **Output File Names:**  
  - Old: `<sample>_FILTERED.fastq`  
  - New: `<sample>_R1_FILTERED.fastq.gz` and `<sample>_R2_FILTERED.fastq.gz`

- **Summary Format:**  
  - Old: CSV/TXT with per-read tables  
  - New: CSV/TSV per sample with read counts, filtering reasons, and pass rates

- **Threading & Progress Bars:**  
  - Multi-threaded processes replace older multiprocessing code.  
  - Progress bars now stacked for all samples.

- **Dry-Run Behavior:**  
  - Outputs go to memory/dummy buffers instead of being skipped silently.

- **Error Handling:**  
  - Mismatched R1/R2 files terminate processing immediately with clear error messages.  
  - Improves reproducibility and avoids silent failures.

---

## **Example Migration Command**

**Old command (v1):**

python fastfilter.py -i /data/project/cutadapt -o /data/project/fastfilter -l 30  

**New command (v2):**

python fastfilter_pe.py -i /data/project/cutadapt -o /data/project/fastfilter -l 30 -s 30 -p 25 -j 4  

- Adds explicit CPU threads (`-j 4`)  
- Adds Phred score filter (`-s 30`)  
- Adds homopolymer maximum length (`-p 25`)  

---

## **Recommendations**

1. Test **FastFilter2** on a small dataset first using `--dryrun`.  
2. Validate that STAR or downstream aligners accept gzipped outputs.  
3. Review `fastfilter_summary.csv` for thresholds, pass rates, and QC statistics.  
4. Adjust `-l`, `-s`, `-p` parameters according to your experimental design.  

---

## **References**

- FastFilter original implementation: [fastfilter GitHub](https://github.com/GamaPintoLab/fastfilter/blob/main/README.md)  
- FastFilter2 GitHub: [fastfilter2](https://github.com/GamaPintoLab/fastfilter2/blob/main/README.md)  
- STAR aligner: [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
