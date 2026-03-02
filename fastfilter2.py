#!/usr/bin/env python3

"""
FastFilter2: High-performance paired-end FASTQ filter for STAR
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon

Features:
- Filters sequences by length, quality, homopolymer runs, and N characters
- Supports .fastq and .fastq.gz input
- Produces .fastq.gz output compatible with STAR
- Multiprocessing with batch processing for maximum CPU utilization
- Dry-run option
- TSV summary of all samples
- Stacked progress bar per sample using tqdm
"""

import argparse
import itertools
import gzip
import csv
import time
from pathlib import Path
from multiprocessing import cpu_count
from tqdm.contrib.concurrent import process_map
from Bio.SeqIO.QualityIO import FastqPhredIterator

# --------------------------- Default Parameters --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
CHUNK_SIZE = 100_000  # Reads per batch for faster processing

# --------------------------- Sequence Analysis --------------------------- #
def analyze_sequence(record, min_len, min_score, homopolymer_coeff):
    """Return True if a sequence passes all filters."""
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])

    if len(seq) < min_len:
        return False
    mean_q = sum(quals) / len(quals) if quals else 0
    if mean_q < min_score:
        return False
    if max((len(list(g)) for _, g in itertools.groupby(seq)), default=0) > homopolymer_coeff:
        return False
    if "N" in seq or "." in seq:
        return False
    return True

def process_chunk(records, min_len, min_score, homopolymer_coeff):
    """Process a chunk of records. Returns filtered records as FASTQ strings."""
    return [rec.format("fastq") for rec in records if analyze_sequence(rec, min_len, min_score, homopolymer_coeff)]

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    """Open FASTQ or gzipped FASTQ file for reading."""
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path, dry_run=False):
    """Open gzipped FASTQ output or /dev/null for dry-run."""
    if dry_run:
        return open("/dev/null", "w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- File Parsing --------------------------- #
def parse_file_pair(args):
    r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, dry_run = args
    sample_name = r1_file.stem.replace("_R1", "")
    out_r1 = open_output(output_dir / f"{r1_file.stem}_FILTERED", dry_run)
    out_r2 = open_output(output_dir / f"{r2_file.stem}_FILTERED", dry_run)

    total_reads = 0
    passed_pairs = 0
    passed_r1 = 0
    passed_r2 = 0

    with open_input(r1_file) as f1, open_input(r2_file) as f2, out_r1, out_r2:
        r1_iter = FastqPhredIterator(f1)
        r2_iter = FastqPhredIterator(f2)
        batch = []

        for rec1, rec2 in itertools.zip_longest(r1_iter, r2_iter):
            if rec1 is None or rec2 is None:
                raise RuntimeError(f"Unequal reads in {r1_file.name} and {r2_file.name}")
            if rec1.id != rec2.id:
                raise RuntimeError(f"Read ID mismatch:\n{rec1.id}\n{rec2.id}")

            batch.append((rec1, rec2))
            if len(batch) >= CHUNK_SIZE:
                r1_records = [r for r,_ in batch]
                r2_records = [r for _,r in batch]
                filtered_r1 = process_chunk(r1_records, min_len, min_score, homopolymer_coeff)
                filtered_r2 = process_chunk(r2_records, min_len, min_score, homopolymer_coeff)
                out_r1.writelines(filtered_r1)
                out_r2.writelines(filtered_r2)

                # Stats
                for i, (r1,r2) in enumerate(batch):
                    total_reads += 1
                    pass1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
                    pass2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
                    passed_r1 += pass1
                    passed_r2 += pass2
                    if pass1 and pass2:
                        passed_pairs += 1
                batch.clear()

        # Process remaining reads
        if batch:
            r1_records = [r for r,_ in batch]
            r2_records = [r for _,r in batch]
            filtered_r1 = process_chunk(r1_records, min_len, min_score, homopolymer_coeff)
            filtered_r2 = process_chunk(r2_records, min_len, min_score, homopolymer_coeff)
            out_r1.writelines(filtered_r1)
            out_r2.writelines(filtered_r2)

            for i, (r1,r2) in enumerate(batch):
                total_reads += 1
                pass1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
                pass2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
                passed_r1 += pass1
                passed_r2 += pass2
                if pass1 and pass2:
                    passed_pairs += 1

    return {
        "sample": sample_name,
        "input_reads": total_reads,
        "passed_pairs": passed_pairs,
        "passed_r1": passed_r1,
        "passed_r2": passed_r2
    }

# --------------------------- Summary --------------------------- #
def write_summary(results, output_dir):
    summary_path = output_dir / "filtering_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample","Input_Reads","Passed_Pairs","Passed_R1","Passed_R2","Percent_Pairs_Passed"])
        for r in results:
            percent = (r["passed_pairs"]/r["input_reads"]*100) if r["input_reads"]>0 else 0
            writer.writerow([r["sample"], r["input_reads"], r["passed_pairs"], r["passed_r1"], r["passed_r2"], f"{percent:.2f}"])
    print(f"Summary written to: {summary_path}")

# --------------------------- Main --------------------------- #
def main():
    parser = argparse.ArgumentParser(description="FastFilter3: Multiprocessing paired-end FASTQ filter")
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-j","--cpus", type=int, default=cpu_count())
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect paired-end files
    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}

    missing_r2 = set(r1_dict) - set(r2_dict)
    missing_r1 = set(r2_dict) - set(r1_dict)
    if missing_r2: raise RuntimeError(f"Missing R2 files for: {missing_r2}")
    if missing_r1: raise RuntimeError(f"Missing R1 files for: {missing_r1}")

    paired_files = [(r1_dict[k], r2_dict[k], output_dir, args.minlen, args.min_score, args.homopolymerlen, args.dryrun)
                    for k in sorted(r1_dict.keys())]

    start_time = time.time()
    results = process_map(parse_file_pair, paired_files, max_workers=args.cpus, chunksize=1)
    write_summary(results, output_dir)
    elapsed = (time.time() - start_time)/60
    print(f"Filtering completed in {elapsed:.2f} minutes")

if __name__ == "__main__":
    main()
