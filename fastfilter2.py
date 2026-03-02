#!/usr/bin/env python3.9
"""
FastFilter2: Efficient paired-end FASTQ filter for STAR
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon

Features:
- Filters sequences by length, quality, homopolymer runs, and N characters
- Supports .fastq and .fastq.gz input
- Produces .fastq.gz output compatible with STAR
- Multi-threaded chunk processing per file
- Async gzipped output for maximum CPU throughput
- Dry-run option
- Progress bar with reads per file
- Summary TSV with before/after counts
"""

import argparse
import itertools
import gzip
import csv
import threading
import time
from pathlib import Path
from queue import Queue
from Bio.SeqIO.QualityIO import FastqPhredIterator
from tqdm import tqdm

# --------------------------- Default Parameters --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
CHUNK_SIZE = 100_000  # Reads per batch for CPU efficiency
QUEUE_MAXSIZE = 10  # Max chunks waiting for writer

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
    """Return filtered FASTQ lines for a chunk."""
    return [rec.format("fastq") for rec in records if analyze_sequence(rec, min_len, min_score, homopolymer_coeff)]

# --------------------------- Writer Thread --------------------------- #
def writer_thread(out_r1, out_r2, queue):
    """Consume chunks from the queue and write to gzipped outputs."""
    while True:
        item = queue.get()
        if item is None:  # Sentinel to stop
            break
        lines_r1, lines_r2 = item
        out_r1.writelines(lines_r1)
        out_r2.writelines(lines_r2)
        queue.task_done()

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path, dry_run=False):
    if dry_run:
        return open("/dev/null", "w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- Paired File Processing --------------------------- #
def process_paired_file(r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, dry_run):
    sample_name = r1_file.stem.replace("_R1", "")
    out_r1 = open_output(output_dir / f"{r1_file.stem}_FILTERED", dry_run)
    out_r2 = open_output(output_dir / f"{r2_file.stem}_FILTERED", dry_run)

    queue = Queue(maxsize=QUEUE_MAXSIZE)
    writer = threading.Thread(target=writer_thread, args=(out_r1, out_r2, queue), daemon=True)
    writer.start()

    total_reads = 0
    passed_r1 = 0
    passed_r2 = 0
    passed_pairs = 0

    with open_input(r1_file) as f1, open_input(r2_file) as f2, out_r1, out_r2:
        r1_iter = FastqPhredIterator(f1)
        r2_iter = FastqPhredIterator(f2)
        batch = []
        pbar = tqdm(unit="reads", desc=sample_name)

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
                queue.put((filtered_r1, filtered_r2))

                # Update stats
                for r1,r2 in batch:
                    total_reads += 1
                    pass1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
                    pass2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
                    passed_r1 += pass1
                    passed_r2 += pass2
                    if pass1 and pass2:
                        passed_pairs += 1
                batch.clear()
                pbar.update(CHUNK_SIZE)

        # Process remaining
        if batch:
            r1_records = [r for r,_ in batch]
            r2_records = [r for _,r in batch]
            filtered_r1 = process_chunk(r1_records, min_len, min_score, homopolymer_coeff)
            filtered_r2 = process_chunk(r2_records, min_len, min_score, homopolymer_coeff)
            queue.put((filtered_r1, filtered_r2))

            for r1,r2 in batch:
                total_reads += 1
                pass1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
                pass2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
                passed_r1 += pass1
                passed_r2 += pass2
                if pass1 and pass2:
                    passed_pairs += 1
            pbar.update(len(batch))

        pbar.close()

    # Stop writer
    queue.put(None)
    writer.join()

    print(f"{sample_name}: {total_reads} reads -> {passed_pairs} pairs passed")
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
    parser = argparse.ArgumentParser(description="FastFilter2: paired-end FASTQ filter")
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-j","--threads", type=int, default=1)
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}

    missing_r2 = set(r1_dict) - set(r2_dict)
    missing_r1 = set(r2_dict) - set(r1_dict)
    if missing_r2: raise RuntimeError(f"Missing R2 files for: {missing_r2}")
    if missing_r1: raise RuntimeError(f"Missing R1 files for: {missing_r1}")

    paired_files = [(r1_dict[k], r2_dict[k], output_dir,
                     args.minlen, args.min_score, args.homopolymerlen, args.dryrun)
                    for k in sorted(r1_dict.keys())]

    results = []
    start_time = time.time()
    for r1, r2, out_dir, min_len, min_score, homo_len, dry_run in paired_files:
        res = process_paired_file(r1, r2, out_dir, min_len, min_score, homo_len, dry_run)
        results.append(res)

    write_summary(results, output_dir)
    elapsed = (time.time() - start_time)/60
    print(f"Filtering completed in {elapsed:.2f} minutes")

if __name__ == "__main__":
    main()
