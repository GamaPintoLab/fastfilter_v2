#!/usr/bin/env python3.9

"""
FastFilter2: High-performance paired-end FASTQ filter for STAR
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon
"""

import argparse
import gzip
import itertools
from pathlib import Path
import threading
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from Bio.SeqIO.QualityIO import FastqPhredIterator

# --------------------------- Defaults --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
CHUNK_SIZE = 10_000  # smaller chunk to reduce freezing
THREADS_DEFAULT = 1

# --------------------------- Sequence Filtering --------------------------- #
def analyze_sequence(record, min_len, min_score, homopolymer_coeff):
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])
    if len(seq) < min_len: return False
    if sum(quals)/len(quals) < min_score if quals else 0: return False
    if max(len(list(g)) for _,g in itertools.groupby(seq)) > homopolymer_coeff: return False
    if "N" in seq or "." in seq: return False
    return True

def process_chunk(records, min_len, min_score, homopolymer_coeff, threads):
    """Filter a chunk using multiple threads."""
    if threads <= 1:
        return [rec.format("fastq") for rec in records if analyze_sequence(rec, min_len, min_score, homopolymer_coeff)]
    filtered = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(analyze_sequence, r, min_len, min_score, homopolymer_coeff) for r in records]
        for rec, f in zip(records, futures):
            if f.result():
                filtered.append(rec.format("fastq"))
    return filtered

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path, dry_run=False):
    if dry_run:
        return open("/dev/null", "w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- Processing --------------------------- #
def process_paired_file(r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, threads, dry_run):
    sample_name = r1_file.stem.replace("_R1","")
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

        pbar = tqdm(unit="reads", desc=sample_name)

        for rec1, rec2 in zip(r1_iter, r2_iter):
            batch.append((rec1, rec2))
            if len(batch) >= CHUNK_SIZE:
                r1_records = [r for r,_ in batch]
                r2_records = [r for _,r in batch]

                filtered_r1 = process_chunk(r1_records, min_len, min_score, homopolymer_coeff, threads)
                filtered_r2 = process_chunk(r2_records, min_len, min_score, homopolymer_coeff, threads)

                out_r1.writelines(filtered_r1)
                out_r2.writelines(filtered_r2)

                # Stats
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

        # Remaining reads
        if batch:
            r1_records = [r for r,_ in batch]
            r2_records = [r for _,r in batch]
            filtered_r1 = process_chunk(r1_records, min_len, min_score, homopolymer_coeff, threads)
            filtered_r2 = process_chunk(r2_records, min_len, min_score, homopolymer_coeff, threads)
            out_r1.writelines(filtered_r1)
            out_r2.writelines(filtered_r2)
            for r1,r2 in batch:
                total_reads += 1
                pass1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
                pass2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
                passed_r1 += pass1
                passed_r2 += pass2
                if pass1 and pass2: passed_pairs += 1
            pbar.update(len(batch))

        pbar.close()

    print(f"{sample_name}: {passed_pairs}/{total_reads} pairs passed ({passed_r1} R1 / {passed_r2} R2)")
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
    with open(summary_path,"w") as f:
        f.write("Sample\tInput_Reads\tPassed_Pairs\tPassed_R1\tPassed_R2\tPercent_Pairs_Passed\n")
        for r in results:
            percent = r["passed_pairs"]/r["input_reads"]*100 if r["input_reads"]>0 else 0
            f.write(f"{r['sample']}\t{r['input_reads']}\t{r['passed_pairs']}\t{r['passed_r1']}\t{r['passed_r2']}\t{percent:.2f}\n")
    print(f"Summary written to {summary_path}")

# --------------------------- Main --------------------------- #
def main():
    import argparse
    import time
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-j","--threads", type=int, default=THREADS_DEFAULT)
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter2"
    output_dir.mkdir(exist_ok=True)

    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}

    paired_files = [(r1_dict[k], r2_dict[k], output_dir, args.minlen, args.min_score, args.homopolymerlen, args.threads, args.dryrun) 
                    for k in sorted(r1_dict.keys()) if k in r2_dict]

    import time
    start = time.time()
    results = []
    for pair in paired_files:
        res = process_paired_file(*pair)
        results.append(res)

    write_summary(results, output_dir)
    elapsed = (time.time() - start)/60
    print(f"Filtering completed in {elapsed:.2f} minutes")

if __name__ == "__main__":
    main()
