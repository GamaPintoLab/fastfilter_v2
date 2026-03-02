#!/usr/bin/env python3.9
import argparse
import gzip
import multiprocessing
from pathlib import Path
import time
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqPhredIterator
import csv

# -------------------- Defaults -------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30

# -------------------- Globals -------------------- #
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
seq_dir = None
output_dir = None
dryrun = False


# -------------------- Utilities -------------------- #
def find_homopolymers(seq: str) -> dict:
    """Return dict indicating if sequence contains homopolymers of length homopolymer_coeff"""
    return {
        "A": int("A" * homopolymer_coeff in seq),
        "T": int("T" * homopolymer_coeff in seq),
        "G": int("G" * homopolymer_coeff in seq),
        "C": int("C" * homopolymer_coeff in seq),
    }


def analyze_sequence(record):
    """Analyze sequence and return pass/fail and exclusion reasons"""
    qual_values = record.letter_annotations.get("phred_quality", [])
    qual_score = sum(qual_values)/len(qual_values) if qual_values else 0
    seq = str(record.seq)
    homopolymers = find_homopolymers(seq)
    homopolymer_exist = any(homopolymers.values())
    n_count = seq.count("N")
    dot_count = seq.count(".")

    meets_criteria = (
        len(seq) >= min_seq_len and
        n_count == 0 and
        dot_count == 0 and
        not homopolymer_exist and
        qual_score >= min_score
    )

    exclusion_reasons = {
        "too_short": len(seq) < min_seq_len,
        "found_n": n_count > 0,
        "found_dot": dot_count > 0,
        "found_homopolymer": homopolymer_exist,
        "low_score": qual_score < min_score
    }

    return meets_criteria, exclusion_reasons


# -------------------- Paired-end parser -------------------- #
def parse_paired_file(r1_path: Path, r2_path: Path):
    """Parse a pair of R1/R2 FASTQ files, write filtered .fastq.gz, return stats"""
    base_name = r1_path.stem.replace("_R1", "")
    out_r1 = output_dir / f"{r1_path.stem}_FILTERED.fastq.gz"
    out_r2 = output_dir / f"{r2_path.stem}_FILTERED.fastq.gz"

    total_reads = 0
    good_reads = 0
    exclusion_counts_r1 = {k:0 for k in ["too_short","found_n","found_dot","found_homopolymer","low_score"]}
    exclusion_counts_r2 = exclusion_counts_r1.copy()

    with gzip.open(out_r1, "wt") as f_out1, gzip.open(out_r2, "wt") as f_out2:
        with open(r1_path) as f1, open(r2_path) as f2:
            r1_iter = FastqPhredIterator(f1)
            r2_iter = FastqPhredIterator(f2)
            for rec1, rec2 in zip(r1_iter, r2_iter):
                total_reads += 1
                passes1, reasons1 = analyze_sequence(rec1)
                passes2, reasons2 = analyze_sequence(rec2)

                if passes1 and passes2:
                    good_reads += 1
                    SeqIO.write(rec1, f_out1, "fastq")
                    SeqIO.write(rec2, f_out2, "fastq")

                for k, v in reasons1.items():
                    exclusion_counts_r1[k] += v
                for k, v in reasons2.items():
                    exclusion_counts_r2[k] += v

    return {
        "file": base_name,
        "total_reads": total_reads,
        "good_reads": good_reads,
        **{f"R1_{k}": v for k,v in exclusion_counts_r1.items()},
        **{f"R2_{k}": v for k,v in exclusion_counts_r2.items()}
    }


# -------------------- Argument parser -------------------- #
def parse_arguments():
    global min_seq_len, homopolymer_coeff, min_score, seq_dir, output_dir, dryrun
    parser = argparse.ArgumentParser(description="Paired-end FASTQ filter with .fastq.gz output")
    parser.add_argument("-l", "--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p", "--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s", "--minscore", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i", "--inputdir", type=str, required=True)
    parser.add_argument("-o", "--outputdir", type=str)
    parser.add_argument("--dryrun", action="store_true")
    args = parser.parse_args()

    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.minscore
    dryrun = args.dryrun
    seq_dir = Path(args.inputdir)
    output_dir = Path(args.outputdir) if args.outputdir else seq_dir / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)


# -------------------- Main -------------------- #
def main():
    start_time = time.time()
    parse_arguments()

    fastq_r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    fastq_r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))
    paired_files = list(zip(fastq_r1_files, fastq_r2_files))

    if not paired_files:
        raise ValueError("No paired-end FASTQ files found in input directory.")

    # Multiprocessing
    pool = multiprocessing.Pool()
    results = pool.starmap(parse_paired_file, paired_files)
    pool.close()
    pool.join()

    # Write summary CSV
    summary_file = output_dir / "fastfilter_summary.csv"
    with open(summary_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    elapsed = (time.time() - start_time)/60
    print(f"Filtering finished in {elapsed:.2f} min")
    print(f"Summary saved to {summary_file}")


if __name__ == "__main__":
    main()
