#!/usr/bin/env python3.9

#!/usr/bin/env python3

"""
fastfilter.py - version 2.0
Author: Lucas da Costa Monteiro <ldmonteiro@fc.ul.pt>
PI: Margarida Gama-Carvalho <mhcarvalho@ciencias.ulisboa.pt>

RNA Systems Biology Lab
BioISI - Biosystems and Integrative Sciences Institute
Department of Chemistry and Biochemistry
Faculty of Sciences, University of Lisbon

(C) 2022-2026
"""

import argparse
import csv
from collections import Counter
from datetime import datetime
import multiprocessing
import re
import sys
import time
from pathlib import Path

# isal (Intel ISA-L) provides SIMD-accelerated gzip decompression/compression.
# 2-4x faster than stdlib gzip on large files. Falls back transparently if not installed.
# Install with: pip install isal
try:
    from isal import igzip as gzip
except ImportError:
    import gzip

from tqdm import tqdm

# Defaults
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
LOG_INTERVAL = 5_000_000   # checkpoint print every N reads when not using tqdm

# Compiled once at import time; used in the inner sync loop millions of times
_ID_SUFFIX_RE = re.compile(r'/[12]$')

# Module-level globals (overwritten by parse_arguments at runtime;
# also re-set in each worker via _init_worker after multiprocessing spawn)
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
max_n = 0
single_end = False
compresslevel = 6
_HOMOPOLYMERS: tuple = ()   # precomputed once per process in _init_worker / main()


def parse_arguments() -> None:
    """Parse runtime arguments and populate global variables."""
    parser = argparse.ArgumentParser(description="Filter FASTQ files.")
    parser.add_argument(
        "-l", "--minlen", metavar="50", type=int, default=MIN_LENGTH_DEFAULT,
        help="sequence length threshold",
    )
    parser.add_argument(
        "-p", "--homopolymerlen", metavar="3", type=int, default=HOMOPOLYMER_COEFF_DEFAULT,
        help="homopolymer length",
    )
    parser.add_argument(
        "-s", "--min-score", type=int, metavar="30", default=MIN_SCORE_DEFAULT,
        help="Minimum mean Phred quality score per read. Default: 30.",
    )
    parser.add_argument(
        "-r1", "--r1-files", nargs="+", type=Path, metavar="sample_R1.fastq",
        help="R1 (forward) FASTQ files for paired-end analysis.",
    )
    parser.add_argument(
        "-r2", "--r2-files", nargs="+", type=Path, metavar="sample_R2.fastq",
        help="R2 (reverse) FASTQ files for paired-end analysis.",
    )
    parser.add_argument(
        "-r", "--reads", nargs="+", type=Path, metavar="sample.fastq",
        help="FASTQ files for single-end analysis.",
    )
    parser.add_argument(
        "-o", "--output-dir", type=str,
        metavar="/data/working_directory/project/fastfilter",
        help="Output directory. Created if it does not exist. "
             "Defaults to a 'fastfilter/' subdirectory next to the first input file.",
    )
    parser.add_argument("-j", "--cpus", type=int, default=1, help="Number of parallel CPUs")
    parser.add_argument(
        "-Z", action="store_true",
        help="Use fast compression (level 1) for .gz output. Default is level 6.",
    )
    parser.add_argument(
        "-n", "--max-n", type=int, default=0, metavar="0",
        help="Maximum number of N bases allowed per read. Default: 0 (no Ns tolerated).",
    )
    args = parser.parse_args()

    global min_seq_len, homopolymer_coeff, min_score, max_n
    global r1_files, r2_files, reads_files
    global output_dir, num_cpus, compresslevel

    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.min_score
    max_n = args.max_n
    r1_files = args.r1_files or []
    r2_files = args.r2_files or []
    reads_files = args.reads or []
    output_dir = Path(args.output_dir) if args.output_dir else None
    num_cpus = args.cpus
    compresslevel = 1 if args.Z else 6


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ts() -> str:
    """Current time as [HH:MM:SS] for checkpoint prints."""
    return datetime.now().strftime("[%H:%M:%S]")


def _length_stats(counter: Counter) -> dict:
    """Compute length statistics from a Counter mapping length → read count."""
    if not counter:
        return {"min": 0, "min_count": 0, "max": 0, "max_count": 0,
                "mean": 0.0, "median": 0.0}
    total = sum(counter.values())
    r_min = min(counter)
    r_max = max(counter)
    r_mean = round(sum(length * count for length, count in counter.items()) / total, 2)
    sorted_lens = sorted(counter.keys())
    # lo_idx / hi_idx are 1-based positions of the two middle values.
    # For odd N they are equal (single middle); for even N they differ by 1.
    # Median = (value_at_lo + value_at_hi) / 2  — correct for both cases.
    lo_idx = (total + 1) // 2   # ceiling(N/2)
    hi_idx = (total + 2) // 2   # ceiling((N+1)/2)
    lo_val = hi_val = None
    cumulative = 0
    for length in sorted_lens:
        cumulative += counter[length]
        if lo_val is None and cumulative >= lo_idx:
            lo_val = length
        if cumulative >= hi_idx:
            hi_val = length
            break
    r_median = (lo_val + hi_val) / 2
    return {
        "min": r_min, "min_count": counter[r_min],
        "max": r_max, "max_count": counter[r_max],
        "mean": r_mean, "median": r_median,
    }


def _open_output(path: Path, is_gz: bool):
    """Open an output FASTQ file for writing, with or without gzip."""
    if is_gz:
        return gzip.open(path, "wt", compresslevel=compresslevel)
    return open(path, "w")


# ---------------------------------------------------------------------------
# FASTQ I/O — BioPython-free
# ---------------------------------------------------------------------------

def _fastq_iter(fh):
    """Fast 4-line FASTQ iterator.

    Yields (header, seq, qual) where header is the full first line
    without the leading '@'. Using a local readline binding avoids
    repeated attribute lookups in the hot loop.

    Replaces BioPython's FastqPhredIterator. No format validation is
    performed — assumes well-formed FASTQ from a sequencer.
    """
    readline = fh.readline
    while True:
        header = readline()
        if not header:
            break
        seq  = readline().rstrip().upper()
        readline()              # discard '+' line
        qual = readline().rstrip()
        yield header[1:].rstrip(), seq, qual


def _write_fastq(name: str, seq: str, qual: str, handle) -> None:
    """Write one FASTQ record directly to an open file handle."""
    handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


# ---------------------------------------------------------------------------
# Worker initializer — fixes globals under multiprocessing spawn
# ---------------------------------------------------------------------------

def _init_worker(cfg: dict) -> None:
    """Set filtering globals in each spawned worker process.

    With multiprocessing.spawn, child processes re-import the module from
    scratch and only see module-level defaults — NOT the values set by
    parse_arguments() in the parent. This initializer runs once per worker
    after spawn and restores the correct runtime values.
    """
    global min_seq_len, homopolymer_coeff, min_score, max_n, compresslevel
    global _HOMOPOLYMERS
    min_seq_len       = cfg["min_seq_len"]
    homopolymer_coeff = cfg["homopolymer_coeff"]
    min_score         = cfg["min_score"]
    max_n             = cfg["max_n"]
    compresslevel     = cfg["compresslevel"]
    _HOMOPOLYMERS     = tuple(c * homopolymer_coeff for c in "ATGC")


# ---------------------------------------------------------------------------
# Filtering logic
# ---------------------------------------------------------------------------

def find_homopolymers(seq: str) -> bool:
    """Return True if seq contains a homopolymer run >= homopolymer_coeff."""
    return any(h in seq for h in _HOMOPOLYMERS)


def analyze_sequence(seq: str, qual: str) -> dict:
    """Evaluate one read against all filter criteria.

    Takes raw string fields from _fastq_iter instead of a BioPython SeqRecord.
    Quality score is computed directly from the ASCII quality string:
        mean_phred = mean(ASCII values) - 33
    This avoids building a list of integers per read.
    """
    seq_len = len(seq)

    if seq_len == 0:
        qual_score = 0.0
    else:
        qual_bytes = qual.encode("ascii")
        qual_score = sum(qual_bytes) / seq_len - 33

    n_count   = seq.count("N")
    dot_count = seq.count(".")

    should_include  = True
    too_short       = False
    found_homopolymer = False
    low_score       = False

    if seq_len < min_seq_len:
        should_include = False
        too_short = True
    exceeded_n = n_count > max_n
    if exceeded_n:
        should_include = False
    if dot_count > 0:
        should_include = False
    if find_homopolymers(seq):
        should_include = False
        found_homopolymer = True
    if qual_score < min_score:
        should_include = False
        low_score = True

    return {
        "meets_criteria":    should_include,
        "length":            seq_len,
        "n_count":           n_count,
        "exceeded_n":        exceeded_n,
        "dot_count":         dot_count,
        "too_short":         too_short,
        "found_homopolymer": found_homopolymer,
        "low_score":         low_score,
    }


# ---------------------------------------------------------------------------
# Per-sample workers
# ---------------------------------------------------------------------------

def parse_file_single_end(
    r1_filename: Path,
    output_dir: Path,
    position: int = 0,
    use_tqdm: bool = True,
) -> dict:
    stem1 = r1_filename.stem
    if r1_filename.suffix == ".gz":
        stem1 = Path(stem1).stem
    filename1_no_extension = stem1

    if r1_filename.suffix == ".gz":
        inner_ext = Path(r1_filename.stem).suffix
        filtered1_name = filename1_no_extension + ".filtered" + inner_ext + ".gz"
    else:
        filtered1_name = filename1_no_extension + ".filtered" + r1_filename.suffix

    total_seqs    = 0
    good_sequences = 0
    len_counter1  = Counter()
    exclusion_reasons_r1 = {
        "too_short": 0, "found_dot": 0, "found_n": 0,
        "found_homopolymer": 0, "low_score": 0,
    }

    opener = gzip.open if r1_filename.suffix == ".gz" else open
    is_gz  = r1_filename.suffix == ".gz"

    pbar = tqdm(
        desc=f"Processing {filename1_no_extension}",
        unit=" reads", position=position, leave=True,
        bar_format="{l_bar}{bar} | {n:,} reads [{elapsed}, {rate_fmt}]",
    ) if use_tqdm else None

    with opener(r1_filename, "rt") as f_r1, \
         _open_output(output_dir / filtered1_name, is_gz) as out_f:

        for name, seq, qual in _fastq_iter(f_r1):
            total_seqs += 1
            res = analyze_sequence(seq, qual)
            len_counter1[res["length"]] += 1

            if res["meets_criteria"]:
                _write_fastq(name, seq, qual, out_f)
                good_sequences += 1
            else:
                exclusion_reasons_r1["too_short"]        += bool(res["too_short"])
                exclusion_reasons_r1["found_dot"]        += bool(res["dot_count"])
                exclusion_reasons_r1["found_n"]          += bool(res["exceeded_n"])
                exclusion_reasons_r1["found_homopolymer"]+= bool(res["found_homopolymer"])
                exclusion_reasons_r1["low_score"]        += bool(res["low_score"])

            if pbar:
                pbar.update(1)
            elif total_seqs % LOG_INTERVAL == 0:
                print(f"{_ts()} {filename1_no_extension}: "
                      f"{total_seqs:,} reads processed — {good_sequences:,} passed so far",
                      flush=True)

    if pbar:
        pbar.close()

    print(f"{_ts()} {r1_filename.name}: "
          f"Finished. {good_sequences:,} / {total_seqs:,} passed.", flush=True)

    return {
        "file":               filename1_no_extension,
        "r1_stem":            filename1_no_extension,
        "r1_file":            r1_filename.name,
        "total_seqs":         total_seqs,
        "good_sequences_count": good_sequences,
        "good_reads_r1":      good_sequences,
        "exclusion_reasons_r1": exclusion_reasons_r1,
        "r1_len_stats":       _length_stats(len_counter1),
    }


def parse_file(
    r1_filename: Path,
    r2_filename: Path,
    output_dir: Path,
    position: int = 0,
    use_tqdm: bool = True,
) -> dict:
    """Run the parsing workflow for a pair of R1/R2 FASTQ files."""
    stem1 = r1_filename.stem
    if r1_filename.suffix == ".gz":
        stem1 = Path(stem1).stem
    filename1_no_extension = stem1
    stem2 = r2_filename.stem
    if r2_filename.suffix == ".gz":
        stem2 = Path(stem2).stem
    filename2_no_extension = stem2

    if r1_filename.suffix == ".gz":
        inner_ext1 = Path(r1_filename.stem).suffix
        filtered1_name = filename1_no_extension + ".filtered" + inner_ext1 + ".gz"
    else:
        filtered1_name = filename1_no_extension + ".filtered" + r1_filename.suffix

    if r2_filename.suffix == ".gz":
        inner_ext2 = Path(r2_filename.stem).suffix
        filtered2_name = filename2_no_extension + ".filtered" + inner_ext2 + ".gz"
    else:
        filtered2_name = filename2_no_extension + ".filtered" + r2_filename.suffix

    total_seqs    = 0
    good_reads_r1 = 0
    good_reads_r2 = 0
    good_sequences = 0
    len_counter1  = Counter()
    len_counter2  = Counter()
    exclusion_reasons_r1 = {
        "too_short": 0, "found_dot": 0, "found_n": 0,
        "found_homopolymer": 0, "low_score": 0,
    }
    exclusion_reasons_r2 = {
        "too_short": 0, "found_dot": 0, "found_n": 0,
        "found_homopolymer": 0, "low_score": 0,
    }

    opener_r1 = gzip.open if r1_filename.suffix == ".gz" else open
    opener_r2 = gzip.open if r2_filename.suffix == ".gz" else open
    is_gz1    = r1_filename.suffix == ".gz"
    is_gz2    = r2_filename.suffix == ".gz"

    pbar = tqdm(
        desc=f"Processing {filename1_no_extension}",
        unit=" reads", position=position, leave=True,
        bar_format="{l_bar}{bar} | {n:,} reads [{elapsed}, {rate_fmt}]",
    ) if use_tqdm else None

    with opener_r1(r1_filename, "rt") as f_r1, \
         opener_r2(r2_filename, "rt") as f_r2, \
         _open_output(output_dir / filtered1_name, is_gz1) as out_f1, \
         _open_output(output_dir / filtered2_name, is_gz2) as out_f2:

        # Streaming ID-based sync of R1/R2 pairs
        r1_iter = _fastq_iter(f_r1)
        r2_iter = _fastq_iter(f_r2)
        r1_next = next(r1_iter, None)
        r2_next = next(r2_iter, None)

        while r1_next is not None and r2_next is not None:

            id1 = _ID_SUFFIX_RE.sub('', r1_next[0].split()[0])
            id2 = _ID_SUFFIX_RE.sub('', r2_next[0].split()[0])

            while id1 < id2:
                r1_next = next(r1_iter, None)
                if r1_next is None:
                    break
                id1 = _ID_SUFFIX_RE.sub('', r1_next[0].split()[0])

            while id2 < id1:
                r2_next = next(r2_iter, None)
                if r2_next is None:
                    break
                id2 = _ID_SUFFIX_RE.sub('', r2_next[0].split()[0])

            if r1_next is None or r2_next is None:
                continue

            name1, seq1, qual1 = r1_next
            name2, seq2, qual2 = r2_next
            total_seqs += 1

            res1 = analyze_sequence(seq1, qual1)
            res2 = analyze_sequence(seq2, qual2)

            len_counter1[res1["length"]] += 1
            len_counter2[res2["length"]] += 1

            good_reads_r1 += res1["meets_criteria"]
            good_reads_r2 += res2["meets_criteria"]

            if res1["meets_criteria"] and res2["meets_criteria"]:
                _write_fastq(name1, seq1, qual1, out_f1)
                _write_fastq(name2, seq2, qual2, out_f2)
                good_sequences += 1
            else:
                exclusion_reasons_r1["too_short"]        += bool(res1["too_short"])
                exclusion_reasons_r1["found_dot"]        += bool(res1["dot_count"])
                exclusion_reasons_r1["found_n"]          += bool(res1["exceeded_n"])
                exclusion_reasons_r1["found_homopolymer"]+= bool(res1["found_homopolymer"])
                exclusion_reasons_r1["low_score"]        += bool(res1["low_score"])

                exclusion_reasons_r2["too_short"]        += bool(res2["too_short"])
                exclusion_reasons_r2["found_dot"]        += bool(res2["dot_count"])
                exclusion_reasons_r2["found_n"]          += bool(res2["exceeded_n"])
                exclusion_reasons_r2["found_homopolymer"]+= bool(res2["found_homopolymer"])
                exclusion_reasons_r2["low_score"]        += bool(res2["low_score"])

            r1_next = next(r1_iter, None)
            r2_next = next(r2_iter, None)

            if pbar:
                pbar.update(1)
            elif total_seqs % LOG_INTERVAL == 0:
                print(f"{_ts()} {filename1_no_extension}: "
                      f"{total_seqs:,} reads processed — {good_sequences:,} passed so far",
                      flush=True)

        # Zero cost when counts match (both None). On mismatch, drain the
        # leftover reads to report an exact count, then abort.
        leftover_r1 = 0
        leftover_r2 = 0
        if r1_next is not None:
            leftover_r1 = 1 + sum(1 for _ in r1_iter)
        if r2_next is not None:
            leftover_r2 = 1 + sum(1 for _ in r2_iter)
        if leftover_r1 or leftover_r2:
            raise RuntimeError(
                f"Read count mismatch — {r1_filename.name} had {leftover_r1} extra read(s), "
                f"{r2_filename.name} had {leftover_r2} extra read(s) "
                f"after {total_seqs:,} pairs processed. "
                f"Verify both files are from the same run and neither is truncated."
            )

    if pbar:
        pbar.close()

    print(f"{_ts()} {r1_filename.name}: "
          f"Finished. {good_sequences:,} / {total_seqs:,} passed.", flush=True)

    return {
        "file":               filename1_no_extension.replace("_R1", ""),
        "r1_stem":            filename1_no_extension,
        "r2_stem":            filename2_no_extension,
        "r1_file":            r1_filename.name,
        "r2_file":            r2_filename.name,
        "total_seqs":         total_seqs,
        "good_sequences_count": good_sequences,
        "exclusion_reasons_r1": exclusion_reasons_r1,
        "good_reads_r1":      good_reads_r1,
        "exclusion_reasons_r2": exclusion_reasons_r2,
        "good_reads_r2":      good_reads_r2,
        "r1_len_stats":       _length_stats(len_counter1),
        "r2_len_stats":       _length_stats(len_counter2),
    }


# ---------------------------------------------------------------------------
# Summary output
# ---------------------------------------------------------------------------

def _write_sample_summary(path: Path, header: list, entry: list) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        w.writerow(["metric", "value"])
        for key, val in zip(header, entry):
            w.writerow([key, val])


def generate_summary(result) -> None:
    """Write one <stem>.summary.csv per input file."""
    global single_end, elapsed_min

    shared_tail = ["min_length", "homopolymer_len", "min_score", "max_n_allowed", "elapsed_min"]

    for r in result:
        total  = r["total_seqs"]
        passed = r["good_sequences_count"]
        failed = total - passed
        pct    = round(passed / total * 100, 2) if total > 0 else 0.0
        tail   = [min_seq_len, homopolymer_coeff, min_score, max_n, round(elapsed_min, 2)]

        # --- R1 summary ---
        r1_header = ["sample", "r1_file", "total_reads", "passed_reads", "failed_reads", "pct_pairs_passed"]
        r1_entry  = [r["file"], r["r1_file"], total, passed, failed, pct]

        if not single_end:
            r1_pass_rate        = round(r["good_reads_r1"] / total * 100, 2) if total > 0 else 0.0
            lost_due_to_r1_fail = r["good_reads_r2"] - passed
            failed_both         = total - r["good_reads_r1"] - r["good_reads_r2"] + passed
            r1_header += ["r1_pass_rate", "lost_due_to_r1_fail", "failed_both"]
            r1_entry  += [r1_pass_rate, lost_due_to_r1_fail, failed_both]

        r1_header += [
            "r1_too_short", "r1_n", "r1_dot", "r1_homopolymer", "r1_low_score",
            "r1_len_min", "r1_len_min_count",
            "r1_len_max", "r1_len_max_count",
            "r1_len_mean", "r1_len_median",
        ]
        r1s = r["r1_len_stats"]
        r1_entry += [
            r["exclusion_reasons_r1"]["too_short"],
            r["exclusion_reasons_r1"]["found_n"],
            r["exclusion_reasons_r1"]["found_dot"],
            r["exclusion_reasons_r1"]["found_homopolymer"],
            r["exclusion_reasons_r1"]["low_score"],
            r1s["min"], r1s["min_count"],
            r1s["max"], r1s["max_count"],
            r1s["mean"], r1s["median"],
        ]
        r1_header += shared_tail
        r1_entry  += tail
        _write_sample_summary(output_dir / (r["r1_stem"] + ".summary.csv"), r1_header, r1_entry)

        if single_end:
            continue

        # --- R2 summary ---
        r2_pass_rate        = round(r["good_reads_r2"] / total * 100, 2) if total > 0 else 0.0
        lost_due_to_r2_fail = r["good_reads_r1"] - passed
        failed_both         = total - r["good_reads_r1"] - r["good_reads_r2"] + passed

        r2_header = [
            "sample", "r2_file",
            "total_reads", "passed_reads", "failed_reads", "pct_pairs_passed",
            "r2_pass_rate", "lost_due_to_r2_fail", "failed_both",
            "r2_too_short", "r2_n", "r2_dot", "r2_homopolymer", "r2_low_score",
            "r2_len_min", "r2_len_min_count",
            "r2_len_max", "r2_len_max_count",
            "r2_len_mean", "r2_len_median",
        ] + shared_tail
        r2s = r["r2_len_stats"]
        r2_entry = [
            r["file"], r["r2_file"], total, passed, failed, pct,
            r2_pass_rate, lost_due_to_r2_fail, failed_both,
            r["exclusion_reasons_r2"]["too_short"],
            r["exclusion_reasons_r2"]["found_n"],
            r["exclusion_reasons_r2"]["found_dot"],
            r["exclusion_reasons_r2"]["found_homopolymer"],
            r["exclusion_reasons_r2"]["low_score"],
            r2s["min"], r2s["min_count"],
            r2s["max"], r2s["max_count"],
            r2s["mean"], r2s["median"],
        ] + tail
        _write_sample_summary(output_dir / (r["r2_stem"] + ".summary.csv"), r2_header, r2_entry)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    multiprocessing.set_start_method("spawn")
    start = time.time()
    parse_arguments()

    global output_dir, single_end, elapsed_min, _HOMOPOLYMERS
    _HOMOPOLYMERS = tuple(c * homopolymer_coeff for c in "ATGC")

    # Validate argument values
    if max_n < 0:
        print("Error: --max-n must be >= 0.")
        sys.exit(1)
    if num_cpus < 1:
        print("Error: --cpus must be >= 1.")
        sys.exit(1)

    # Validate inputs and determine mode
    if r1_files and r2_files:
        if len(r1_files) != len(r2_files):
            print(f"Error: {len(r1_files)} -r1 file(s) but {len(r2_files)} -r2 file(s).")
            sys.exit(1)
        first_file = r1_files[0]
        num_samples = len(r1_files)
    elif reads_files:
        single_end = True
        first_file = reads_files[0]
        num_samples = len(reads_files)
    else:
        print("Error: supply -r1/-r2 for paired-end or -r for single-end.")
        sys.exit(1)

    # Check all input files exist before doing anything else
    missing = []
    if not single_end:
        for r1, r2 in zip(r1_files, r2_files):
            if not r1.exists():
                missing.append(str(r1))
            if not r2.exists():
                missing.append(str(r2))
    else:
        for rf in reads_files:
            if not rf.exists():
                missing.append(str(rf))
    if missing:
        print("Error: input file(s) not found:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    if output_dir is None:
        output_dir = first_file.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Cap workers at number of samples; warn if user asked for more
    effective_cpus = min(num_cpus, num_samples)
    if num_cpus > num_samples:
        print(f"Note: {num_cpus} CPUs requested but only {num_samples} sample(s). "
              f"Using {effective_cpus}.")

    # tqdm only when stdout is a live terminal and one worker is running;
    # otherwise use periodic timestamped prints (safe for screen / nohup / logs)
    use_tqdm = sys.stdout.isatty() and effective_cpus == 1

    # Report which gzip backend is active
    gzip_backend = "isal (fast)" if "isal" in gzip.__name__ else "stdlib"

    # Startup summary
    mode_str = "single-end" if single_end else "paired-end"
    print(f"\nfastfilter — {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Output dir  : {output_dir}")
    print(f"  Min length  : {min_seq_len}")
    print(f"  Min score   : {min_score}")
    print(f"  Homopolymer : {homopolymer_coeff}")
    print(f"  Max N       : {max_n}")
    print(f"  Compression : level {compresslevel}  [{gzip_backend}]")
    print(f"  Mode        : {mode_str} | {num_samples} sample(s) | {effective_cpus} CPU(s)")
    print()
    if not single_end:
        for i, (r1, r2) in enumerate(zip(r1_files, r2_files), 1):
            print(f"  [{i}] {r1.name}  +  {r2.name}")
    else:
        for i, r in enumerate(reads_files, 1):
            print(f"  [{i}] {r.name}")
    print()

    worker_cfg = {
        "min_seq_len":       min_seq_len,
        "homopolymer_coeff": homopolymer_coeff,
        "min_score":         min_score,
        "max_n":             max_n,
        "compresslevel":     compresslevel,
    }

    with multiprocessing.Pool(
        processes=effective_cpus,
        initializer=_init_worker,
        initargs=(worker_cfg,),
    ) as pool:
        if not single_end:
            args_for_pool = [
                (r1, r2, output_dir, pos, use_tqdm)
                for pos, (r1, r2) in enumerate(zip(r1_files, r2_files))
            ]
            return_values = pool.starmap(parse_file, args_for_pool)
        else:
            args_for_pool = [
                (f, output_dir, pos, use_tqdm)
                for pos, f in enumerate(reads_files)
            ]
            return_values = pool.starmap(parse_file_single_end, args_for_pool)

    finish = time.time()
    elapsed_min = (finish - start) / 60
    print(f"\n{_ts()} All done. Ran in {elapsed_min:.2f} min.")
    generate_summary(return_values)


if __name__ == "__main__":
    main()
