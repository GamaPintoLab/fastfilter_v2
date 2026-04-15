# Changelog

All notable changes to fastfilter are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [2.0.0] — 2026

Complete rewrite of the original `fastfilter.py`. All filtering logic preserved and validated
against the original — read counts match exactly on all tested datasets.

### Added
- Explicit `-r1` / `-r2` file arguments for paired-end mode replacing directory-based auto-discovery (`-i`)
- `-r` argument for single-end mode
- `-n` / `--max-n` flag: configurable maximum N bases per read (default `0`, preserving original behaviour)
- Per-sample `.summary.csv` in vertical `metric,value` format for readability
- Summary fields: `pct_pairs_passed`, `r1_pass_rate`, `r2_pass_rate`, `lost_due_to_r1_fail`, `lost_due_to_r2_fail`, `failed_both`, read length min/max/mean/median
- All filter parameters recorded in every summary file (`min_length`, `homopolymer_len`, `min_score`, `max_n_allowed`, `elapsed_min`)
- Startup summary block: mode, thresholds, sample list, gzip backend, CPU count
- Timestamped checkpoint prints every 5 M reads for screen / nohup sessions
- tqdm progress bar for interactive single-worker runs
- Runtime read-count mismatch detection: aborts with a clear error if R1 and R2 file lengths differ
- Input file existence check before any processing begins
- Argument validation: `--max-n >= 0`, `--cpus >= 1`
- gzip backend detection: uses [isal](https://github.com/pycompression/python-isal) (Intel ISA-L, 2–4× faster) with silent stdlib fallback
- `-Z` flag: compression level 1 for fast output

### Changed
- **BioPython removed** — replaced `FastqPhredIterator` with a custom 4-line FASTQ iterator (`_fastq_iter`) and `SeqIO.write` with a direct string writer (`_write_fastq`)
- **Pandas removed** — read length statistics use `collections.Counter` instead of `pd.Series`; memory reduced from ~4.7 GB to kilobytes on 167 M reads
- Homopolymer strings precomputed once per process in `_init_worker` — eliminates 4 string allocations per read call
- Sequence bases uppercased at read time (normalisation)
- Median calculation corrected for even-N datasets (averages the two middle values)
- Worker globals fixed for `multiprocessing.spawn`: `_init_worker` restores runtime thresholds in each spawned worker (previously workers used module defaults regardless of CLI flags)
- Output files use `"w"` mode — re-runs overwrite cleanly
- Pool wrapped in context manager for guaranteed cleanup
- `rstrip("/12")` ID stripping replaced by compiled `re.compile(r'/[12]$')` regex
- `failed_r1_only` / `failed_r2_only` renamed to `lost_due_to_r1_fail` / `lost_due_to_r2_fail` for clarity
- N-count exclusion reason now correctly reflects `n_count > max_n` (not `n_count > 0`), preventing false counts when `--max-n > 0`
- Worker pool size capped at number of samples; excess CPU warning printed
- Shebang changed to `#!/usr/bin/env python3` for portability

### Removed
- `-i` / `--sequences-dir` directory input and auto-discovery logic
- `inquirer` interactive prompt dependency
- BioPython dependency (`biopython`)
- Pandas dependency (`pandas`)
- `dryrun` flag and logic
- `query_seqDir()` function
- Dead code: unreachable `id1 != id2` branch in the sync guard

---

## [1.0.0] — 2022

Initial release (`fastfilter.py`).

- Paired-end FASTQ filtering via `-i` input directory
- Auto-discovery of `_R1_` / `_R2_` file pairs
- Filters: minimum length, mean Phred quality score, homopolymer runs, N bases, dot characters
- BioPython (`FastqPhredIterator`, `SeqIO.write`) for FASTQ I/O
- Pandas for read length statistics
- Per-sample CSV and overview summary
- Multiprocessing support via `multiprocessing.Pool`
