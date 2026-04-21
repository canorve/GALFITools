#!/usr/bin/env python3
"""
Batch-run GALFIT over a list of input files.

Each non-empty, non-comment line in the list file must contain the path to a
GALFIT input file, for example:

    /obj1/z/galfit.init
    /obj2/z/galfit.init
    /obj3/z/galfit.init

Features
--------
- Serial execution with --jobs 1
- Parallel execution with --jobs N using concurrent.futures
- Safe subprocess execution without shell=True
- Ignores empty lines and comment lines starting with '#'
- Checks that each input file exists
- Continues even if some jobs fail
- Captures return code, stdout, and stderr
- Optional verbose output
- Optional CSV summary report

Exit status
-----------
- 0 if all jobs succeed
- non-zero if any job fails or input validation fails
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Sequence


_PRINT_LOCK = threading.Lock()


@dataclass(frozen=True)
class JobResult:
    """Store the result of one GALFIT execution."""

    input_file: Path
    success: bool
    returncode: int
    stdout: str
    stderr: str
    error_message: str | None = None


def log(message: str) -> None:
    """Print a timestamped log message in a thread-safe way."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with _PRINT_LOCK:
        print(f"[{timestamp}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Batch-run GALFIT over input files listed in a text file."
    )
    parser.add_argument(
        "list_file",
        help="Text file containing one GALFIT input-file path per line.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel workers to use (default: 1).",
    )
    parser.add_argument(
        "--galfit-bin",
        default="galfit",
        help='Path to GALFIT executable (default: "galfit").',
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print stdout/stderr for every GALFIT run.",
    )
    parser.add_argument(
        "--summary-csv",
        type=str,
        default=None,
        help="Write a CSV summary report to this path.",
    )
    return parser.parse_args()


def read_list_file(list_file: Path) -> List[Path]:
    """
    Read the list file and return the input paths.

    Empty lines and lines starting with '#' are ignored.
    Environment variables and '~' are expanded.
    Relative paths are resolved against the current working directory.
    """
    paths: List[Path] = []

    with list_file.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            expanded = os.path.expandvars(os.path.expanduser(line))
            path = Path(expanded)

            if not path.is_absolute():
                path = Path.cwd() / path

            paths.append(path)

    return paths


def validate_input_files(paths: Iterable[Path]) -> tuple[List[Path], List[JobResult]]:
    """
    Split paths into existing files and missing files.

    Returns
    -------
    valid_files
        Existing input files.
    invalid_results
        Failure results for paths that do not exist.
    """
    valid_files: List[Path] = []
    invalid_results: List[JobResult] = []

    for path in paths:
        if path.is_file():
            valid_files.append(path)
        else:
            invalid_results.append(
                JobResult(
                    input_file=path,
                    success=False,
                    returncode=-1,
                    stdout="",
                    stderr="",
                    error_message="Input file does not exist or is not a regular file.",
                )
            )

    return valid_files, invalid_results


def print_streams(result: JobResult) -> None:
    """Print stdout/stderr of one job in a readable block."""
    with _PRINT_LOCK:
        print("=" * 80, flush=True)
        print(f"File: {result.input_file}", flush=True)
        print(f"Return code: {result.returncode}", flush=True)

        if result.error_message:
            print(f"Error: {result.error_message}", flush=True)

        if result.stdout.strip():
            print("---- stdout ----", flush=True)
            print(result.stdout.rstrip(), flush=True)

        if result.stderr.strip():
            print("---- stderr ----", flush=True)
            print(result.stderr.rstrip(), flush=True)

        print("=" * 80, flush=True)


def run_galfit(input_file: Path, galfit_bin: str, verbose: bool = False) -> JobResult:
    """
    Run GALFIT on one input file.

    GALFIT is executed from the directory that contains the input file, and the
    input filename is passed as a basename. This is usually the safest approach
    because GALFIT commonly writes outputs relative to the working directory.
    """
    workdir = input_file.parent
    filename = input_file.name

    log(f"START  | {input_file}")

    try:
        completed = subprocess.run(
            [galfit_bin, filename],
            cwd=str(workdir),
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        result = JobResult(
            input_file=input_file,
            success=False,
            returncode=-1,
            stdout="",
            stderr="",
            error_message=f'GALFIT executable not found: "{galfit_bin}"',
        )
        log(f"FAILED | {input_file} | GALFIT executable not found")
        if verbose:
            print_streams(result)
        return result
    except Exception as exc:
        result = JobResult(
            input_file=input_file,
            success=False,
            returncode=-1,
            stdout="",
            stderr="",
            error_message=f"Unexpected error while running GALFIT: {exc}",
        )
        log(f"FAILED | {input_file} | unexpected exception")
        if verbose:
            print_streams(result)
        return result

    success = completed.returncode == 0

    if success:
        log(f"OK     | {input_file}")
    else:
        log(f"FAILED | {input_file} | return code {completed.returncode}")

    result = JobResult(
        input_file=input_file,
        success=success,
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        error_message=None,
    )

    if verbose:
        print_streams(result)

    return result


def run_serial(
    files: Sequence[Path], galfit_bin: str, verbose: bool
) -> List[JobResult]:
    """Run GALFIT sequentially."""
    results: List[JobResult] = []

    for index, input_file in enumerate(files, start=1):
        log(f"Progress: {index}/{len(files)}")
        result = run_galfit(input_file, galfit_bin, verbose=verbose)
        results.append(result)

    return results


def run_parallel(
    files: Sequence[Path],
    galfit_bin: str,
    jobs: int,
    verbose: bool,
) -> List[JobResult]:
    """Run GALFIT in parallel using a thread pool."""
    results: List[JobResult] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        future_to_file = {
            executor.submit(run_galfit, input_file, galfit_bin, verbose): input_file
            for input_file in files
        }

        completed_count = 0
        total = len(files)

        for future in concurrent.futures.as_completed(future_to_file):
            input_file = future_to_file[future]
            completed_count += 1

            try:
                result = future.result()
            except Exception as exc:
                result = JobResult(
                    input_file=input_file,
                    success=False,
                    returncode=-1,
                    stdout="",
                    stderr="",
                    error_message=f"Worker crashed unexpectedly: {exc}",
                )

            results.append(result)
            log(f"Progress: {completed_count}/{total}")

    return results


def print_failure_details(results: Iterable[JobResult]) -> None:
    """Print detailed information only for failed jobs."""
    failures = [result for result in results if not result.success]

    if not failures:
        return

    print("\nDetailed failure report:", flush=True)
    print("=" * 80, flush=True)

    for result in failures:
        print(f"File       : {result.input_file}", flush=True)
        print(f"Success    : {result.success}", flush=True)
        print(f"Return code: {result.returncode}", flush=True)

        if result.error_message:
            print(f"Error      : {result.error_message}", flush=True)

        if result.stdout.strip():
            print("---- stdout ----", flush=True)
            print(result.stdout.rstrip(), flush=True)

        if result.stderr.strip():
            print("---- stderr ----", flush=True)
            print(result.stderr.rstrip(), flush=True)

        print("-" * 80, flush=True)


def print_summary(results: Sequence[JobResult]) -> None:
    """Print a final execution summary."""
    total = len(results)
    successes = sum(result.success for result in results)
    failures = total - successes

    print("\nSummary", flush=True)
    print("=" * 80, flush=True)
    print(f"Total jobs   : {total}", flush=True)
    print(f"Succeeded    : {successes}", flush=True)
    print(f"Failed       : {failures}", flush=True)

    if failures:
        print("\nFailed files:", flush=True)
        for result in results:
            if not result.success:
                reason = result.error_message or f"return code {result.returncode}"
                print(f"  - {result.input_file} ({reason})", flush=True)


def write_summary_csv(results: Sequence[JobResult], csv_path: Path) -> None:
    """Write a CSV summary report."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["input_file", "success", "returncode", "error_message", "stdout", "stderr"]
        )

        for result in results:
            writer.writerow(
                [
                    str(result.input_file),
                    int(result.success),
                    result.returncode,
                    result.error_message or "",
                    result.stdout,
                    result.stderr,
                ]
            )


def batch_galfit() -> int:
    """Program entry point."""
    args = parse_args()

    list_file = Path(os.path.expandvars(os.path.expanduser(args.list_file)))
    jobs = args.jobs
    galfit_bin = args.galfit_bin
    verbose = args.verbose
    summary_csv = args.summary_csv

    if jobs < 1:
        print("Error: --jobs must be at least 1.", file=sys.stderr)
        return 2

    if not list_file.is_file():
        print(f'Error: list file does not exist: "{list_file}"', file=sys.stderr)
        return 2

    log(f"Reading list file: {list_file}")
    all_paths = read_list_file(list_file)

    if not all_paths:
        print("Error: no valid entries found in the list file.", file=sys.stderr)
        return 2

    valid_files, invalid_results = validate_input_files(all_paths)

    if invalid_results:
        for result in invalid_results:
            log(f"MISSING | {result.input_file}")

    if not valid_files:
        print("Error: none of the listed GALFIT input files exist.", file=sys.stderr)
        print_summary(invalid_results)
        return 1

    log(f"Found {len(valid_files)} existing input file(s).")
    if invalid_results:
        log(f"Found {len(invalid_results)} missing input file(s).")

    if jobs == 1:
        log("Running in sequential mode.")
        run_results = run_serial(valid_files, galfit_bin, verbose=verbose)
    else:
        log(f"Running in parallel mode with {jobs} workers.")
        run_results = run_parallel(valid_files, galfit_bin, jobs=jobs, verbose=verbose)

    results = invalid_results + run_results

    if not verbose:
        print_failure_details(results)

    print_summary(results)

    if summary_csv is not None:
        csv_path = Path(os.path.expandvars(os.path.expanduser(summary_csv)))
        write_summary_csv(results, csv_path)
        log(f"CSV summary written to: {csv_path}")

    all_success = all(result.success for result in results)
    return 0 if all_success else 1
