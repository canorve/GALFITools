#!/usr/bin/env python3
"""
Batch-run GALFIT over input files listed in a text file.

Files located in the same directory are always executed sequentially.
Files located in different directories may be executed concurrently.

Each non-empty, non-comment line must contain the path to a GALFIT input file:

    /obj1/z/galfit.init
    /obj1/z/galfit2.init
    /obj2/z/galfit.init

Features
--------
- Serial execution with --jobs 1
- Parallel execution with --jobs N
- Files in the same directory never run simultaneously
- Safe subprocess execution without shell=True
- Ignores empty lines and lines beginning with '#'
- Checks whether every input file exists
- Continues processing after failures
- Captures return code, stdout, and stderr
- Optional verbose output
- Optional CSV summary

Exit status
-----------
- 0 if every job succeeds
- 1 if one or more jobs fail
- 2 for command-line or input-list errors
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence


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


@dataclass
class ProgressTracker:
    """Maintain a thread-safe count of completed jobs."""

    total: int
    completed: int = 0
    lock: threading.Lock = field(
        default_factory=threading.Lock,
        repr=False,
    )

    def advance(self, input_file: Path) -> None:
        """Increment and print the completed-job count."""
        with self.lock:
            self.completed += 1
            completed = self.completed

        log(f"Progress: {completed}/{self.total} completed | {input_file}")


def log(message: str) -> None:
    """Print a timestamped message in a thread-safe way."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with _PRINT_LOCK:
        print(f"[{timestamp}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Batch-run GALFIT over files listed in a text file. "
            "Files in the same directory are always run sequentially."
        )
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
        help=("Maximum number of directories processed concurrently " "(default: 1)."),
    )

    parser.add_argument(
        "--galfit-bin",
        default="galfit",
        help='GALFIT executable or path to it (default: "galfit").',
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print stdout and stderr for every GALFIT execution.",
    )

    parser.add_argument(
        "--summary-csv",
        type=str,
        default=None,
        help="Write a CSV execution report to the specified file.",
    )

    return parser.parse_args()


def expand_path(path_text: str) -> Path:
    """Expand environment variables and the user home directory."""
    expanded = os.path.expandvars(os.path.expanduser(path_text))
    return Path(expanded)


def normalize_executable(executable: str) -> str:
    """
    Normalize an executable path.

    Executable names such as 'galfit' are left unchanged so they can be found
    through PATH. Explicit relative paths such as './bin/galfit' are converted
    to absolute paths because subprocesses run from each input directory.
    """
    expanded = os.path.expandvars(os.path.expanduser(executable))

    contains_separator = os.sep in expanded

    if os.altsep is not None:
        contains_separator = contains_separator or os.altsep in expanded

    if not contains_separator:
        return expanded

    executable_path = Path(expanded)

    if not executable_path.is_absolute():
        executable_path = Path.cwd() / executable_path

    return str(executable_path.resolve())


def read_list_file(list_file: Path) -> list[Path]:
    """
    Read GALFIT input paths from a text file.

    Empty lines and lines beginning with '#' are ignored. Relative paths are
    interpreted relative to the current working directory.
    """
    paths: list[Path] = []

    with list_file.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            path = expand_path(line)

            if not path.is_absolute():
                path = Path.cwd() / path

            paths.append(path)

    return paths


def validate_input_files(
    paths: Iterable[Path],
) -> tuple[list[Path], list[JobResult]]:
    """
    Separate valid files from missing or invalid paths.

    Returns
    -------
    valid_files
        Input paths that exist and are regular files.
    invalid_results
        Failure results corresponding to invalid paths.
    """
    valid_files: list[Path] = []
    invalid_results: list[JobResult] = []

    for path in paths:
        if path.is_file():
            valid_files.append(path)
            continue

        invalid_results.append(
            JobResult(
                input_file=path,
                success=False,
                returncode=-1,
                stdout="",
                stderr="",
                error_message=("Input file does not exist or is not a regular file."),
            )
        )

    return valid_files, invalid_results


def group_files_by_directory(
    files: Sequence[Path],
) -> dict[Path, list[Path]]:
    """
    Group input files according to their resolved parent directory.

    The file order within each directory follows the order in the list file.
    """
    groups: dict[Path, list[Path]] = {}

    for input_file in files:
        directory = input_file.parent.resolve()
        groups.setdefault(directory, []).append(input_file)

    return groups


def print_streams(result: JobResult) -> None:
    """Print captured output for one GALFIT job."""
    with _PRINT_LOCK:
        print("=" * 80, flush=True)
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

        print("=" * 80, flush=True)


def run_galfit(
    input_file: Path,
    galfit_bin: str,
    verbose: bool = False,
) -> JobResult:
    """
    Run GALFIT on one input file.

    GALFIT is executed from the directory containing the input file. Only the
    input filename is passed to GALFIT:

        galfit galfit.init

    Parameters
    ----------
    input_file
        GALFIT initial-parameters file.
    galfit_bin
        GALFIT executable name or absolute path.
    verbose
        Print captured stdout and stderr when True.

    Returns
    -------
    JobResult
        Execution status and captured output.
    """
    working_directory = input_file.parent
    input_filename = input_file.name

    log(f"START  | {input_file}")

    try:
        completed = subprocess.run(
            [galfit_bin, input_filename],
            cwd=working_directory,
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

        log(f"FAILED | {input_file} | executable not found")

        if verbose:
            print_streams(result)

        return result

    except PermissionError:
        result = JobResult(
            input_file=input_file,
            success=False,
            returncode=-1,
            stdout="",
            stderr="",
            error_message=(f'Permission denied when executing GALFIT: "{galfit_bin}"'),
        )

        log(f"FAILED | {input_file} | permission denied")

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
            error_message=f"Unexpected execution error: {exc}",
        )

        log(f"FAILED | {input_file} | unexpected error")

        if verbose:
            print_streams(result)

        return result

    success = completed.returncode == 0

    result = JobResult(
        input_file=input_file,
        success=success,
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        error_message=None,
    )

    if success:
        log(f"OK     | {input_file}")
    else:
        log(f"FAILED | {input_file} | " f"return code {completed.returncode}")

    if verbose:
        print_streams(result)

    return result


def run_directory_group(
    directory: Path,
    files: Sequence[Path],
    galfit_bin: str,
    verbose: bool,
    progress: ProgressTracker,
) -> list[JobResult]:
    """
    Run all GALFIT files from one directory sequentially.

    This function is submitted as one executor task. Therefore, no two input
    files from the same directory can run simultaneously.
    """
    results: list[JobResult] = []

    log(f"DIRECTORY START | {directory} | " f"{len(files)} input file(s)")

    for input_file in files:
        try:
            result = run_galfit(
                input_file=input_file,
                galfit_bin=galfit_bin,
                verbose=verbose,
            )
        except Exception as exc:
            # Defensive fallback. run_galfit already handles expected errors.
            result = JobResult(
                input_file=input_file,
                success=False,
                returncode=-1,
                stdout="",
                stderr="",
                error_message=f"Unexpected worker error: {exc}",
            )

            log(f"FAILED | {input_file} | worker error")

        results.append(result)
        progress.advance(input_file)

    log(f"DIRECTORY DONE  | {directory}")

    return results


def run_serial(
    files: Sequence[Path],
    galfit_bin: str,
    verbose: bool,
) -> list[JobResult]:
    """Run every GALFIT input file sequentially."""
    results: list[JobResult] = []
    progress = ProgressTracker(total=len(files))

    for input_file in files:
        result = run_galfit(
            input_file=input_file,
            galfit_bin=galfit_bin,
            verbose=verbose,
        )

        results.append(result)
        progress.advance(input_file)

    return results


def run_parallel(
    files: Sequence[Path],
    galfit_bin: str,
    jobs: int,
    verbose: bool,
) -> list[JobResult]:
    """
    Run directories concurrently while serializing files within each directory.

    The effective parallelism cannot exceed the number of unique directories.
    """
    grouped_files = group_files_by_directory(files)
    progress = ProgressTracker(total=len(files))
    results: list[JobResult] = []

    effective_workers = min(jobs, len(grouped_files))

    log(f"Found {len(grouped_files)} unique input directories.")
    log(f"Using {effective_workers} concurrent directory worker(s).")

    with concurrent.futures.ThreadPoolExecutor(
        max_workers=effective_workers
    ) as executor:
        future_to_group = {
            executor.submit(
                run_directory_group,
                directory,
                directory_files,
                galfit_bin,
                verbose,
                progress,
            ): (directory, directory_files)
            for directory, directory_files in grouped_files.items()
        }

        for future in concurrent.futures.as_completed(future_to_group):
            directory, directory_files = future_to_group[future]

            try:
                directory_results = future.result()
            except Exception as exc:
                log(
                    f"DIRECTORY FAILED | {directory} | "
                    f"unexpected worker failure: {exc}"
                )

                directory_results = [
                    JobResult(
                        input_file=input_file,
                        success=False,
                        returncode=-1,
                        stdout="",
                        stderr="",
                        error_message=(
                            "Directory worker failed unexpectedly: " f"{exc}"
                        ),
                    )
                    for input_file in directory_files
                ]

            results.extend(directory_results)

    return results


def print_failure_details(results: Iterable[JobResult]) -> None:
    """Print detailed output for failed jobs."""
    failures = [result for result in results if not result.success]

    if not failures:
        return

    print("\nDetailed failure report", flush=True)
    print("=" * 80, flush=True)

    for result in failures:
        print(f"File       : {result.input_file}", flush=True)
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
    """Print the final execution summary."""
    total = len(results)
    succeeded = sum(result.success for result in results)
    failed = total - succeeded

    print("\nSummary", flush=True)
    print("=" * 80, flush=True)
    print(f"Total jobs: {total}", flush=True)
    print(f"Succeeded : {succeeded}", flush=True)
    print(f"Failed    : {failed}", flush=True)

    if failed:
        print("\nFailed files:", flush=True)

        for result in results:
            if result.success:
                continue

            reason = result.error_message or f"return code {result.returncode}"

            print(
                f"  - {result.input_file} ({reason})",
                flush=True,
            )


def write_summary_csv(
    results: Sequence[JobResult],
    csv_path: Path,
) -> None:
    """Write execution results to a CSV file."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle)

        writer.writerow(
            [
                "input_file",
                "directory",
                "success",
                "returncode",
                "error_message",
                "stdout",
                "stderr",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    str(result.input_file),
                    str(result.input_file.parent),
                    int(result.success),
                    result.returncode,
                    result.error_message or "",
                    result.stdout,
                    result.stderr,
                ]
            )


def batch_galfit() -> int:
    """Run the command-line program."""
    args = parse_args()

    if args.jobs < 1:
        print(
            "Error: --jobs must be at least 1.",
            file=sys.stderr,
        )
        return 2

    list_file = expand_path(args.list_file)

    if not list_file.is_absolute():
        list_file = Path.cwd() / list_file

    if not list_file.is_file():
        print(
            f'Error: list file does not exist: "{list_file}"',
            file=sys.stderr,
        )
        return 2

    galfit_bin = normalize_executable(args.galfit_bin)

    log(f"Reading list file: {list_file}")

    try:
        listed_paths = read_list_file(list_file)
    except OSError as exc:
        print(
            f"Error reading list file: {exc}",
            file=sys.stderr,
        )
        return 2

    if not listed_paths:
        print(
            "Error: no valid entries were found in the list file.",
            file=sys.stderr,
        )
        return 2

    valid_files, invalid_results = validate_input_files(listed_paths)

    for result in invalid_results:
        log(f"MISSING | {result.input_file}")

    if not valid_files:
        print(
            "Error: none of the listed GALFIT input files exist.",
            file=sys.stderr,
        )
        print_summary(invalid_results)
        return 1

    log(f"Existing input files: {len(valid_files)}")

    if invalid_results:
        log(f"Missing input files: {len(invalid_results)}")

    if args.jobs == 1:
        log("Running in sequential mode.")

        execution_results = run_serial(
            files=valid_files,
            galfit_bin=galfit_bin,
            verbose=args.verbose,
        )
    else:
        log(f"Running in parallel mode with up to " f"{args.jobs} workers.")

        execution_results = run_parallel(
            files=valid_files,
            galfit_bin=galfit_bin,
            jobs=args.jobs,
            verbose=args.verbose,
        )

    results = invalid_results + execution_results

    if not args.verbose:
        print_failure_details(results)

    print_summary(results)

    if args.summary_csv is not None:
        summary_path = expand_path(args.summary_csv)

        if not summary_path.is_absolute():
            summary_path = Path.cwd() / summary_path

        try:
            write_summary_csv(results, summary_path)
        except OSError as exc:
            print(
                f"Error writing CSV summary: {exc}",
                file=sys.stderr,
            )
            return 1

        log(f"CSV summary written to: {summary_path}")

    return 0 if all(result.success for result in results) else 1


if __name__ == "__main__":
    sys.exit(main())
