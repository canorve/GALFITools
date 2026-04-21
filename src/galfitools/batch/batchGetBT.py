#!/usr/bin/env python3
"""
Template to process a list of GALFIT files and compute B/T quantities.

This script reads a text file containing one GALFIT file path per line,
changes to the directory of each file, applies `get_bt()` to the file,
stores the results, and writes them to an output CSV file.

The stored path is relative to the directory where the program is launched.

Python 3.11+
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import csv
import os
import sys
import argparse

from galfitools.galout.getBT import getBT


# ============================================================================
# Configuration
# ============================================================================

DEFAULT_ENCODING = "utf-8"
DEFAULT_OUTPUT_FILE = "bt_results.csv"
ROUND_DECIMALS = 2


# ============================================================================
# Data structures
# ============================================================================


@dataclass(slots=True)
class BTResult:
    """Store the processing result for one GALFIT file."""

    relative_path: str
    bulge_total: float | None
    totmag: float | None
    n_components: int | None
    success: bool
    message: str = ""


def read_file_list(list_file: Path, encoding: str = DEFAULT_ENCODING) -> list[Path]:
    """
    Read a text file containing one file path per line.

    Empty lines and lines starting with '#' are ignored.

    Parameters
    ----------
    list_file : Path
        Path to the file containing the list of input file paths.
    encoding : str, optional
        File encoding.

    Returns
    -------
    list[Path]
        List of resolved input paths.
    """
    if not list_file.is_file():
        raise FileNotFoundError(f"List file not found: {list_file}")

    paths: list[Path] = []

    with list_file.open("r", encoding=encoding) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            path = Path(line).expanduser()

            if not path.is_absolute():
                path = (list_file.parent / path).resolve()

            paths.append(path)

    if not paths:
        raise ValueError(f"No valid file paths found in: {list_file}")

    return paths


def make_relative_path(file_path: Path, base_dir: Path) -> str:
    """
    Return the file path relative to the base directory.

    If the file is not inside the base directory, return the full path.

    Parameters
    ----------
    file_path : Path
        Absolute path to the file.
    base_dir : Path
        Base directory used to compute the relative path.

    Returns
    -------
    str
        Relative path if possible, otherwise absolute path as string.
    """
    try:
        return str(file_path.relative_to(base_dir))
    except ValueError:
        return str(file_path)


def process_file(
    file_path: Path, base_dir: Path, dis: float, num_comp: int
) -> BTResult:
    """
    Change to the file directory, process the file, and return the result.

    Parameters
    ----------
    file_path : Path
        GALFIT file to process.
    base_dir : Path
        Directory from which relative paths are computed.
    dis : float
        Maximum distance among components.
    num_comp : int
        Number of component where the center will be taken.

    Returns
    -------
    BTResult
        Result object with rounded outputs and status information.
    """
    relative_path = make_relative_path(file_path, base_dir)

    if not file_path.exists():
        return BTResult(
            relative_path=relative_path,
            bulge_total=None,
            totmag=None,
            n_components=None,
            success=False,
            message="File does not exist.",
        )

    original_cwd = Path.cwd()

    try:
        os.chdir(file_path.parent)

        bulge_total, totmag, n_components = getBT(file_path.name, dis, num_comp)

        return BTResult(
            relative_path=relative_path,
            bulge_total=round(float(bulge_total), ROUND_DECIMALS),
            totmag=round(float(totmag), ROUND_DECIMALS),
            n_components=int(n_components),
            success=True,
            message="OK",
        )

    except Exception as exc:
        return BTResult(
            relative_path=relative_path,
            bulge_total=None,
            totmag=None,
            n_components=None,
            success=False,
            message=f"{type(exc).__name__}: {exc}",
        )

    finally:
        os.chdir(original_cwd)


def process_files(
    file_paths: Iterable[Path],
    base_dir: Path,
    dis: float,
    num_comp: int,
) -> list[BTResult]:
    """
    Process all files in the input iterable.

    Parameters
    ----------
    file_paths : Iterable[Path]
        Files to process.
    base_dir : Path
        Directory from which relative paths are computed.
    dis : float
        Maximum distance among components.
    num_comp : int
        Number of component where the center will be taken.

    Returns
    -------
    list[BTResult]
        Collected results for all files.
    """
    results: list[BTResult] = []

    for file_path in file_paths:
        results.append(process_file(file_path, base_dir, dis, num_comp))

    return results


def write_results(results: Iterable[BTResult], output_file: Path) -> None:
    """
    Write processing results to a CSV file.

    Parameters
    ----------
    results : Iterable[BTResult]
        Results to write.
    output_file : Path
        Destination CSV file.
    """
    with output_file.open("w", newline="", encoding=DEFAULT_ENCODING) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "relative_path",
                "bulge_total",
                "totmag",
                "n_components",
                "success",
                "message",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.relative_path,
                    result.bulge_total,
                    result.totmag,
                    result.n_components,
                    result.success,
                    result.message,
                ]
            )


def print_summary(results: Iterable[BTResult]) -> None:
    """
    Print a short summary of the processing results.

    Parameters
    ----------
    results : Iterable[BTResult]
        Results to summarize.
    """
    results_list = list(results)
    total = len(results_list)
    ok = sum(result.success for result in results_list)
    failed = total - ok

    print(f"Processed files: {total}")
    print(f"Successful: {ok}")
    print(f"Failed: {failed}")


def mainGetBT() -> int:
    """
    Run the script.

    Returns
    -------
    int
        Exit status code.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Read a list of GALFIT files, compute B/T quantities for each file, "
            "and write the results to a CSV file."
        )
    )
    parser.add_argument(
        "list_file",
        type=Path,
        help="Text file containing one GALFIT file path per line.",
    )
    parser.add_argument(
        "-n",
        "--num_comp",
        default=1,
        type=int,
        help="Component number used to define the center.",
    )
    parser.add_argument(
        "-d",
        "--dis",
        type=float,
        default=3,
        help="Maximum distance among components.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=f"Output CSV file. Default: {DEFAULT_OUTPUT_FILE}",
    )

    args = parser.parse_args()

    base_dir = Path.cwd().resolve()
    list_file = args.list_file.expanduser().resolve()
    dis = args.dis
    num_comp = args.num_comp
    output_file = (
        args.output.expanduser().resolve()
        if args.output is not None
        else (base_dir / DEFAULT_OUTPUT_FILE)
    )

    try:
        file_paths = read_file_list(list_file)
        results = process_files(file_paths, base_dir, dis, num_comp)
        write_results(results, output_file)
        print_summary(results)
        print(f"Results written to: {output_file}")
        return 0

    except Exception as exc:
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


# if __name__ == "__main__":
#    raise SystemExit(main())
