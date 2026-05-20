#!/usr/bin/env python3
"""
script to process a list of GALFIT files and convert Sersic bar components
to Ferrer components.

This script reads a text file containing one GALFIT file path per line,
changes to the directory of each file, applies `Sersic2Ferrer()` to the file,
stores the results, and writes them to an output CSV file.

The stored path is relative to the directory where the program is launched.

"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from galfitools.galin.sersic2ferrer import Sersic2Ferrer

DEFAULT_ENCODING = "utf-8"
DEFAULT_OUTPUT_FILE = "batchFerrerOut.csv"
DEFAULT_SUFFIX = "_ferrer"
DEFAULT_ALPHA = False
DEFAULT_BETA = False


@dataclass(slots=True)
class Sersic2FerrerResult:
    """Store the processing result for one GALFIT file."""

    relative_input_path: str
    relative_output_path: str | None
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


def make_output_name(file_path: Path, suffix: str) -> str:
    """
    Build the output GALFIT file name for one input file.

    Parameters
    ----------
    file_path : Path
        Input file path.
    suffix : str
        Suffix added before the file extension.

    Returns
    -------
    str
        Output file name.
    """
    return f"{file_path.stem}{suffix}{file_path.suffix}"


def process_file(
    file_path: Path,
    base_dir: Path,
    alpha: bool,
    beta: bool,
    suffix: str,
) -> Sersic2FerrerResult:
    """
    Change to the file directory, process the file, and return the result.

    Parameters
    ----------
    file_path : Path
        GALFIT file to process.
    base_dir : Path
        Directory from which relative paths are computed.
    alpha : bool
        If True, leave the Ferrer alpha parameter free.
    beta : bool
        If True, leave the Ferrer beta parameter free.
    suffix : str
        Suffix used to build the output file name.

    Returns
    -------
    Sersic2FerrerResult
        Result object with status information.
    """
    relative_input_path = make_relative_path(file_path, base_dir)

    if not file_path.exists():
        return Sersic2FerrerResult(
            relative_input_path=relative_input_path,
            relative_output_path=None,
            success=False,
            message="File does not exist.",
        )

    original_cwd = Path.cwd()

    try:
        os.chdir(file_path.parent)

        output_name = make_output_name(file_path, suffix)
        output_path = file_path.parent / output_name

        Sersic2Ferrer(
            galfit_file=file_path.name,
            alpha=alpha,
            beta=beta,
            fileout=output_name,
        )

        return Sersic2FerrerResult(
            relative_input_path=relative_input_path,
            relative_output_path=make_relative_path(output_path, base_dir),
            success=True,
            message="OK",
        )

    except Exception as exc:
        return Sersic2FerrerResult(
            relative_input_path=relative_input_path,
            relative_output_path=None,
            success=False,
            message=f"{type(exc).__name__}: {exc}",
        )

    finally:
        os.chdir(original_cwd)


def process_files(
    file_paths: Iterable[Path],
    base_dir: Path,
    alpha: bool,
    beta: bool,
    suffix: str,
) -> list[Sersic2FerrerResult]:
    """
    Process all files in the input iterable.

    Parameters
    ----------
    file_paths : Iterable[Path]
        Files to process.
    base_dir : Path
        Directory from which relative paths are computed.
    alpha : bool
        If True, leave the Ferrer alpha parameter free.
    beta : bool
        If True, leave the Ferrer beta parameter free.
    suffix : str
        Suffix used to build the output file name.

    Returns
    -------
    list[Sersic2FerrerResult]
        Collected results for all files.
    """
    results: list[Sersic2FerrerResult] = []

    for file_path in file_paths:
        results.append(
            process_file(
                file_path=file_path,
                base_dir=base_dir,
                alpha=alpha,
                beta=beta,
                suffix=suffix,
            )
        )

    return results


def write_results(
    results: Iterable[Sersic2FerrerResult],
    output_file: Path,
) -> None:
    """
    Write processing results to a CSV file.

    Parameters
    ----------
    results : Iterable[Sersic2FerrerResult]
        Results to write.
    output_file : Path
        Destination CSV file.
    """
    with output_file.open("w", newline="", encoding=DEFAULT_ENCODING) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "relative_input_path",
                "relative_output_path",
                "success",
                "message",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.relative_input_path,
                    result.relative_output_path,
                    result.success,
                    result.message,
                ]
            )


def print_summary(results: Iterable[Sersic2FerrerResult]) -> None:
    """
    Print a short summary of the processing results.

    Parameters
    ----------
    results : Iterable[Sersic2FerrerResult]
        Results to summarize.
    """
    results_list = list(results)
    total = len(results_list)
    ok = sum(result.success for result in results_list)
    failed = total - ok

    print(f"Processed files: {total}")
    print(f"Successful: {ok}")
    print(f"Failed: {failed}")


def build_parser() -> argparse.ArgumentParser:
    """
    Build and return the command-line argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Read a list of GALFIT files, convert Sersic bar components to "
            "Ferrer components, and write the results to a CSV file."
        )
    )
    parser.add_argument(
        "list_file",
        type=Path,
        help="Text file containing one GALFIT file path per line.",
    )
    parser.add_argument(
        "--alpha",
        action="store_true",
        help="Leave the Ferrer alpha parameter free.",
    )
    parser.add_argument(
        "--beta",
        action="store_true",
        help="Leave the Ferrer beta parameter free.",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default=DEFAULT_SUFFIX,
        help=f"Suffix added to each output GALFIT file. Default: {DEFAULT_SUFFIX}",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT_FILE),
        help=f"Output CSV file. Default: {DEFAULT_OUTPUT_FILE}",
    )

    return parser


def mainbatchSersic2Ferrer() -> int:
    """
    Run the script.

    Returns
    -------
    int
        Exit status code.
    """
    parser = build_parser()
    args = parser.parse_args()

    base_dir = Path.cwd().resolve()
    list_file = args.list_file.expanduser().resolve()
    output_file = args.output.expanduser().resolve()

    try:
        file_paths = read_file_list(list_file)
        results = process_files(
            file_paths=file_paths,
            base_dir=base_dir,
            alpha=args.alpha,
            beta=args.beta,
            suffix=args.suffix,
        )
        write_results(results, output_file)
        print_summary(results)
        print(f"Results written to: {output_file}")
        return 0

    except Exception as exc:
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(mainSersic2Ferrer())
