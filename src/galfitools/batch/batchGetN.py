#!/usr/bin/env python3
"""
Process a list of GALFIT files and estimate the Sersic index.

This script reads a text file containing one GALFIT file path per line,
changes to the directory of each file, applies `getN()` to the file,
stores the results, and writes them to an output CSV file.

The stored path is relative to the directory where the program is launched.

Python 3.11+
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from galfitools.galout.getN import getN


# ============================================================================
# Configuration
# ============================================================================

DEFAULT_ENCODING = "utf-8"
DEFAULT_OUTPUT_FILE = "sersic_index_results.csv"

DEFAULT_DIS = 10.0
DEFAULT_FRAC = 0.9
DEFAULT_ANGLE: float | None = None
DEFAULT_NUM_COMP = 1
DEFAULT_PLOT = False
DEFAULT_CONST = 0.0

ROUND_DECIMALS = 2


# ============================================================================
# Data structures
# ============================================================================


@dataclass(slots=True)
class NResult:
    """Store the processing result for one GALFIT file."""

    relative_path: str
    sersic: float | None
    mean_ns: float | None
    std_ns: float | None
    totmag: float | None
    n_components: int | None
    theta: float | None
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
    file_path: Path,
    base_dir: Path,
    dis: float,
    frac: float,
    angle: float | None,
    num_comp: int,
    plot: bool,
    const: float,
) -> NResult:
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
    frac : float
        Fraction of light.
    angle : float | None
        Angle of the major axis of the galaxy.
    num_comp : int
        Number of component used to define the center.
    plot : bool
        If True, create the Sersic-index plot.
    const : float
        Constant subtracted from the plot.

    Returns
    -------
    NResult
        Result object with rounded outputs and status information.
    """
    relative_path = make_relative_path(file_path, base_dir)

    if not file_path.exists():
        return NResult(
            relative_path=relative_path,
            sersic=None,
            mean_ns=None,
            std_ns=None,
            totmag=None,
            n_components=None,
            theta=None,
            success=False,
            message="File does not exist.",
        )

    original_cwd = Path.cwd()

    try:
        os.chdir(file_path.parent)

        sersic, mean_ns, std_ns, totmag, n_components, theta = getN(
            galfitFile=file_path.name,
            dis=dis,
            frac=frac,
            angle=angle,
            num_comp=num_comp,
            plot=plot,
            const=const,
        )

        return NResult(
            relative_path=relative_path,
            sersic=round(float(sersic), ROUND_DECIMALS),
            mean_ns=round(float(mean_ns), ROUND_DECIMALS),
            std_ns=round(float(std_ns), ROUND_DECIMALS),
            totmag=round(float(totmag), ROUND_DECIMALS),
            n_components=int(n_components),
            theta=round(float(theta), ROUND_DECIMALS),
            success=True,
            message="OK",
        )

    except Exception as exc:
        return NResult(
            relative_path=relative_path,
            sersic=None,
            mean_ns=None,
            std_ns=None,
            totmag=None,
            n_components=None,
            theta=None,
            success=False,
            message=f"{type(exc).__name__}: {exc}",
        )

    finally:
        os.chdir(original_cwd)


def process_files(
    file_paths: Iterable[Path],
    base_dir: Path,
    dis: float,
    frac: float,
    angle: float | None,
    num_comp: int,
    plot: bool,
    const: float,
) -> list[NResult]:
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
    frac : float
        Fraction of light.
    angle : float | None
        Angle of the major axis of the galaxy.
    num_comp : int
        Number of component used to define the center.
    plot : bool
        If True, create the Sersic-index plot.
    const : float
        Constant subtracted from the plot.

    Returns
    -------
    list[NResult]
        Collected results for all files.
    """
    results: list[NResult] = []

    for file_path in file_paths:
        results.append(
            process_file(
                file_path=file_path,
                base_dir=base_dir,
                dis=dis,
                frac=frac,
                angle=angle,
                num_comp=num_comp,
                plot=plot,
                const=const,
            )
        )

    return results


def write_results(results: Iterable[NResult], output_file: Path) -> None:
    """
    Write processing results to a CSV file.

    Parameters
    ----------
    results : Iterable[NResult]
        Results to write.
    output_file : Path
        Destination CSV file.
    """
    with output_file.open("w", newline="", encoding=DEFAULT_ENCODING) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "relative_path",
                "sersic",
                "mean_ns",
                "std_ns",
                "totmag",
                "n_components",
                "theta",
                "success",
                "message",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.relative_path,
                    result.sersic,
                    result.mean_ns,
                    result.std_ns,
                    result.totmag,
                    result.n_components,
                    result.theta,
                    result.success,
                    result.message,
                ]
            )


def print_summary(results: Iterable[NResult]) -> None:
    """
    Print a short summary of the processing results.

    Parameters
    ----------
    results : Iterable[NResult]
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
            "Read a list of GALFIT files, estimate Sersic-index quantities "
            "for each file, and write the results to a CSV file."
        )
    )
    parser.add_argument(
        "list_file",
        type=Path,
        help="Text file containing one GALFIT file path per line.",
    )
    parser.add_argument(
        "-d",
        "--dis",
        type=float,
        default=3,
        help=f"Maximum distance among components. Default: 3 ",
    )
    parser.add_argument(
        "-f",
        "--frac",
        type=float,
        default=0.5,
        help=f"Fraction of light. Default: .5 ",
    )
    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        default=None,
        help=(
            "Angle of the major axis of the galaxy. "
            "If omitted, the angle of the last component is used."
        ),
    )
    parser.add_argument(
        "-n",
        "--num_comp",
        type=int,
        default=1,
        help=f"Component number used to define the center. Default: 1",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="Create the Sersic-index plot.",
    )
    parser.add_argument(
        "--const",
        type=float,
        default=0,
        help=f"Constant subtracted from the plot. Default: 0",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default="getn.csv",
        help=f"Output CSV file. Default: getn.csv",
    )

    return parser


def mainbatchGetN() -> int:
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
    output_file = (
        args.output.expanduser().resolve()
        if args.output is not None
        else (base_dir / DEFAULT_OUTPUT_FILE)
    )

    try:
        file_paths = read_file_list(list_file)
        results = process_files(
            file_paths=file_paths,
            base_dir=base_dir,
            dis=args.dis,
            frac=args.frac,
            angle=args.angle,
            num_comp=args.num_comp,
            plot=args.plot,
            const=args.const,
        )
        write_results(results, output_file)
        print_summary(results)
        print(f"Results written to: {output_file}")
        return 0

    except Exception as exc:
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


# if __name__ == "__main__":
#    raise SystemExit(main())
