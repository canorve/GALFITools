#!/usr/bin/env python3
"""
Template to process a list of GALFIT files and compute the effective radius
or another light-fraction radius.

This script reads a text file containing one GALFIT file path per line,
changes to the directory of each file, applies getReComp to the file,
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

from galfitools.galout.getRads import getReComp


DEFAULT_ENCODING = "utf-8"
DEFAULT_OUTPUT_FILE = "out.csv"
ROUND_DECIMALS = 2


@dataclass(slots=True)
class ReCompResult:
    """Store the processing result for one GALFIT file."""

    relative_path: str
    effrad: float | None
    totmag: float | None
    meanme: float | None
    me: float | None
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
    eff: float,
    angle: float | None,
    num_comp: int,
    mecorr: float,
) -> ReCompResult:
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
    eff : float
        Fraction of total light.
    angle : float | None
        Position angle of the major axis of the galaxy.
    num_comp : int
        Number of component used to define the center.
    mecorr : float
        Surface-brightness correction for universe expansion.

    Returns
    -------
    ReCompResult
        Result object with rounded outputs and status information.
    """
    relative_path = make_relative_path(file_path, base_dir)

    if not file_path.exists():
        return ReCompResult(
            relative_path=relative_path,
            effrad=None,
            totmag=None,
            meanme=None,
            me=None,
            n_components=None,
            theta=None,
            success=False,
            message="File does not exist.",
        )

    original_cwd = Path.cwd()

    try:
        os.chdir(file_path.parent)

        effrad, totmag, meanme, me, n_components, theta = getReComp(
            galfit_file=file_path.name,
            dis=dis,
            eff=eff,
            angle=angle,
            num_comp=num_comp,
            mecorr=mecorr,
        )

        return ReCompResult(
            relative_path=relative_path,
            effrad=round(float(effrad), ROUND_DECIMALS),
            totmag=round(float(totmag), ROUND_DECIMALS),
            meanme=round(float(meanme), ROUND_DECIMALS),
            me=round(float(me), ROUND_DECIMALS),
            n_components=int(n_components),
            theta=round(float(theta), ROUND_DECIMALS),
            success=True,
            message="OK",
        )

    except Exception as exc:
        return ReCompResult(
            relative_path=relative_path,
            effrad=None,
            totmag=None,
            meanme=None,
            me=None,
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
    eff: float,
    angle: float | None,
    num_comp: int,
    mecorr: float,
) -> list[ReCompResult]:
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
    eff : float
        Fraction of total light.
    angle : float | None
        Position angle of the major axis of the galaxy.
    num_comp : int
        Number of component used to define the center.
    mecorr : float
        Surface-brightness correction for universe expansion.

    Returns
    -------
    list[ReCompResult]
        Collected results for all files.
    """
    results: list[ReCompResult] = []

    for file_path in file_paths:
        results.append(
            process_file(
                file_path=file_path,
                base_dir=base_dir,
                dis=dis,
                eff=eff,
                angle=angle,
                num_comp=num_comp,
                mecorr=mecorr,
            )
        )

    return results


def write_results(results: Iterable[ReCompResult], output_file: Path) -> None:
    """
    Write processing results to a CSV file.

    Parameters
    ----------
    results : Iterable[ReCompResult]
        Results to write.
    output_file : Path
        Destination CSV file.
    """
    with output_file.open("w", newline="", encoding=DEFAULT_ENCODING) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "relative_path",
                "effrad",
                "totmag",
                "meanme",
                "me",
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
                    result.effrad,
                    result.totmag,
                    result.meanme,
                    result.me,
                    result.n_components,
                    result.theta,
                    result.success,
                    result.message,
                ]
            )


def print_summary(results: Iterable[ReCompResult]) -> None:
    """
    Print a short summary of the processing results.

    Parameters
    ----------
    results : Iterable[ReCompResult]
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
            "Read a list of GALFIT files, compute the effective radius or "
            "another light-fraction radius, and write the results to a CSV file."
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
        help="Maximum distance among components. Default: 3",
    )
    parser.add_argument(
        "--eff",
        type=float,
        default=0.5,
        help="Fraction of total light. Must be between 0 and 1. Default: 0.5",
    )
    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        default=None,
        help=(
            "Position angle of the major axis of the galaxy. "
            "If omitted, the angle of the last component is used."
        ),
    )
    parser.add_argument(
        "-n",
        "--num_comp",
        type=int,
        default=1,
        help="Component number used to define the center. Default: 1",
    )
    parser.add_argument(
        "--mecorr",
        type=float,
        default=0.0,
        help="Surface-brightness correction for universe expansion. Default: 0.0",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("out.csv"),
        help="Output CSV file. Default: out.csv",
    )

    return parser


def mainbathcGetReComp() -> int:
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

    if not (0 < args.eff <= 1):
        print("Error: --eff must be a value between 0 and 1.", file=sys.stderr)
        return 1

    try:
        file_paths = read_file_list(list_file)
        results = process_files(
            file_paths=file_paths,
            base_dir=base_dir,
            dis=args.dis,
            eff=args.eff,
            angle=args.angle,
            num_comp=args.num_comp,
            mecorr=args.mecorr,
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
