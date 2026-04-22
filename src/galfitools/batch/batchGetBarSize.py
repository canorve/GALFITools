#! /usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional
import csv
import os
import sys
import argparse

from galfitools.galout.getBarSize import getBarSize

from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import getBreak2, getKappa2

# ============================================================================
# Configuration
# ============================================================================

DEFAULT_ENCODING = "utf-8"
DEFAULT_OUTPUT_FILE = "bar_sizes.csv"
BAR_SIZE_DECIMALS = 2

# ============================================================================
# Data structures
# ============================================================================


@dataclass
class BarResult:
    """Store the processing result for one input file."""

    relative_path: Path
    input_galfit: Path
    bar_size: float | None  # pixels
    bar_size_arcsec: float | None
    success: bool
    message: str = ""


# ============================================================================
# Core domain logic
# ============================================================================


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

    Raises
    ------
    FileNotFoundError
        If the list file does not exist.
    ValueError
        If no valid file paths are found.
    """

    if not list_file.is_file():
        raise FileNotFoundError(f"List file not found: {list_file}")

    paths: list[Path] = []

    with list_file.open("r", encoding=encoding) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            path = Path(line).expanduser()

            # Resolve relative paths with respect to the list file location.
            if not path.is_absolute():
                path = (list_file.parent / path).resolve()

            paths.append(path)

    if not paths:
        raise ValueError(f"No valid file paths found in: {list_file}")

    return paths


def process_file(
    file_path: Path, base_dir: Path, dis, numcomp, plot, ranx, out, red
) -> BarResult:
    """
    Change to the file directory, process the file, and return the result.

    Parameters
    ----------
    file_path : Path
        File to process.

    base_dir : Path
        Directory from which relative paths are computed.

    rest of parameters: parameters given to getBarSize

    Returns
    -------
    BarResult
        Result object with status, computed bar size, and message.
    """
    relative_path = make_relative_path(file_path, base_dir)

    if not file_path.exists():
        return BarResult(
            relative_path=relative_path,
            input_galfit=None,
            bar_size=None,
            bar_size_arcsec=None,
            success=False,
            message="File does not exist.",
        )

    original_cwd = Path.cwd()
    try:
        os.chdir(file_path.parent)
        plate = Galfit(file_path.name).ReadHead().scale
        inputimage = Galfit(file_path.name).ReadHead().inputimage

        bar, N, theta = getBarSize(file_path.name, dis, numcomp, plot, ranx, out, red)
        bar_rounded = round(float(bar), BAR_SIZE_DECIMALS)
        bar_arcsec_rounded = round(float(bar * plate), BAR_SIZE_DECIMALS)
        return BarResult(
            relative_path=relative_path,
            input_galfit=inputimage,
            bar_size=bar_rounded,
            bar_size_arcsec=bar_arcsec_rounded,
            success=True,
            message="OK",
        )

    except Exception as exc:
        return BarResult(
            relative_path=relative_path,
            input_galfit=None,
            bar_size=None,
            bar_size_arcsec=None,
            success=False,
            message=f"{type(exc).__name__}: {exc}",
        )

    finally:
        os.chdir(original_cwd)


def process_files(
    file_paths: Iterable[Path], base_dir: Path, dis, numcomp, plot, ranx, out, red
) -> list[BarResult]:
    """
    Process all files in the input iterable.

    Parameters
    ----------
    file_paths : Iterable[Path]
        Files to process.
    rest of parameters: parameters given to GetBarSize

    Returns
    -------
    list[BarResult]
        Collected results for all files.
    """
    results: list[BarResult] = []

    for file_path in file_paths:
        result = process_file(file_path, base_dir, dis, numcomp, plot, ranx, out, red)
        results.append(result)

    return results


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


def write_results(results: Iterable[BarResult], output_file: Path) -> None:
    """
    Write processing results to a CSV file.

    Parameters
    ----------
    results : Iterable[BarResult]
        Results to write.
    output_file : Path
        Destination CSV file.
    """

    with output_file.open("w", newline="", encoding=DEFAULT_ENCODING) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "relative_path",
                "input_galfit",
                "bar_size",
                "bar_size_arcsec",
                "success",
                "message",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    str(result.relative_path),
                    str(result.input_galfit),
                    result.bar_size,
                    result.bar_size_arcsec,
                    result.success,
                    result.message,
                ]
            )


def print_summary(results: Iterable[BarResult]) -> None:
    """
    Print a short summary of the processing results.

    Parameters
    ----------
    results : Iterable[BarResult]
        Results to summarize.
    """
    results_list = list(results)
    total = len(results_list)
    ok = sum(result.success for result in results_list)
    failed = total - ok

    print(f"Processed files: {total}")
    print(f"Successful: {ok}")
    print(f"Failed: {failed}")


# ============================================================================
# Main entry point
# ============================================================================
def getMulBarSize(
    list_file: str,
    dis: int,
    numcomp: int,
    plot: bool,
    ranx: list,
    out: str,
    red: bool,
    output_file: str,
) -> float:
    """gets the bar length of the spiral galaxies


    It runs getBarSize to estimate the barlength over
    multiple files given the path of each one in a file.

    It assumes the bar model is the second component of the
    GALFIT file. Bar model can be a Sersic or Ferrer function.
    The rest of components must be Sersic (or related) models.

    Parameters
    ----------
    file: str
          file containing the paths of every GALFIT file

    dis : int
         maximum list among components
    numcomp : int
              Number of component where it'll obtain center
              of all components. in other words it selects
              the galaxy that contains the bar if simultaneous
              fitting of galaxies was used.
    plot : bool
            If True, it draws plots of the break and kappa radius
    ranx : list
        range of search (xmin to xmax) for the kappa radius and break
        radius. If None, it will search in a range of r=1 to 2.5*Re
        of effetive radius of the bar model.
    out : str
         Name of the output file for the DS9 ellipse region marking
         the bar.
    output_file : str
         Name of the output file for the output barlength results

    red : bool
            If True, draws DS9 region ellipse as red color

    Returns
    -------
    int
        Exit status code.

    See also
    --------
    getBarSize: get the bar length

    """

    base_dir = Path.cwd().resolve()
    list_file = Path(list_file).expanduser().resolve()
    output_file = Path(output_file).expanduser().resolve()

    try:

        file_paths = read_file_list(list_file)
        results = process_files(
            file_paths, base_dir, dis, numcomp, plot, ranx, out, red
        )
        write_results(results, output_file)
        print_summary(results)
        print(f"Results written to: {output_file.name}")
        return 0

    except Exception as exc:
        print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


def mainbatchGetBarSize(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBarSize: gets the bar size from Sersic and Ferrer models"
    )
    p.add_argument("InputFile", help="file containing a list of file path GALFIT files")
    p.add_argument(
        "-d", "--dis", type=int, default=3, help="maximum distance among components"
    )
    p.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=1,
        help="number of component to be selected",
    )
    p.add_argument(
        "-o",
        "--out",
        type=str,
        default="barlength.reg",
        help="output DS9 ellipse region",
    )
    p.add_argument(
        "-co", "--output", type=str, default="barlenghts.csv", help="output csv file"
    )
    p.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="plots a file of kappa and break radius",
    )
    p.add_argument(
        "-r",
        "--red",
        action="store_true",
        help="If activated, DS9 region ellipse is red",
    )
    p.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="range of radius to search for barlength",
    )
    a = p.parse_args(argv)

    ret = getMulBarSize(
        a.InputFile, a.dis, a.numcomp, a.plot, a.ranx, a.out, a.red, a.output
    )
    print("done ")

    return 0
