#! /usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import csv
import os
import sys


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


@dataclass(slots=True)
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


def getBarSize(
    galfitFile: str,
    dis: int,
    num_comp: int,
    plot: bool,
    ranx: list,
    out: str,
    red: bool,
) -> float:
    """gets the bar size of the spiral galaxies

    It takes the average of Kappa radius (maximum curvature) and
    Break radius (maximum of double derivative) to estimate the
    bar size of the three composed model of bulge, bar, and disk.

    It assumes the bar model is the second component of the
    GALFIT file. Bar model can be a Sersic or Ferrer function.
    The rest of components must be Sersic (or related) models.

    Parameters
    ----------
    galfitFile : str
    dis : int
    num_comp : int
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
    red : bool
            If True, draws DS9 region ellipse as red color

    Returns
    -------
    rbar : float
           bar size in pixels
    N : int
        number of components of the galaxy
    theta : float
        angular position of the galactic's bar

    See also
    --------
    getBreak2 : get the break radius
    getKappa2 : get the kappa radius


    """

    galfit = Galfit(galfitFile)

    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps = SelectGal(comps, dis, num_comp)

    maskgal = comps.Active == 1
    # it assumes bar is positioned as second galfit component
    # i.e. [1]
    theta = comps.PosAng[maskgal][1]

    AxRat = comps.AxRat[maskgal][1]
    X = comps.PosX[maskgal][1]
    Y = comps.PosY[maskgal][1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute bar size")
        print("exiting..")
        sys.exit(1)

    #########################
    # computing the slope
    #########################
    if ranx:  # pragma: no cover
        (xmin, xmax) = ranx[0], ranx[1]
    else:
        Re = comps.Rad[maskgal][
            1
        ]  # it assumes bar is positioned as second galfit component

        # it assumes bar size is in this range. Hopefully it founds the solution there:
        xmin = 1
        xmax = 2.5 * Re

        ranx = [xmin, xmax]

    rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)

    rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)

    # bar size is just the average of these two radius

    rbar = (rbreak + rkappa) / 2

    # now it creates the ellipse region file
    fout = open(out, "w")

    if red:
        color = "red"
    else:
        color = "blue"

    line = "# Region file format: DS9 version 4.1 \n"
    fout.write(line)
    linea = "global color=" + color + " dashlist=8 3 width=2 "
    lineb = 'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
    linec = "fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"
    line = linea + lineb + linec
    fout.write(line)
    line = "physical\n"
    fout.write(line)

    rbarminor = rbar * AxRat

    elline = "ellipse({:.2f}, {:.2f}, {:.2f}, {:.2f} {:.2f}) \n".format(
        X, Y, rbarminor, rbar, theta
    )
    fout.write(elline)

    fout.close()

    return rbar, N, theta


"""
Template to process a list of file paths.

This script reads a text file containing one file path per line, changes to the
directory of each file, applies `get_bar_size()` to the file, stores the results,
and optionally writes the collected results to an output file.

Python 3.11+
"""


def get_bar_size(file_path: Path) -> float:
    """
    Compute the bar size for a given file.

    Replace the body of this function with the actual implementation.

    Parameters
    ----------
    file_path : Path
        Path to the file to be processed.

    Returns
    -------
    float
        Computed bar size.
    """
    # Replace this placeholder with the real implementation.
    # Example:
    # return GetBarSize(str(file_path))
    raise NotImplementedError("Implement `get_bar_size()` with your real logic.")


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
