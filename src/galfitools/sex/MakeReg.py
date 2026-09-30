#! /usr/bin/env python

from pathlib import Path

import numpy as np
from galfitools.sex.MakeMask import CheckFlag


def makeReg(
    sexfile: str,
    scale: float,
    region_dir="kron_regions",
    region_id=None,
) -> None:
    """Creates Ds9 ellipse regfiles from a catalog of SExtractor

    It creates a mask file for GALFIT using information from a
    SExtractor catalog. It includes masking of saturated regions.

    Parameters
    ----------
    sexfile : str
            name of the Sextractor catalog

    scale: float,
            Scale factor by which the ellipse will be enlarged or diminished.

    region_dir : str or pathlib.Path or None
            Directory for ellipse_<catalog ID>.reg files. Created automatically.
            Default: "kron_regions" relative to the working directory.
            Set to None to disable region export. Existing matching files are
            overwritten; unrelated files are retained.
    region_id : int or None
            Export only this catalog ID, or all masked objects if None.
            Objects excluded by saturation checks are not exported.

    Returns
    -------
    None

    """

    MakeRegs(
        sexfile,
        scale,
        0,
        region_dir=region_dir,
        region_id=region_id,
    )  # offset set to 0 for now


def MakeRegs(
    catfile,
    scale,
    offset,
    region_dir="kron_regions",
    region_id=None,
):
    """Creates ellipse masks for every object of the SExtractor catalog

    Parameters
    ----------
    maskimage: str
            name of the mask file. This file should already exists
    catfile: str,
            SExtractor catalog
    scale: float,
            Scale factor by which the ellipse will be enlarged or diminished.
    offset: float
            constant to be added to the ellipse size
    regfile: str
            DS9 region file containing the saturated region.

    region_dir : str or pathlib.Path or None
            Directory for ellipse_<catalog ID>.reg files. Created automatically.
            Default: "kron_regions" relative to the working directory.
            Set to None to disable region export. Existing matching files are
            overwritten; unrelated files are retained.
    region_id : int or None
            Export only this catalog ID, or all masked objects if None.
            Objects excluded by saturation checks are not exported.

    Returns
    -------
    None

    """

    checkflag = 0
    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag

    (
        n,
        alpha,
        delta,
        xx,
        yy,
        mg,
        kr,
        fluxrad,
        ia,
        ai,
        e,
        theta,
        bkgd,
        idx,
        flg,
    ) = np.genfromtxt(catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    Rkron = scale * ai * kr + offset

    mask = Rkron < 1
    if mask.any():  # pragma: no cover
        Rkron[mask] = 1

    print("Creating Ds9 ellipse region for every object \n")

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = CheckFlag(flg[idx], flagsat, maxflag)

        if checkflag is False:

            idn = n[idx]
            R = Rkron[idx]
            x = xx[idx]
            y = yy[idx]
            q = 1 - e[idx]
            bim = q * R
            angle = theta[idx]

            if region_id is None or val == region_id:
                selected_dir = region_dir
            else:
                selected_dir = None

            if selected_dir is not None:

                output_dir = Path(region_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                region_path = output_dir / f"ellipse_{int(idn)}.reg"
                region_path.write_text(
                    "# Region file format: DS9 version 4.1\n"
                    "global color=green width=1\n"
                    "image\n"
                    f"ellipse({x + 1:.10g},{y + 1:.10g},{R:.10g},{bim:.10g},"
                    f"{angle:.10g}) # text={{{int(idn)}}}\n",
                    encoding="utf-8",
                )

    return 0
