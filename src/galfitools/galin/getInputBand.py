#! /usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np
import argparse
import scipy
import os
from galfitools.galin.galfit import (
    GalComps,
    Galfit,
    GalHead,
    GetRadAng,
    SelectGal,
    conver2Sersic,
    conver2Ferrer,
    numComps,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)

from galfitools.galout.PhotDs9 import photDs9


def getInputBand(
    galfitFile: str,
    ds9region: str,
    numcomp: int,
    zeropoint=22.5,
    fileout="newfit.init",
):
    """computes initial conditions from another fitting.


    computes initial parameters conditions from another
    galfit file and a DS9 ellipse region

    If you have already a fit from another filter/band
    this can be used to compute the initial conditions
    for a different wavelenght. This programs
    uses the magnitude computed from the region enclosed
    in the DS9 ellipse to redistribute the magnitude
    among components  based on the previous magnitude
    distribution from original galfit file.


    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    ds9region: str
             DS9 ellipse region to compute the magnitude

    numcomp: int
           number of component from GALFIT file.

    zeropoint: float
             zeropoint of the new band


    fileout: str
            name of the GALFIT output file

    Returns
    -------
    galcomps: GalComps
              data class defined in galfitools.galin.galfit


    See Also
    -------
    bestInputParam: selects best input from a list of galfit file


    Notes
    -----

    """
    dis = 3  # maximum distance among components

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    galcomps = SelectGal(galcomps, dis, numcomp)

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute bulge_total ratio")
        print("exiting..")
        sys.exit(1)

    maskgal = comps.Active == 1

    # assuming the parameter here is magnitude. Does not apply for nuker, ferrer
    # and other SB models which uses mag/''^2 instead of mag

    comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

    totFlux = comps.Flux[maskgal].sum()

    comps.PerLight[maskgal] = comps.Flux[maskgal] / totFlux

    # calling photds9 to redistribute the new magnitude

    sky = galsky.sky
    RegFile = ds9region
    ImageFile = head.inputimage
    # zeropoint = head.mgzpt
    mask = head.maskimage
    plate = head.scale
    magds9, sb, exptime = photDs9(ImageFile, RegFile, mask, zeropoint, plate, sky)

    fluxds9 = 10 ** ((zeropoint - magds9) / 2.5)  # it uses same magzpt

    comps.Flux[maskgal] = fluxds9 * comps.PerLight[maskgal]

    comps.Mag[maskgal] = -2.5 * np.log10(comps.Flux[maskgal]) + zeropoint

    head.mgzpt = zeropoint  # new zeropoint

    fout = open(fileout, "w")

    galPrintHeader(fout, head)
    index = 0
    for index, item in enumerate(comps.N):

        galPrintComp(fout, index + 1, index, comps)

    galPrintSky(fout, index + 1, galsky)

    fout.close()

    return comps


def get_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Compute initial GALFIT conditions from another band "
            "using a GALFIT file and a DS9 ellipse region."
        )
    )

    parser.add_argument(
        "galfitFile",
        type=str,
        help="Input GALFIT file.",
    )

    parser.add_argument(
        "ds9region",
        type=str,
        help="DS9 region file used to compute the magnitude.",
    )

    parser.add_argument(
        "numcomp",
        type=int,
        help="Number of the GALFIT component.",
    )

    parser.add_argument(
        "-zp",
        "--zeropoint",
        type=float,
        default=22.5,
        help="Zeropoint of the new band. Default: 22.5.",
    )

    parser.add_argument(
        "-o",
        "--fileout",
        type=str,
        default="newfit.init",
        help="Output GALFIT file. Default: newfit.init.",
    )

    return parser


def read_args(args=None) -> argparse.Namespace:
    """Read command-line arguments."""

    parser = get_parser()
    return parser.parse_args(args)


def mainGetInputBand(args=None):
    """Command-line interface for getInputBand."""

    args = read_args(args)

    galcomps = getInputBand(
        args.galfitFile,
        args.ds9region,
        args.numcomp,
        zeropoint=args.zeropoint,
        fileout=args.fileout,
    )

    return galcomps


if __name__ == "__main__":
    mainGetInputBand()
