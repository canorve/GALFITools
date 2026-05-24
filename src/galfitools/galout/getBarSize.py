#! /usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional
import csv
import os
import sys


from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import getBreak2, getKappa2, getSlope
from galfitools.galout.getRads import getReComp

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
    scale=1.0,
    method="break_kappa",
) -> tuple[float, int, float]:
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

    scale: float
            constant to multiply the bar length. Default =1

    method: str
            indicates which method is used to measure the
            bar length. Options include 'break_kappa', 'break'
            'kappa','re', 'all'. Default='break_kappa'

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

    rbar = 0

    options = ["break_kappa", "break", "kappa", "re", "all"]
    if not (method in options):
        print("option not found. Setting to break_kappa")
        print("options available: break_kappa, break, kappa, re, all")
        method = "break_kappa"

    print(f"method used: {method}")

    if method == "break_kappa":
        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)

        rbar = scale * ((rbreak + rkappa) / 2)

    if method == "break":

        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar = scale * rbreak

    if method == "kappa":

        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar = scale * rkappa

    if method == "re":
        rbar = scale * comps.Rad[maskgal][1]

        # slope = 3.12
        # rgam, N, theta = getSlope(galfitFile, dis, slope, theta, num_comp, plot, ranx)
        # rbar = rgam

    if method == "all":

        # angle = comps.PosAng[maskgal][1]
        # fracrad = 0.5
        # EffRad, totmag, meanme, me, N, theta = getReComp(
        #    galfitFile, 3, fracrad, angle, num_comp
        # )

        # rbar0 = EffRad

        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar1 = rkappa

        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar2 = rbreak

        rbar3 = comps.Rad[maskgal][1]

        # rbar = (rbar0+rbar1 +rbar2+rbar3)/4
        rbar = scale * (rbar1 + rbar2 + rbar3) / 3

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
