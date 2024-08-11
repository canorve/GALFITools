#! /usr/bin/env python

import sys

from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import getBreak2, getKappa2


def getBarSize(
    galfitFile: str, dis: int, num_comp: int, plot: bool, ranx: list, out: str
) -> float:
    """gets the Kappa radius (maximum curvature) from a set of
    Sersics using another method"""

    galfit = Galfit(galfitFile)

    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps = SelectGal(comps, dis, num_comp)

    maskgal = comps.Active == 1

    theta = comps.PosAng[maskgal][
        1
    ]  # it assumes bar is positioned as second galfit component

    AxRat = comps.AxRat[maskgal][1]
    X = comps.PosX[maskgal][1]
    Y = comps.PosY[maskgal][1]

    N = numComps(comps, "all")

    if N == 0:
        print("not enough number of components to compute bar size")
        print("exiting..")
        sys.exit(1)

    #########################
    # computing the slope
    #########################
    if ranx:
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

    line = "# Region file format: DS9 version 4.1 \n"
    fout.write(line)
    linea = "global color=red dashlist=8 3 width=2 "
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
