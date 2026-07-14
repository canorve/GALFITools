#!/usr/bin/env python3

import numpy as np
import os
from galfitools.galin.galfit import GalComps
from galfitools.galin.galfit import GalHead
from galfitools.galin.galfit import GalSky
from galfitools.galin.galfit import galPrintHeader
from galfitools.galin.galfit import galPrintComp
from galfitools.galin.galfit import galPrintSky

from galfitools.galout.getPeak import getPeak
from galfitools.galout.PhotDs9 import photDs9

from galfitools.galin.std import GetInfoEllip
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import Ds9ell2Kronell


def getSersic(
    image: str,
    regfile: str,
    center: bool,
    maskfile: str,
    zeropoint: float,
    sky: float,
    noprint: bool,
    bulgetot: float,
    bards9: str,
    plate: float,
    out: str,
    galfit_out=None,
    freeser=False,
    consbulge=False,
    nser=2,
    bulgebarat=1,
) -> GalComps:
    """Obtains the initial parameters for GALFIT

    given an DS9 ellipse region, it prints the initial parameters
    of the surface brightness model in GALFIT format.


    Parameters
    ----------
    image : float
            FITS image of the galaxy

    regfile : str
            DS9 ellipse region file

    center : bool
            If True, it uses the geometrical center of the
            DS9 ellipse. Otherwise it will use the peak pixel
            as a center of the galaxy.

    maskfile : str
              mask file

    zeropoint : float
               magnitude zero point

    sky : float
          Sky/background value

    noprint : bool
            if False, avoids to print to STDOUT


    nser : sersic index initial parameter. Default = 2

    bulgetot : float
            Estimates the bulge-to-total ratio of the galaxy to
            determine the initial parameters for the bulge/disk model.
            If set to None, it computes the initial parameters for
            a single Sersic model. Here it is assumed that bar
            is part of the bulge (see bulgebarat below)

    bards9 : str
            If a DS9 ellipse region is provided (different from
            the regfile defined above), it will be used to estimate
            the initial parameters of the bar. In this case, the
            initial parameters of the bulge, bar, and disk will
            be printed. The bulgetot parameter must be specified
            for this to function.

    bulgebarat : str
            If bards9 is activated, this indicated the
            bulge/bar flux ratio. It divides the magnitude
            that correspond to the bulge in bulgetot into
            bulge and bar. Default = 1


    plate: float
          plate scale

    out : str
          output GALFIT file. This is a file where the surface
          brightness model is formated as a GALFIT file without
          the header.

    galfit_out : str
          Name of the output GALFIT cube fits. This is the name
          of the cube image after the fit.

    freeser : bool
            keeps the Sersic index of the seccond component
            as free

    consbulge: bool
            add constraints to the bulge q > 0.6
            and for bar q < 0.6 if bards9 is activated



    Returns
    -------
    galcomps : GalComps data class defined in galfit.py containing
               the initial parameters

    """

    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)

    mag, sb, exptime = photDs9(image, regfile, maskfile, zeropoint, plate, sky)

    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)

    (ncol, nrow) = GetAxis(image)

    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos, ypos, rx, ry, angle)
    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)

    # enlarging size of the fitting region by 6

    xsize = 6 * (xmax - xmin)
    ysize = 6 * (ymax - ymin)

    xmax = round(xpos + xsize / 2)
    xmin = round(xpos - xsize / 2)

    ymax = round(ypos + ysize / 2)
    ymin = round(ypos - ysize / 2)

    # correcting for fitting regions outside of image

    if xmin < 1:
        xmin = 1

    if xmax > ncol:
        xmax = ncol

    if ymin < 1:
        ymin = 1

    if ymax > nrow:
        ymax = nrow

    if rx >= ry:
        Re = rx / 2  # wild guess
    else:
        Re = ry / 2  # same

    # computing the bulge_bar_total_ratio:
    # from bulge bar ratio
    barbulgerat = 1 / bulgebarat
    totbulgebar = 1 + barbulgerat
    bulgebartot = 1 / totbulgebar

    if bards9:
        Xbar, Ybar, AxRatbar, PAbar = getPeak(image, bards9, center, maskfile)
        magbar, sb, exptimebar = photDs9(image, bards9, maskfile, zeropoint, plate, sky)
        objbar, xposbar, yposbar, rxbar, rybar, anglebar = GetInfoEllip(bards9)

    if bulgetot:
        Fluxtot = 10 ** (-mag / 2.5)
        FluxBulge = Fluxtot * bulgetot
        FluxDisk = Fluxtot - FluxBulge
        mag = -2.5 * np.log10(FluxBulge)
        mag2 = -2.5 * np.log10(FluxDisk)

        if bards9:
            Fluxbar = FluxBulge * bulgebartot
            FluxBulge = FluxBulge * bulgebartot
            mag = -2.5 * np.log10(FluxBulge)
            magbar = -2.5 * np.log10(Fluxbar)

            if rxbar >= rybar:
                Rebar = rxbar / 2  # wild guess
            else:
                Rebar = rybar / 2  # same

            Rebulge = Rebar * bulgebartot
        else:
            Rebulge = Re * bulgetot

    skip = 0
    N = 1

    fileconst = "constraints.txt"

    # store in GalHead, GalComps and GalSky data class
    galcomps = GalComps()
    galhead = GalHead()
    galsky = GalSky()

    # sky setup
    galsky.sky = sky

    name, extension = os.path.splitext(image)

    # header setup
    galhead.inputimage = image

    if galfit_out is None:
        galhead.outimage = name + "-out.fits"
    else:
        galhead.outimage = galfit_out

    galhead.sigimage = "sigma.fits"
    galhead.psfimage = "psf.fits"
    galhead.maskimage = maskfile

    galhead.constraints = fileconst
    galhead.mgzpt = zeropoint

    galhead.convx = 100
    galhead.convy = 100

    galhead.xmin = xmin
    galhead.xmax = xmax
    galhead.ymin = ymin
    galhead.ymax = ymax
    galhead.scale = plate
    galhead.scaley = plate

    fserout = open(out, "w")

    galPrintHeader(fserout, galhead)

    idxcount = 0  # index for components

    if bulgetot:
        if not (noprint):
            print(
                "# The initial parameters for the Sersic component based on "
                + "the DS9 ellipse region are: "
            )
            print("")
            print(
                "# WARNING: these are initial parameters. True values will be "
                + "computed by GALFIT"
            )

            printTerminal(1, X, Y, mag, Rebulge, nser, 1, 0, skip, 1)

            if bards9:

                if freeser is True:

                    printTerminal(2, X, Y, magbar, Rebar, 0.5, AxRatbar, PAbar, skip, 1)

                else:

                    printTerminal(2, X, Y, magbar, Rebar, 0.5, AxRatbar, PAbar, skip, 0)

            printTerminal(3, X, Y, mag2, Re, 1, AxRat, PA, skip, 0)

        print("# parameter constraints file: ", fileconst)

        fout = open(fileconst, "w")

        if bards9:

            print("# 1_2_3    x    offset ")
            print("# 1_2_3    y    offset ")
            constlinex = " 1_2_3    x    offset \n"
            constliney = " 1_2_3    y    offset \n"
            fout.write(constlinex)
            fout.write(constliney)

            if consbulge:

                print("# 1    q  0.5 to 1  ")
                print("# 2    q  0 to 0.65 ")
                constlinebulge = " 1   q  0.5 to 1 \n"
                constlinebar = " 2   q  0 to 0.65 \n"
                fout.write(constlinebulge)
                fout.write(constlinebar)

        else:

            print("# 1_2    x    offset ")
            print("# 1_2    y    offset ")

            constlinex = " 1_2    x    offset \n"
            constliney = " 1_2    y    offset \n"
            fout.write(constlinex)
            fout.write(constliney)

            if consbulge:

                print("# 1    q  0.5 to 1  ")
                # print("# 1    n  0.1 to 10  ")
                constlinebulge = " 1   q  0.5 to 1 \n"
                # constlinesersic = " 1   n  0.1 to 10 \n"
                fout.write(constlinebulge)
                # fout.write(constlinesersic)

        fout.close()

        # first component: bulge

        galcomps = copy2Galcomps(
            galcomps, N, "sersic", X, Y, mag, Rebulge, nser, 1, 0, skip, 1
        )

        galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
        idxcount += 1

        if bards9:

            # alternative component: bar
            N = N + 1

            if freeser is True:
                galcomps = copy2Galcomps(
                    galcomps,
                    N,
                    "sersic",
                    X,
                    Y,
                    magbar,
                    Rebar,
                    0.5,
                    AxRatbar,
                    PAbar,
                    skip,
                    1,
                )
            else:
                galcomps = copy2Galcomps(
                    galcomps,
                    N,
                    "sersic",
                    X,
                    Y,
                    magbar,
                    Rebar,
                    0.5,
                    AxRatbar,
                    PAbar,
                    skip,
                    0,
                )

            galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
            idxcount += 1

        # second component: disk
        N = N + 1

        galcomps = copy2Galcomps(
            galcomps, N, "sersic", X, Y, mag2, Re, 1, AxRat, PA, skip, 0
        )

        galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
        idxcount += 1

    else:

        if not (noprint):
            print(
                "# The initial parameters for the Sersic component based on "
                + "the DS9 ellipse region are: "
            )
            print("")
            print(
                "# WARNING: these are initial parameters. True values will be"
                + " computed by GALFIT"
            )
            print("")

            printTerminal(1, X, Y, mag, Re, nser, AxRat, PA, skip, 1)

        print("# parameter constraints file: ", fileconst)

        fout = open(fileconst, "w")

        constline = "  \n"
        fout.write(constline)

        fout.close()

        galcomps = copy2Galcomps(
            galcomps, N, "sersic", X, Y, mag, Re, nser, AxRat, PA, skip, 1
        )

        galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
        idxcount += 1

    galPrintSky(fserout, idxcount + 1, galsky)

    fserout.close()

    return galcomps


def printTerminal(num, X, Y, mag, Re, nser, axis, pa, skip, freeser):
    """print component to terminal"""

    print("")

    print("# Component number: {}".format(num))
    print("0) sersic # Component type")
    print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
    print("3) {:.2f}   1       # Integrated magnitude ".format(mag))
    print("4) {:.2f}   1       # R_e (effective radius) ".format(Re))
    print("5) {:.2f}   {}       # Sersic index n  ".format(nser, freeser))
    print("6) 0.0000   0      # ----  ")
    print("7) 0.0000   0      # ----  ")
    print("8) 0.0000   0      # ----  ")
    print("9) {:.2f}   1       # Axis Ratio (b/a)  ".format(axis))
    print("10) {:.2f}  1       # Position angle (PA)  ".format(pa))
    print("Z) {}  # Skip this model in output image?  ".format(skip))

    print("")


def copy2Galcomps(galcomps, N, comp, X, Y, mag, Re, nser, axrat, pa, skip, freeser):
    """copy component parameter value to GalComps"""

    galcomps.N = np.append(galcomps.N, N)
    galcomps.NameComp = np.append(galcomps.NameComp, comp)
    galcomps.PosX = np.append(galcomps.PosX, X)
    galcomps.PosY = np.append(galcomps.PosY, Y)

    galcomps.Mag = np.append(galcomps.Mag, mag)
    galcomps.Rad = np.append(galcomps.Rad, Re)
    galcomps.Exp = np.append(galcomps.Exp, nser)
    galcomps.Exp2 = np.append(galcomps.Exp2, 0)
    galcomps.Exp3 = np.append(galcomps.Exp3, 0)
    galcomps.AxRat = np.append(galcomps.AxRat, axrat)
    galcomps.PosAng = np.append(galcomps.PosAng, pa)
    galcomps.skip = np.append(galcomps.skip, skip)

    # free parameters
    galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
    galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
    galcomps.MagFree = np.append(galcomps.MagFree, 1)
    galcomps.RadFree = np.append(galcomps.RadFree, 1)
    galcomps.ExpFree = np.append(galcomps.ExpFree, freeser)
    galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
    galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
    galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
    galcomps.PosAngFree = np.append(galcomps.PosAngFree, 1)

    return galcomps
