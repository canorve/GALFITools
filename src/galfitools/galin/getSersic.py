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
    image: float,
    regfile: str,
    center: bool,
    maskfile: str,
    zeropoint: float,
    sky: float,
    noprint: float,
    bulgetot: float,
    bards9: str,
    out: str,
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
    noprint : float
            if False, avoids to print to STDOUT

    bulgetot : float
            Estimates the bulge-to-total ratio of the galaxy to
            determine the initial parameters for the bulge/disk model.
            If set to None, it computes the initial parameters for
            a single Sersic model.

    bards9 : str
            If a DS9 ellipse region is provided (different from
            the regfile defined above), it will be used to estimate
            the initial parameters of the bar. In this case, the
            initial parameters of the bulge, bar, and disk will
            be printed. The bulgetot parameter must be specified
            for this to function.

    out : str
          output GALFIT file. This is a file where the surface
          brightness model is formated as a GALFIT file without
          the header.

    Returns
    -------
    galcomps : GalComps data class defined in galfit.py containing
               the initial parameters

    """

    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)

    mag, exptime = photDs9(image, regfile, maskfile, zeropoint, sky)

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

    if bards9:

        Xbar, Ybar, AxRatbar, PAbar = getPeak(image, bards9, center, maskfile)
        magbar, exptimebar = photDs9(image, bards9, maskfile, zeropoint, sky)
        objbar, xposbar, yposbar, rxbar, rybar, anglebar = GetInfoEllip(bards9)

    if bulgetot:
        Fluxtot = 10 ** (-mag / 2.5)
        FluxBulge = Fluxtot * bulgetot
        FluxDisk = Fluxtot - FluxBulge
        mag = -2.5 * np.log10(FluxBulge)
        mag2 = -2.5 * np.log10(FluxDisk)

        if bards9:
            Fluxbar = 10 ** (-magbar / 2.5)  # assuming bar flux containts bulge flux
            FluxDisk = Fluxtot - Fluxbar
            FluxBulge = Fluxbar * 0.3  # wild guess
            Fluxbar = Fluxbar * 0.7  # wild guess
            mag = -2.5 * np.log10(FluxBulge)
            mag2 = -2.5 * np.log10(FluxDisk)
            magbar = -2.5 * np.log10(Fluxbar)

            if rxbar >= rybar:
                Rebar = rxbar / 2  # wild guess
            else:
                Rebar = rybar / 2  # same

    if rx >= ry:
        Re = rx / 2  # wild guess
    else:
        Re = ry / 2  # same

    # wild guesses for n and Re
    n = 2
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
    galhead.outimage = name + "-out.fits"
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
    galhead.scale = 0.262  # plate scale for DESI
    galhead.scaley = 0.262  # plate scale for DESI

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

            print("")

            print("# Component number: 1")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
            print("3) {:.2f}   1       # Integrated magnitude ".format(mag))
            print(
                "4) {:.2f}   1       # R_e (effective radius) ".format(Re * bulgetot)
            )  # wild guess
            print("5) {:.2f}   1       # Sersic index n  ".format(n))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}   1       # Axis Ratio (b/a)  ".format(1))
            print("10) {:.2f}  1       # Position angle (PA)  ".format(0))
            print("Z) {}  # Skip this model in output image?  ".format(skip))

            print("")

            if bards9:
                print("# Component number: 1b bar")
                print("0) sersic # Component type")
                print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
                print("3) {:.2f}    1      # Integrated magnitude ".format(magbar))
                print("4) {:.2f}    1      # R_e (effective radius) ".format(Rebar))
                print("5) {:.2f}    0      # Sersic index n  ".format(0.5))
                print("6) 0.0000   0      # ----  ")
                print("7) 0.0000   0      # ----  ")
                print("8) 0.0000   0      # ----  ")
                print("9) {:.2f}   1       # Axis Ratio (b/a)  ".format(AxRatbar))
                print("10) {:.2f}  1       # Position angle (PA)  ".format(PAbar))
                print("Z) {}  # Skip this model in output image?  ".format(skip))
                print("")

            print("# Component number: 2")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
            print("3) {:.2f}    1      # Integrated magnitude ".format(mag2))
            print("4) {:.2f}    1      # R_e (effective radius) ".format(Re))
            print("5) {:.2f}    0      # Sersic index n  ".format(1))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}    1      # Axis Ratio (b/a)  ".format(AxRat))
            print("10) {:.2f}   1      # Position angle (PA)  ".format(PA))
            print("Z) {}  # Skip this model in output image?  ".format(skip))

            print("")
            print("# parameter constraints file: ", fileconst)

            fout = open(fileconst, "w")

            if bards9:

                print("# 1_2_3    x    offset ")
                print("# 1_2_3    y    offset ")
                constlinex = " 1_2_3    x    offset \n"
                constliney = " 1_2_3    y    offset \n"
                fout.write(constlinex)
                fout.write(constliney)

            else:

                print("# 1_2    x    offset ")
                print("# 1_2    y    offset ")

                constlinex = " 1_2    x    offset \n"
                constliney = " 1_2    y    offset \n"
                fout.write(constlinex)
                fout.write(constliney)

            fout.close()

        # first component: bulge
        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)

        galcomps.Mag = np.append(galcomps.Mag, mag)
        galcomps.Rad = np.append(galcomps.Rad, Re * bulgetot)
        galcomps.Exp = np.append(galcomps.Exp, n)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, 1)
        galcomps.PosAng = np.append(galcomps.PosAng, 0)
        galcomps.skip = np.append(galcomps.skip, skip)

        # free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 1)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
        galcomps.PosAngFree = np.append(galcomps.PosAngFree, 1)
        galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
        idxcount += 1

        if bards9:

            # alternative component: bar
            N = N + 1
            galcomps.PosX = np.append(galcomps.PosX, X)
            galcomps.PosY = np.append(galcomps.PosY, Y)
            galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
            galcomps.N = np.append(galcomps.N, N)

            galcomps.Mag = np.append(galcomps.Mag, magbar)
            galcomps.Rad = np.append(galcomps.Rad, Rebar)
            galcomps.Exp = np.append(galcomps.Exp, 0.5)
            galcomps.Exp2 = np.append(galcomps.Exp2, 0)
            galcomps.Exp3 = np.append(galcomps.Exp3, 0)
            galcomps.AxRat = np.append(galcomps.AxRat, AxRatbar)
            galcomps.PosAng = np.append(galcomps.PosAng, PAbar)
            galcomps.skip = np.append(galcomps.skip, skip)

            # free parameters
            galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
            galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
            galcomps.MagFree = np.append(galcomps.MagFree, 1)
            galcomps.RadFree = np.append(galcomps.RadFree, 1)
            galcomps.ExpFree = np.append(galcomps.ExpFree, 0)
            galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
            galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
            galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
            galcomps.PosAngFree = np.append(galcomps.PosAngFree, 1)

            galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
            idxcount += 1

        # second component: disk
        N = N + 1
        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)

        galcomps.Mag = np.append(galcomps.Mag, mag2)
        galcomps.Rad = np.append(galcomps.Rad, Re)
        galcomps.Exp = np.append(galcomps.Exp, 1)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
        galcomps.PosAng = np.append(galcomps.PosAng, PA)
        galcomps.skip = np.append(galcomps.skip, skip)

        # free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 0)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
        galcomps.PosAngFree = np.append(galcomps.PosAngFree, 1)

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

            print("# Component number: 1")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1  # Position x, y".format(X, Y))
            print("3) {:.2f}    1      # Integrated magnitude ".format(mag))
            print("4) {:.2f}    1      # R_e (effective radius) ".format(Re))
            print("5) {:.2f}    1      # Sersic index n  ".format(n))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}    1      # Axis Ratio (b/a)  ".format(AxRat))
            print("10) {:.2f}   1      # Position angle (PA)  ".format(PA))
            print("Z) {}  # Skip this model in output image?  ".format(skip))

        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)

        galcomps.Mag = np.append(galcomps.Mag, mag)
        galcomps.Rad = np.append(galcomps.Rad, Re)
        galcomps.Exp = np.append(galcomps.Exp, n)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
        galcomps.PosAng = np.append(galcomps.PosAng, PA)
        galcomps.skip = np.append(galcomps.skip, skip)

        # free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 1)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
        galcomps.PosAngFree = np.append(galcomps.PosAngFree, 1)

        galPrintComp(fserout, idxcount + 1, idxcount, galcomps)
        idxcount += 1

    galPrintSky(fserout, idxcount + 1, galsky)

    fserout.close()

    return galcomps
