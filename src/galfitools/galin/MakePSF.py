#!/usr/bin/env python3


import re
import subprocess as sp
import sys
import shutil

from galfitools.galin.galfit import (
    Galfit,
    GalHead,
    GalSky,
    galfitLastFit,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)
from galfitools.galin.getStar import getStar
from galfitools.galin.std import Ds9ell2Kronellv2
from galfitools.galin.std import GetInfoEllip
from galfitools.mge.mge2galfit import mge2gal
from galfitools.galout.getRads import getReComp

from galfitools.galin.MaskDs9 import maskDs9


def makePSF(
    image: str,
    regfile: str,
    center: bool,
    psfout: str,
    sigma: str,
    imsize: int,
    sky: float,
    twist: bool,
    numgauss: int,
    fracrad=0.95,
) -> None:
    """makes PSF model

    It creates a PSF model image model for GALFIT.
    It model a star from an image using the MGE method

    Parameters
    ----------
    image : str
            name of the FITS image
    regfile: str
            DS9 ellipse region file containing the star.
    center: bool
            If True, the geometrical center of the DS9 ellipse
            region will be used. Otherwise, the center with the
            maximum count level of the DS9 ellipse region will be used.
    psfout: str
            name of the PSF output file
    sigma: str
            If available, the sigma image, otherwise set to None
    imsize: int
            size of the PSF image
    sky: float
            value of sky background
    twist: bool
            If true use the twist option for the MGE method
    numgauss: int
            maximum number of gaussians to be used.


    Returns
    -------
    None

    """

    #######################
    # reading options from galfit header

    galfitFile = "star.init"

    inputimage = "star.fits"

    sigfile = "psfsig.fits"

    # creation of GalfitFile
    galhead = GalHead()

    galsky = GalSky()

    galhead.inputimage = inputimage
    galhead.outimage = "star.out.fits"
    galhead.sigimage = sigfile
    galhead.constraints = "constar.txt"
    galhead.xmin = 1
    galhead.xmax = imsize
    galhead.ymin = 1
    galhead.ymax = imsize
    galhead.convx = 101
    galhead.convy = 101
    galhead.mgzpt = 22.5
    galhead.scale = 0.262
    galhead.scaley = 0.262

    galsky.sky = sky

    fout = open(galfitFile, "w")

    galPrintHeader(fout, galhead)

    galPrintSky(fout, 1, galsky)

    fout.close()

    #################

    # inputimage of galfit header will be the output of getStar
    getStar(image, regfile, imsize, center, sky, inputimage, sigma, sigfile)

    # checking the center position
    even = False

    if imsize % 2 == 0:  # pragma: no cover
        even = True

    if even:
        lx = imsize / 2

        xpeak = int(lx + 1)
        ypeak = int(lx + 1)

    else:  # pragma: no cover
        lx = imsize / 2

        xpeak = int(lx + 0.5)
        ypeak = int(lx + 0.5)

    ####

    psf = 0
    gauss = None
    freeser = False
    freesky = False

    xypos = [xpeak, ypeak]

    # calling to these two functions to obtain eps and theta
    obj_trash, xpos_trash, ypos_trash, rx_trash, ry_trash, angle_trash = GetInfoEllip(
        regfile
    )
    xx_trash, yy_trash, Rkron_trash, theta, eps = Ds9ell2Kronellv2(
        xpos_trash, ypos_trash, rx_trash, ry_trash, angle_trash
    )

    ellip = eps
    posang = theta
    regmgefile = None
    outname = mge2gal(
        galfitFile,
        regmgefile,
        center,
        psf,
        twist,
        gauss,
        freeser,
        freesky,
        numgauss,
        xypos=xypos,
        ellip=ellip,
        posang=posang,
    )

    print("calling GALFIT to fit MGE")

    rungal = "galfit  {}".format("mseGALFIT.txt")
    errgal = sp.run(
        [rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    )

    # check if galfit failed
    if re.search("Doh!", errgal.stdout):  # pragma: no cover
        print("ERROR: GALFIT has been unable to find a solution")
        print("Try to reduce the number of gaussians with -ng option")
        print("Exiting now")
        sys.exit(1)

    # check if galfit quit
    if re.search("QUIT", errgal.stdout):  # pragma: no cover
        print("ERROR: GALFIT has unexpectedly quit.")
        print(
            "Probably the cause is trying to constraint a component which doesn't exist "
        )
        print("Exiting now")
        sys.exit(1)

    lastfit_file = galfitLastFit(".")

    # estimating radius at % of light determined fracrad

    EffRad, totmag, meanme, me, N, theta = getReComp(lastfit_file, 3, fracrad, None, 1)

    EffRadb, totmag, meanme, me, N, theta = getReComp(
        lastfit_file, 3, fracrad, theta + 90, 1
    )

    if EffRadb < EffRad:
        axrat = EffRadb / EffRad
    else:
        axrat = EffRad / EffRadb

    makeDS9ellip("psfell.reg", xpeak, ypeak, EffRad, axrat, theta + 90)
    ###

    # making star PSF
    psfout2 = "psf_star.fits"
    shutil.copyfile(inputimage, psfout2)
    maskDs9(
        psfout2,
        "psfell.reg",
        0,
        None,
        False,
        0,
        skymean=None,
        skystd=None,
        invert=True,
    )

    #########
    galfit = Galfit(lastfit_file)
    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    # printing output file to create psf model
    fout1 = open("psfmodel.txt", "w")

    head.outimage = psfout
    head.P = 1

    galPrintHeader(fout1, head)

    index = 0

    for index, item in enumerate(galcomps.N):

        galPrintComp(fout1, index + 1, index, galcomps)

    galsky.skip = 1

    galPrintSky(fout1, index + 1, galsky)

    fout1.close()

    print("calling GALFIT again to create final model")

    rungal = "galfit  {}".format("psfmodel.txt")
    errgal = sp.run(
        [rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    )

    print("Done. PSF model file: ", psfout)
    print("Done. PSF Star model file: ", psfout2)

    print("Check the image, model and residual: ", outname)


def makeDS9ellip(out, X, Y, rad, AxRat, theta):
    """makes DS9 region

    It creates an ellipse DS9 region

    Parameters
    ----------
    out: str
            name of the DS9 output file
    X, Y: float, float
           X, Y position of the center of PSF
    rad: float
         radius of the axis major of ellipse
    AxRat: float
           Axis ratio
    theta: float
           angular position of the ellipse

    Returns
    -------
    0

    """

    # now it creates the ellipse region file
    fout = open(out, "w")

    color = "red"

    line = "# Region file format: DS9 version 4.1 \n"
    fout.write(line)
    linea = "global color=" + color + " dashlist=8 3 width=2 "
    lineb = 'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
    linec = "fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"
    line = linea + lineb + linec
    fout.write(line)
    line = "physical\n"
    fout.write(line)

    radminor = rad * AxRat

    elline = "ellipse({:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}) \n".format(
        X, Y, radminor, rad, theta
    )
    fout.write(elline)

    fout.close()

    return 0


#############################################################################
#  End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
