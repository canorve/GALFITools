#!/usr/bin/env python3


import re
import subprocess as sp
import sys

from galfitools.galin.galfit import (
    Galfit,
    galfitLastFit,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)
from galfitools.galin.getStar import getStar
from galfitools.galin.std import Ds9ell2Kronellv2
from galfitools.galin.std import GetInfoEllip
from galfitools.mge.mge2galfit import mge2gal


def makePSF(
    galfitFile: str,
    image: str,
    regfile: str,
    center: bool,
    psfout: str,
    sigma: str,
    twist: bool,
    numgauss: int,
) -> None:
    """makes PSF model

    It creates a PSF model image model for GALFIT.
    It model a star from an image usign the MGE method

    Parameters
    ----------
    galfitFile : str
                GALFIT initial parameters file
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

    galfit = Galfit(galfitFile)
    head = galfit.ReadHead()
    galsky = galfit.ReadSky()

    inputimage = head.inputimage

    sky = galsky.sky

    sigfile = head.sigimage

    imsize = head.xmax - head.xmin + 1

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

    print("Check the image, model and residual: ", outname)


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
