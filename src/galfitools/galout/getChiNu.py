#!/usr/bin/env python3


import subprocess as sp
import numpy as np

from astropy.io import fits

from galfitools.galin.MaskDs9 import maskDs9
from galfitools.galout.getRads import getReComp

from galfitools.galin.galfit import (
    Galfit,
    readDataImg,
    numParFree,
    SelectGal,
)


from galfitools.galin.MakePSF import makeDS9ellip


def getChiNu(galfitfile, numcomp, fracrad):

    galfit = Galfit(galfitfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    galcomps = SelectGal(galcomps, 3, numcomp)

    maskgal = galcomps.Active == 1

    # axrat = galcomps.AxRat[maskgal][1]

    xpeak = galcomps.PosX[maskgal][1]
    ypeak = galcomps.PosY[maskgal][1]

    xmin = galhead.xmin
    ymin = galhead.ymin

    # correcting for cube galfit image:
    xpeak = xpeak - xmin
    ypeak = ypeak - ymin

    EffRad, totmag, meanme, me, N, theta = getReComp(galfile, 3, fracrad, None, numcomp)

    EffRadb, totmag, meanme, me, N, theta = getReComp(
        lastfit_file, 3, fracrad, theta + 90, 1
    )

    if EffRadb < EffRad:
        axrat = EffRadb / EffRad
    else:
        axrat = EffRad / EffRadb

    makeDS9ellip("chinu.reg", xpeak, ypeak, EffRad, axrat, theta + 90)

    print("calling GALFIT to create sigma")

    rungal = "galfit -outsig -o3  {}".format(galfitfile)
    errgal = sp.run(
        [rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    )

    dataimg = readDataImg(galhead, sigma="sigma.fits")

    sigflux = dataimg.sigma
    galflux = dataimg.img
    modflux = dataimg.model

    resflux = (galflux - modflux) ** 2

    varflux = sigflux**2

    chisquareimg = resflux / varflux

    if dataimg.mask:

        immask = dataimg.mask
        maskm = immask == True
        chisquareimg[maskm] = 0

    hdu = fits.PrimaryHDU(chisquareimg)

    hdu.writeto("chisquare.fits", overwrite=True)

    maskDs9(
        "chisquare.fits",
        "chinu.reg",
        0,
        None,
        False,
        0,
        skymean=None,
        skystd=None,
        invert=True,
    )

    # Open the FITS file
    with fits.open("chisquare.fits") as hdul:
        data = hdul[0].data  # Extract the image data
        # Calculate the total sum
        chisquare = np.sum(data)
        maskchi = data != 0

    pixcountchi = np.size(galflux[maskchi])
    totfreepar = numParFree(galcomps)
    ndof = pixcountchi - totfreepar
    chinu = chisquare / ndof

    return chinu
