#!/usr/bin/env python3


import subprocess as sp
import numpy as np

from astropy.io import fits
from pathlib import Path

from galfitools.galin.MaskDs9 import maskDs9
from galfitools.galout.getRads import getReComp

from galfitools.galin.galfit import (
    Galfit,
    readDataImg,
    numParFree,
    SelectGal,
)


from galfitools.galin.MakePSF import makeDS9ellip


def getChiNu(galfile, numcomp, fracrad=0.98, ds9reg=None, delete=False):
    """
    getChiNu computes the Chinu square

    computes the Chinu square inside a ellipse
    with the fraction of total light given by fracrad.

    returns Chinu, Akaike information Criterion,
    Bayesian information criterion and total
    of free parameters



    """

    galfit = Galfit(galfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    galcomps = SelectGal(galcomps, 3, numcomp)

    filename = Path(galfile)

    root = filename.stem
    extension = filename.suffix.lstrip(".")

    newname = root + "-" + extension + "-sigma.fits"
    output_sigma = newname

    newname = root + "-" + extension + "-chisquare.fits"
    output_chisqr = newname

    maskgal = galcomps.Active == 1
    # axrat = galcomps.AxRat[maskgal][1]

    xpeak = galcomps.PosX[maskgal][1]
    ypeak = galcomps.PosY[maskgal][1]

    xmin = galhead.xmin
    ymin = galhead.ymin

    xmax = galhead.xmax
    ymax = galhead.ymax

    # sigma image
    sigma_im = galhead.sigimage

    # correcting for cube galfit image:
    xpeak = xpeak - xmin
    ypeak = ypeak - ymin

    dis = 3  # maximum distance among components

    EffRad, totmag, meanme, me, N, theta = getReComp(
        galfile, dis, fracrad, None, numcomp
    )

    EffRadb, totmag, meanme, me, N, theta = getReComp(
        galfile, dis, fracrad, theta + 90, numcomp
    )

    if EffRadb < EffRad:
        axrat = EffRadb / EffRad
    else:
        axrat = EffRad / EffRadb
        EffRad, EffRadb = EffRadb, EffRad

    chifile_path = Path("chinu.reg")

    if not chifile_path.exists():
        makeDS9ellip("chinu.reg", xpeak, ypeak, EffRad, axrat, theta + 90)

    if not ds9reg:
        makeDS9ellip("chinu.reg", xpeak, ypeak, EffRad, axrat, theta + 90)

    file_path = Path(sigma_im)

    if file_path.exists():
        print("Creating a sigma image from file ")

        # Open the FITS file
        with fits.open(sigma_im) as hdul:
            data = hdul[0].data
            header = hdul[0].header.copy()

            # Cut the image
            new_data = data[ymin - 1 : ymax, xmin - 1 : xmax]

            # Update header size keywords
            header["NAXIS1"] = new_data.shape[1]
            header["NAXIS2"] = new_data.shape[0]

            # Write the new FITS image
            fits.writeto(output_sigma, new_data, header, overwrite=True)

        print("calling GALFIT to create output image")
        rungal = f"galfit -o2  {galfile}"

        errgal = sp.run(
            rungal,
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            text=True,
        )

        dataimg = readDataImg(galhead, sigma=output_sigma)

    else:
        print("The sigma file does not exist")

        print("calling GALFIT to create sigma")

        rungal = f"galfit -o2 -outsig {galfile}"

        errgal = sp.run(
            rungal,
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            text=True,
        )

        old_file = Path("sigma.fits")
        new_file = Path(output_sigma)

        old_file.rename(new_file)
        dataimg = readDataImg(galhead, sigma=output_sigma)

    sigflux = dataimg.sigma
    galflux = dataimg.img
    modflux = dataimg.model

    resflux = (galflux - modflux) ** 2

    varflux = sigflux**2

    chisquareimg = resflux / varflux

    if dataimg.mask.any():

        immask = dataimg.mask
        maskm = immask == True
        chisquareimg[maskm] = 0

    hdu = fits.PrimaryHDU(chisquareimg)

    hdu.writeto(output_chisqr, overwrite=True)

    if ds9reg:
        maskDs9(
            output_chisqr,
            ds9reg,
            0,
            None,
            False,
            0,
            skymean=None,
            skystd=None,
            invert=True,
        )

    else:
        maskDs9(
            output_chisqr,
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
    with fits.open(output_chisqr) as hdul:
        data = hdul[0].data  # Extract the image data
        # Calculate the total sum
        chisquare = np.sum(data)
        maskchi = data != 0

    pixcountchi = np.size(galflux[maskchi])
    totfreepar = numParFree(galcomps)
    ndof = pixcountchi - totfreepar
    chinu = chisquare / ndof

    aic = chisquare + 2 * totfreepar  # Akaike information criterion
    bic = chisquare + totfreepar * np.log(pixcountchi)  # Bayesian information criterion

    if delete:

        file_path = Path(output_sigma)
        file_path.unlink()

        filechi_path = Path(output_chisqr)
        filechi_path.unlink()

    return chinu, aic, bic, totfreepar
