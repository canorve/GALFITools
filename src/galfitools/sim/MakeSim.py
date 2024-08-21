#! /usr/bin/env python3

import numpy as np
from astropy.io import fits
from galfitools.galin.MaskDs9 import GetAxis


def makeSim(image, GAIN, skymean, skystd, newimage) -> None:
    """Simulates a observed galaxy from a GALFIT model

    Makes a simple artificial galaxy model. It adds Poisson
    noise and sky noise to the GALFIT model.


    Parameters
    ----------
    image : str
            name of the FITS image that contains the galaxy model
    GAIN : float
           CCD's gain e-/ADU
    skymean : float
            value of the mean of sky background
    skystd : float
            value of standard deviation of sky background
    newimage : str
             name of the new simulated image

    Returns
    -------
    None

    """
    sizex, sizey = GetAxis(image)

    hdu = fits.open(image)

    img = hdu[0].data

    eimg = GAIN * img

    noisyimg = np.random.poisson(eimg)

    pimg = noisyimg / GAIN

    hdu[0].data = pimg

    hdu.writeto("poissonimg.fits", overwrite=True)

    sky = np.random.normal(skymean, skystd, (sizex, sizey))

    hdusky = fits.PrimaryHDU()
    hdusky.data = sky

    hdusky.writeto("skynoise.fits", overwrite=True)

    # adds the two images:

    imgsim = pimg + sky

    hdu[0].data = imgsim

    hdu.writeto(newimage, overwrite=True)

    hdu.close()
