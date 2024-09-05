#!/usr/bin/env python3


import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import GetInfoEllip
from galfitools.galin.std import Ds9ell2Kronellv2
from galfitools.galin.std import GetPmax


def getStar(
    image: str,
    regfile: str,
    imsize: int,
    center: bool,
    sky: float,
    imout: str,
    sigma: str,
    sigout: str,
) -> None:
    """extracts a piece of the image file

    Given a DS9 ellipse region file and a specified size,
    it extracts a cutout image that includes the region
    enclosed by the ellipse.

    Parameters
    ----------
    image : str
           name of the image file

    regfile : str
            DS9 ellipse region file
    imsize : int
            size of the new cutout image in pixels
    center : bool
            If True, it takes the geometrical center
            of the DS9 ellipse region. Otherwise, it will
            take the peak's pixels as the coordinate
    sky: float
         value of the background sky
    imout: str
        name of the output image file
    sigma: str
        name of the sigma image file if it exists otherwise None
    sigout: str
        name of the output sigma image file

    Returns
    -------
    None

    """

    ##########################################
    #  Internal flags
    even = False
    ##########################################
    ##########################################

    # root_ext = os.path.splitext(image)
    # timg= root_ext[0]

    mask = np.array([])  # empty for the moment

    i = 0  # index where data is located at hdu

    # if checkCompHDU(image):
    #    i = 1

    hdu = fits.open(image)
    img = hdu[i].data

    # img = img.astype(float)

    hdu.close()

    (ncol, nrow) = GetAxis(image)

    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)
    xx, yy, Rkron, theta, eps = Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle)

    if center:
        print("center of ds9 ellipse region will be used")
        xpeak, ypeak = round(xpos), round(ypos)
    else:
        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta + 90, eps, ncol, nrow)
        xpeak, ypeak = GetPmax(img, mask, xmin, xmax, ymin, ymax)

    print("object found at ", xpeak + 1, ypeak + 1)
    print("Ellipticity, Angle = ", eps, theta)

    if sky:
        print("Sky = ", sky)
    else:
        print("No sky input")

    if imsize % 2 == 0:
        even = True

    if even:
        lx = int(imsize / 2)

        xlo = xpeak - lx + 1
        ylo = ypeak - lx + 1

        xhi = xpeak + lx
        yhi = ypeak + lx

    else:
        lx = int(imsize / 2 + 0.5)

        xlo = xpeak - lx + 2
        ylo = ypeak - lx + 2

        xhi = xpeak + lx
        yhi = ypeak + lx

    GetFits(image, imout, sky, xlo, xhi, ylo, yhi)

    if sigma:
        GetFits(sigma, sigout, 0, xlo, xhi, ylo, yhi)


####################################################
####################################################
#####################################################


def GetFits(
    Image: str, Imageout: str, sky: float, xlo: int, xhi: int, ylo: int, yhi: int
) -> None:
    """Get a piece from the image

    Given an image, and X-minimum, X-maximum,
    Y-minimum, Y-maximum, it extracts a cutout image.
    the X-Y-minimum, X-Y-maximum determines the
    coordinates of the extracted neww image.

    Parameters
    ----------
    Image : str
           name of the image file
    Imageout : str
        name of the output image file

    sky : float
         value of the background sky
    xlo : int
          Pixel with the minimum X value.
    xhi : int
          Pixel with the maximum X value.
    ylo : int
          Pixel with the minimum Y value.
    yhi : int
          Pixel with the maximum Y value.


    Returns
    -------
    None

    # repeated

    """

    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
    # new fits

    i = 0  # index indicated where the data is located
    # if checkCompHDU(Image):
    #    i = 1
    #    newhdu = fits.CompImageHDU()
    # else:
    newhdu = fits.PrimaryHDU()

    hdu = fits.open(Image)
    dat = hdu[i].data[ylo - 1 : yhi, xlo - 1 : xhi]
    head = hdu[i].header

    if sky:
        newhdu.data = dat - sky
    else:
        newhdu.data = dat

    newhdu.header = head

    newhdu.writeto(Imageout, overwrite=True)

    hdu.close()


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
