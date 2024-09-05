#!/usr/bin/env python3


import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import GetInfoEllip
from galfitools.galin.std import Ds9ell2Kronellv2


def getPeak(image: str, regfile: str, center: bool, maskfile: str) -> None:
    """Returns the coordinates of the pixel with the
        highest value of counts from a DS9 ellipse region file


    Parameters
    ----------
    image : str
            name of the image FITS file
    regfile : str
            name of the DS9 region file that containts the ellipse
    center : bool
            if True it uses the geometrical center of the ellipse instead
            of the peak
    maskfile : str
              name of the mask FITS file

    Returns
    -------
    float, float
                X,Y coordinates of the pixel with the highest count value
    axratio : float
             axis ratio of the DS9 ellipse region file

    theta : float



    """

    # root_ext = os.path.splitext(image)
    # timg= root_ext[0]

    if maskfile:

        errmsg = "file {} does not exist".format(maskfile)
        assert os.path.isfile(maskfile), errmsg
        hdu = fits.open(maskfile)
        mask = hdu[0].data
        hdu.close()

    else:
        mask = np.array([])

    # mask = np.array([]) #empty for the moment

    i = 0

    hdu = fits.open(image)
    img = hdu[i].data

    hdu.close()

    (ncol, nrow) = GetAxis(image)

    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)
    xx, yy, Rkron, theta, eps = Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle)

    if center:
        print("center of DS9 ellipse region will be used")
        xpeak, ypeak = round(xpos), round(ypos)
    else:
        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, eps, ncol, nrow)

        xpeak, ypeak = GetPmax(img, mask, xmin, xmax, ymin, ymax)

    # print("object found at ", xpeak + 1, ypeak + 1)

    axratio = 1 - eps

    return xpeak + 1, ypeak + 1, axratio, theta


####################################################
####################################################
#####################################################


def GetPmax(image, mask, xmin, xmax, ymin, ymax):
    """get the pixel coordinates with the highest value

    # repeated
    """

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    chuckimg = image[ymin - 1 : ymax, xmin - 1 : xmax]
    if mask.any():
        chuckmsk = mask[ymin - 1 : ymax, xmin - 1 : xmax]

        invmask = np.logical_not(chuckmsk)

        invmask = invmask * 1

        chuckimg = chuckimg * invmask

    maxy, maxx = np.where(chuckimg == np.max(chuckimg))

    xpos = maxx[0] + xmin - 1
    ypos = maxy[0] + ymin - 1

    return (xpos, ypos)


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
