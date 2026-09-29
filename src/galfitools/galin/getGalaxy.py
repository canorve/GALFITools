#!/usr/bin/env python3


import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import GetInfoEllip
from galfitools.galin.std import Ds9ell2Kronellv2

from galfitools.galin.getStar import GetFits


def getGalaxy(
    image: str,
    regfile: str,
    sky: float,
    imout: str,
    sigma: str,
    sigout: str,
    scale=6,
):
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
    sky: float
         value of the background sky
    imout: str
        name of the output image file
    sigma: str
        name of the sigma image file if it exists otherwise None
    sigout: str
        name of the output sigma image file
    scale: float
         enlarge the image for this factor. Default = 6

    Returns
    -------
    x_pos, y_pos: new pixel coordinate of the galaxy.

    x_cor, y_cor: x, y coordinates tranformation

    """

    (ncol, nrow) = GetAxis(image)

    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)

    # enlarge galaxy image for this factor
    rx = rx * scale
    ry = ry * scale

    # 30 is the minimum size:
    if rx < 30:
        rx = 30

    if ry < 30:
        ry = 30

    xx, yy, Rkron, theta, eps = Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle)

    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta + 90, eps, ncol, nrow)

    GetFits(image, imout, sky, xmin, xmax, ymin, ymax)

    if sigma:
        GetFits(sigma, sigout, 0, xmin, xmax, ymin, ymax)

    x_cor = xmin
    y_cor = ymin

    # computing new x, y positions for stamp
    x_small = xx - x_cor
    y_small = yy - y_cor

    return (x_small, y_small, x_cor, y_cor)


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
