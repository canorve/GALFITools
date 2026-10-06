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
    regFile="ellipse_imout.reg",
    idx=1,
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
    regFile: str
        name of the output Ds9 region file for the new image

    Output
    ------
        imout: new image
        regFile: new Ds9 ellipse in the new image


    Returns
    -------
    x_pos, y_pos: new pixel coordinate of the galaxy.

    x_cor, y_cor: x, y coordinates tranformation

    """

    (ncol, nrow) = GetAxis(image)

    obj, xpos, ypos, rxx, ryy, angle = GetInfoEllip(regfile)

    # enlarge galaxy image for this factor
    rx = rxx * scale
    ry = ryy * scale

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

    # writing output to DS9 File
    writeDs9Ellipse(
        idx,
        x_small,
        y_small,
        rxx,
        ryy,
        angle,
        region_file=regFile,
    )

    return (x_small, y_small, x_cor, y_cor)


def writeDs9Ellipse(objid, xpos, ypos, rx, ry, angle, region_file="ellipse.reg"):
    """Write an ellipse to a DS9 region file.

    Parameters
    ----------
    objid : str or int
        Object identifier, displayed as a label.
    xpos, ypos : float
        Center coordinates in DS9 image pixels (1-based).
    rx, ry : float
        Ellipse semi-axis lengths in pixels.
    angle : float
        Rotation angle in degrees.
    region_file : str or path-like, optional
        Output filename. An existing file is overwritten.
    """
    with open(region_file, "w", encoding="utf-8") as output:
        output.write("# Region file format: DS9 version 4.1\n")
        output.write("image\n")
        output.write(
            f"ellipse({xpos},{ypos},{rx},{ry},{angle}) " f"# text={{{objid}}}\n"
        )


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
