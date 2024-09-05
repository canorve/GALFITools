#!/usr/bin/env python3

import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import GetExpTime
from galfitools.galin.std import Ds9ell2Kronell
from matplotlib.path import Path


def photDs9(ImageFile, RegFile, maskfile, zeropoint, sky):
    """computes the magnitude inside a DS9 region file

    Computes the magnitude inside the region defined by ellipse,
    box or polygon in DS9 region format.

    Parameters
    ----------
    ImageFile : str
                name of the image file
    RegFile : str
            name of the DS9 region file
    maskfile : str
               name of the mask image
    zeropoint : float
                magnitude zero point
    sky : float
         sky background in counts

    Returns
    -------
    mag : float
          magnitude measured from region
    exptime : float
            exposition time from image

    """
    if not os.path.exists(ImageFile):

        print("image filename does not exist!")
        sys.exit()

    if not os.path.exists(RegFile):
        print("%s: reg filename does not exist!" % (sys.argv[2]))
        sys.exit()

    (ncol, nrow) = GetAxis(ImageFile)

    exptime = GetExpTime(ImageFile)

    hdu = fits.open(ImageFile)

    Image = hdu[0].data

    hdu.close()

    # removing sky background from image

    Image = Image - sky

    v0 = []
    v1 = []
    v2 = []
    v3 = []
    v4 = []
    v5 = []

    tupVerts = []
    Pol = []

    f1 = open(RegFile, "r")

    lines = f1.readlines()

    f1.close()

    flag = False
    flagpoly = False
    # reading reg file
    for line in lines:

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True

        if b1[0] == "box":

            x0 = "box"
            x2 = x1[4:]
            flag = True

        if b1[0] == "polygon":

            polv = line.split(")")
            pol, ver = polv[0].split("(")

            points = ver.split(",")

            N = len(points)

            i = 0

            x = np.array([])
            y = np.array([])

            x = x.astype(float)
            y = y.astype(float)

            verts = []

            # make tuple for vertices
            while i < N:

                verts = verts + [
                    (round(float(points[i]) - 1), round(float(points[i + 1]) - 1))
                ]
                i = i + 2

            flagpoly = True

        if flag is True:

            x3 = p[4]
            x4 = x3[:-2]

            v0.append(x0)
            v1.append(float(x2) - 1)
            v2.append(float(p[1]) - 1)
            v3.append(float(p[2]))
            v4.append(float(p[3]))
            v5.append(float(x4))

            flag = False

        if flagpoly is True:

            Pol.append(pol)
            tupVerts.append(verts)
            flagpoly = False

    obj = np.array(v0)
    xpos = np.array(v1)
    ypos = np.array(v2)
    rx = np.array(v3)
    ry = np.array(v4)
    angle = np.array(v5)

    Pol = np.array(Pol)

    totFlux = 0

    if (maskfile == "none") or (maskfile == "None"):
        maskfile = None

    # mask file
    if maskfile:
        errmsg = "file {} does not exist".format(maskfile)
        assert os.path.isfile(maskfile), errmsg

        hdu2 = fits.open(maskfile)
        maskimage = hdu2[0].data
        maskb = np.array(maskimage, dtype=bool)
        invmask = np.logical_not(maskb)
        invmask = invmask * 1
        Image = Image * invmask
        hdu2.close()

    for idx, item in enumerate(obj):

        # get Flux

        if obj[idx] == "ellipse":

            ellFlux = FluxEllip(
                Image, xpos[idx], ypos[idx], rx[idx], ry[idx], angle[idx], ncol, nrow
            )

            totFlux = totFlux + ellFlux

        if obj[idx] == "box":

            boxFlux = FluxBox(
                Image, xpos[idx], ypos[idx], rx[idx], ry[idx], angle[idx], ncol, nrow
            )

            totFlux = totFlux + boxFlux

    for idx, item in enumerate(Pol):

        # converting ellipse from DS9 ell to kron ellipse

        if Pol[idx] == "polygon":

            polFlux = FluxPolygon(Image, tupVerts[idx], ncol, nrow)

            totFlux = totFlux + polFlux

    mag = -2.5 * np.log10(totFlux / exptime) + zeropoint

    return mag, exptime


def FluxEllip(Image, xpos, ypos, rx, ry, angle, ncol, nrow):
    """Gets the flux from an DS9 region ellipse in an image"""

    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos, ypos, rx, ry, angle)
    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)
    flux = FluxKron(Image, xx, yy, Rkron, theta, e, xmin, xmax, ymin, ymax)

    return flux


def FluxPolygon(Image, tupVerts, ncol, nrow):
    """Gets the flux from a DS9 region polygon in an image"""

    x, y = np.meshgrid(
        np.arange(ncol), np.arange(nrow)
    )  # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T
    p = Path(tupVerts)  # make a polygon

    grid = p.contains_points(points)
    mask = grid.reshape(nrow, ncol)  # now you have a mask with points inside a polygon

    flux = Image[mask].sum()

    return flux


def FluxBox(Image, xpos, ypos, rx, ry, angle, ncol, nrow):
    """Gets the flux from a DS9 region box in an image"""

    anglerad = angle * np.pi / 180
    beta = np.pi / 2 - anglerad

    lx = (rx / 2) * np.cos(anglerad) - (ry / 2) * np.cos(beta)
    lx2 = (rx / 2) * np.cos(anglerad) + (ry / 2) * np.cos(beta)
    ly = (rx / 2) * np.sin(anglerad) + (ry / 2) * np.sin(beta)
    ly2 = (rx / 2) * np.sin(anglerad) - (ry / 2) * np.sin(beta)

    v1x = round(xpos - lx)
    v1y = round(ypos - ly)

    v2x = round(xpos - lx2)
    v2y = round(ypos - ly2)

    v3x = round(xpos + lx)
    v3y = round(ypos + ly)

    v4x = round(xpos + lx2)
    v4y = round(ypos + ly2)

    Verts = [(v1x, v1y), (v2x, v2y), (v3x, v3y), (v4x, v4y)]

    x, y = np.meshgrid(
        np.arange(ncol), np.arange(nrow)
    )  # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T

    p = Path(Verts)  # make a polygon

    grid = p.contains_points(points)

    mask = grid.reshape(nrow, ncol)  # now you have a mask with points inside a polygon

    flux = Image[mask].sum()

    return flux


def FluxKron(imagemat, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    """This subroutine obtain the flux from a Kron ellipse delimited
    by box defined by: xmin, xmax, ymin, ymax
    """

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    q = 1 - ell
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1 : ymax + 1, xmin - 1 : xmax + 1]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x) ** 2 + (yell - y) ** 2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist <= dell
    flux = imagemat[ypos[mask], xpos[mask]].sum()

    return flux


#############################################################################
#   End of program
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
