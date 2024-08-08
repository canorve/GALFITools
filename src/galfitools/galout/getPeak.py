#!/usr/bin/env python3


import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.MaskDs9 import GetAxis


def getPeak(image: str, regfile: str, center: bool, maskfile: str) -> None:
    """return the coordinates of the peak given an DS9 region file"""

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


def GetInfoEllip(regfile):

    if not os.path.exists(regfile):
        print("%s: reg filename does not exist!" % (regfile))
        sys.exit()

    f1 = open(regfile, "r")

    lines = f1.readlines()

    f1.close()

    flag = False
    found = False

    # reading reg file
    for line in lines:

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True
            found = True

        if flag is True:
            x3 = p[4]
            x4 = x3[:-2]

            v0 = x0

            v1 = float(x2)
            v2 = float(p[1])
            v3 = float(p[2])
            v4 = float(p[3])
            v5 = float(x4)

            flag = False

    if found:
        obj = v0
        xpos = v1
        ypos = v2
        rx = v3
        ry = v4
        angle = v5

        return obj, xpos, ypos, rx, ry, angle
    else:
        print("ellipse region was not found in file. Exiting.. ")
        sys.exit()

    return 0, 0, 0, 0


def Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle):

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle - 90
        xx = xpos
        yy = ypos
    else:
        q = rx / ry
        e = 1 - q
        Rkron = ry
        theta = angle  # + 90
        xx = xpos
        yy = ypos

    return xx, yy, Rkron, theta, e


def GetSize(x, y, R, theta, ell, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    # k Check
    q = 1 - ell
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    xmin = x - np.sqrt(
        (R ** 2) * (np.cos(theta)) ** 2 + (bim ** 2) * (np.sin(theta)) ** 2
    )

    xmax = x + np.sqrt(
        (R ** 2) * (np.cos(theta)) ** 2 + (bim ** 2) * (np.sin(theta)) ** 2
    )

    ymin = y - np.sqrt(
        (R ** 2) * (np.sin(theta)) ** 2 + (bim ** 2) * (np.cos(theta)) ** 2
    )

    ymax = y + np.sqrt(
        (R ** 2) * (np.sin(theta)) ** 2 + (bim ** 2) * (np.cos(theta)) ** 2
    )

    mask = xmin < 1
    if mask.any():
        xmin = 1

    mask = xmax > ncol
    if mask.any():
        xmax = ncol

    mask = ymin < 1
    if mask.any():
        ymin = 1

    mask = ymax > nrow
    if mask.any():
        ymax = nrow

    return (xmin, xmax, ymin, ymax)


def GetPmax(image, mask, xmin, xmax, ymin, ymax):

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
