#!/usr/bin/env python3


import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.MaskDs9 import GetAxis
from galfitools.mge.mge2galfit import Ds9ell2Kronellv2


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


def GetInfoEllip(regfile: str):
    """gets ellipse information from DS9 region files

    Parameters
    ----------
    regfile : str
                DS9 region file

    Returns
    -------
    obj : str
          object name
    xpos : float, x-center
    ypos : float, y-center
    rx : float, major or minor axis
    ry : float, minor or major axis
    angle : float, angular position


    returns 0, 0, 0, 0, 0, 0 if ellipse region was not found in file

    # repeated

    """

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

        # avoids ds9 regions with Area = 0
        if rx < 1:
            rx = 1
        if ry < 1:
            ry = 1

        return obj, xpos, ypos, rx, ry, angle

    else:
        print("ellipse region was not found in file. Exiting.. ")
        sys.exit()

    return 0, 0, 0, 0


def Ds9ell2Kronell(xpos: float, ypos: float, rx: float, ry: float, angle: float):
    """Converts DS9 ellipse parameters to geometrical parameters

    Parameters
    ----------
    obj : str
          object name
    xpos : float, x-center
    ypos : float, y-center
    rx : float, major or minor axis
    ry : float, minor or major axis
    angle : float, angular position

    Returns
    -------
    xx : float, x center position
    yy : float, y center positon
    Rkron : float, major axis
    theta: float, angular position measured from Y-axis
    e: float, ellipticity

    # repeated
    """

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle + 90
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
    """Get the (x,y) coordinates that encompass the ellipse

    Parameters
    ----------
    x : float, x-center of ellipse
    y : float, y-center of ellipse
    R : float, major axis of ellipse
    theta : float, angular position of ellipse. measured from X-axis
    ell : float, ellipticity
    ncol : number of columns of the image
    nrow : number of rows of the image


    Returns
    -------
    xmin, xmax, ymin, ymax : Minimum and maximum coordinates
    that encompass the ellipse.

    # repeated

    """

    # k Check
    q = 1 - ell
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    xmin = x - np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )

    xmax = x + np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )

    ymin = y - np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
    )

    ymax = y + np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
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
    """gets the peak coordinates

    Given an image, and optionally a mask, it identifies
    the (x, y) pixels where the maximum value is located.

    Parameters
    ----------
    image: 2D-array of the image
    mask: 2D-array of the mask
    xmin, xmax, ymin, ymax : int, int, int, int
                            coordinates of the image section
                            where the maximum will be obtained

    Returns
    -------
    (xpos, ypos) : (x, y) coordinates of the maximum

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
