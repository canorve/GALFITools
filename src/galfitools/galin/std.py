#! /usr/bin/env python


import numpy as np
from galfitools.galin.galfit import GalComps

from astropy.io import fits
import os.path


def MakeImage(newfits, sizex, sizey):
    """create a new blank FITS Image

    Parameters
    ----------
    newfits : str
              name of the new image
    sizex : int
            number of columns of the new image
    sizey : int
            number of rows of the new image

    Returns
    -------
    bool

    """

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        sp.run(
            [runcmd],
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True,
        )

    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex), dtype=np.float64)
    hdu.writeto(newfits, overwrite=True)

    return True


def GetSize(x, y, R, theta, ell, ncol, nrow):
    """Get the (x,y) coordinates that encompass the ellipse

    Parameters
    ----------
    x : float, x-center of ellipse
    y : float, y-center of ellipse
    R : float, major axis of ellipse
    theta : float, angular position of ellipse
    ell : float, ellipticity
    ncol : number of columns of the image
    nrow : number of rows of the image


    Returns
    -------
    xmin, xmax, ymin, ymax : Minimum and maximum coordinates
    that encompass the ellipse.

    """

    # k Check
    q = 1 - ell
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    constx = np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )
    consty = np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
    )

    xmin = x - constx
    xmax = x + constx
    ymin = y - consty
    ymax = y + consty

    mask = xmin < 1
    if mask.any():
        xmin = 1

    mask = xmax > ncol
    if mask.any():
        xmax = ncol - 1

    mask = ymin < 1
    if mask.any():
        ymin = 1

    mask = ymax > nrow
    if mask.any():
        ymax = nrow - 1

    return (round(xmin), round(xmax), round(ymin), round(ymax))


def Ds9ell2Kronell(xpos, ypos, rx, ry, angle):
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


    """

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle
        xx = xpos
        yy = ypos
    else:
        q = rx / ry
        e = 1 - q
        Rkron = ry
        theta = angle + 90
        xx = xpos
        yy = ypos

    return xx, yy, Rkron, theta, e


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


def Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle):
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


    """

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


def GetExpTime(Image):
    """Get exposition time
    from the image header

    """

    try:
        hdu = fits.open(Image)
        exptime = hdu[0].header["EXPTIME"]
        hdu.close()
    except Exception:
        exptime = 1
    return float(exptime)


def GetAxis(Image):
    """Get number of rows and columns from the image"""

    i = 0  # index indicated where the data is located

    # if checkCompHDU(Image):
    #    i = 1

    hdu = fits.open(Image)

    ncol = hdu[i].header["NAXIS1"]
    nrow = hdu[i].header["NAXIS2"]

    hdu.close()

    return ncol, nrow
