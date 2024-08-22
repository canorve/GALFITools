#!/usr/bin/env python3

import os.path
import sys

import numpy as np
from astropy.io import fits
from astropy.io.fits import CompImageHDU
from matplotlib.path import Path


def maskDs9(
    MaskFile: str,
    RegFile: str,
    fill: float,
    image: str,
    bor_flag: bool,
    borValue: int,
    skymean=None,
    skystd=None,
) -> None:
    """Creates masks from DS9 regions

    given DS9 regions such as box, ellipse or polygon it
    creates those regions as masks in a mask image for GALFIT.
    If the mask image does not exists, it creates one.

    It also removes unwanted regions in the original image
    such as saturated regions by providing mean and standard
    deviation of sky.

    Parameters
    ----------

    MaskFile : str
              name of the mask image file
    RegFile : str
             name of the file containing the DS9 regions
    fill : float
            value to fill withing DS9 regions
    image : str
            new of the image file to take the size of the matrix
            to create a new mask file if it does not exists.
    bor_flag : False, optional
               if True, it will mask the border of the image. This is
               for those region where the image matrix is larger than the
               data matrix, e.g.  Hubble images
    borValue : float
               value of the border
    skymean : None, optional
              mean of the sky value to be substituted inside of
              the DS9 region
    skystd : None, optional
             standard deviation of the sky value to be substituted
             inside of the DS9 region

    Returns
    -------
    None

    """
    if not os.path.exists(MaskFile):

        print("%s: image filename does not exist!" % (MaskFile))
        print("Creating a new image file ")

        hdu = fits.PrimaryHDU()

        if image:
            (nx, ny) = GetAxis(image)
        else:
            nx = input("enter numbers of pixels in X ")
            ny = input("enter numbers of pixels in Y ")

            nx = np.int64(nx)
            ny = np.int64(ny)

        Image = np.zeros([ny, nx])
        hdu.data = Image
        hdu.writeto(MaskFile, overwrite=True)

    (ncol, nrow) = GetAxis(MaskFile)

    i = 0  # index of data

    hdu = fits.open(MaskFile)

    Image = hdu[i].data

    if not os.path.exists(RegFile):
        print("%s: reg filename does not exist!" % (sys.argv[2]))
        sys.exit(1)

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

        line = line.split("#")
        line = line[0]

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

    # avoids ds9 regions with Area = 0
    maskrx = rx < 1
    if maskrx.any():
        rx[maskrx] = 1

    maskry = ry < 1
    if maskry.any():
        ry[maskry] = 1
    #

    Pol = np.array(Pol)

    for idx, item in enumerate(obj):

        # converting ellipse from DS9 ell to kron ellipse

        if obj[idx] == "ellipse":

            Image = MakeEllip(
                Image,
                fill,
                xpos[idx],
                ypos[idx],
                rx[idx],
                ry[idx],
                angle[idx],
                ncol,
                nrow,
                skymean=skymean,
                skystd=skystd,
            )

        if obj[idx] == "box":

            Image = MakeBox(
                Image,
                fill,
                xpos[idx],
                ypos[idx],
                rx[idx],
                ry[idx],
                angle[idx],
                ncol,
                nrow,
                skymean=skymean,
                skystd=skystd,
            )

    for idx, item in enumerate(Pol):

        # converting ellipse from DS9 ell to kron ellipse

        if Pol[idx] == "polygon":

            Image = MakePolygon(
                Image, fill, tupVerts[idx], ncol, nrow, skymean=skymean, skystd=skystd
            )

    if image:

        hduim = fits.open(image)
        dataImage = hduim[0].data

        # masking the border in case:
        bor_val = 100
        if bor_flag:
            # print("masking the border")
            bor_mask = dataImage == borValue

            if bor_mask.any():
                Image[bor_mask] = bor_val

        hduim.close()

    # writing mask file

    hdu.data = Image
    hdu.writeto(MaskFile, overwrite=True)
    hdu.close()


def MakeEllip(
    Image, fill, xpos, ypos, rx, ry, angle, ncol, nrow, skymean=None, skystd=None
):
    """Draw an ellipse in an image

    Parameters
    ----------
    Image : ndarray
            the image matrix array
    fill : float
            value to be filled inside the DS9 region
    xpos, ypos : float, float
                coordinates of the center
    rx : float
        major axis or minor axis of the ellipse

    ry : float
        minor axis or major axis of the ellipse

    angle : float
            angular position of the ellipse

    ncol, nrow : int, int
              size of the image
    skymean : float
              mean of the sky background
    skystd : float
             standard deviation of sky


    Returns
    -------

    Image : ndarray
            the image with the ellipse

    """

    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos, ypos, rx, ry, angle)
    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)
    Image = MakeKron(
        Image,
        fill,
        xx,
        yy,
        Rkron,
        theta,
        e,
        xmin,
        xmax,
        ymin,
        ymax,
        ncol,
        nrow,
        skymean=skymean,
        skystd=skystd,
    )

    return Image


def MakePolygon(Image, fill, tupVerts, ncol, nrow, skymean=None, skystd=None):
    """Make a polygon in an image

    Parameters
    ----------
    Image : ndarray
           the image matrix array
    fill : float
           value to be filled inside the DS9 region
    tupVerts : tuple
              vertices of the polygon

    ncol, nrow : int, int
                size of the image
    skymean : float
              mean of the sky background
    skystd : float
            standard deviation of sky


    Returns
    -------
    Image : ndarray
            the new image with the polygon

    """

    x, y = np.meshgrid(
        np.arange(ncol), np.arange(nrow)
    )  # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x, y)).T
    p = Path(tupVerts)  # make a polygon

    grid = p.contains_points(points)
    mask = grid.reshape(nrow, ncol)  # now you have a mask with points inside a polygon

    if skymean:

        sky = np.random.normal(skymean, skystd, (nrow, ncol))
        Image[mask] = sky[mask]

    else:

        Image[mask] = fill

    return Image


def MakeBox(
    Image, fill, xpos, ypos, rx, ry, angle, ncol, nrow, skymean=None, skystd=None
):
    """Make a box in an image

    Parameters
    ----------

    Image : ndarray
           the image matrix array
    fill : float
           value to be filled inside the DS9 region
    xpos, ypos : float, float
                coordinates of the center
    rx : float
        x-longitude of the box

    ry : float
        y-longitude of the box

    angle : float
            angular position of the box

    ncol, nrow : int, int
                size of the image
    skymean : float
              mean of the sky background
    skystd : float
            standard deviation of sky




    Returns
    -------
    Image : ndarray
            the image with the box region

    """

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

    if skymean:

        sky = np.random.normal(skymean, skystd, (nrow, ncol))
        Image[mask] = sky[mask]

    else:

        Image[mask] = fill

    return Image


def GetAxis(Image):
    """Get number of rows and columns from the image"""

    i = 0  # index indicated where the data is located

    if checkCompHDU(Image):
        i = 1

    hdu = fits.open(Image)

    ncol = hdu[i].header["NAXIS1"]
    nrow = hdu[i].header["NAXIS2"]

    hdu.close()

    return ncol, nrow


def checkCompHDU(file):
    """check if fits file is a CompImageHDU

    Notes
    -----
    function not used anymore

    """
    flag = False
    hdul = fits.open(file)

    for i, hdu in enumerate(hdul):

        if isinstance(hdu, CompImageHDU):
            flag = True

    hdul.close()

    return flag


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

    # repeated

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


def MakeKron(
    imagemat,
    idn,
    x,
    y,
    R,
    theta,
    ell,
    xmin,
    xmax,
    ymin,
    ymax,
    ncol,
    nrow,
    skymean=None,
    skystd=None,
):
    """This creates a ellipse in an image

    This function creates an ellipse and fills the pixels inside it
    with the value specified by idn. The ellipse is created within
    a box delimited by xmin, xmax, ymin, and ymax.

    Parameters
    ----------
    imagemat : ndarray, the image matrix
    idn : int the value to fill the ellipse
    x, y : position of the ellipse's center
    R : ellipse major axis
    theta : angular position of the ellipse measured from X-axis
    ell : ellipticity of the ellipse
    xmin, xmax, ymin, ymax : int, int, int, int
            box delimitation of the ellipse

    Returns
    -------
    imagemat : the image with the new ellipse

    # repeated

    """

    # Check

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

    if skymean:

        sky = np.random.normal(skymean, skystd, (nrow, ncol))
        imagemat[ypos[mask], xpos[mask]] = sky[ypos[mask], xpos[mask]]

    else:

        imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


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

    # repeated
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
