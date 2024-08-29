#! /usr/bin/env python

import os
import os.path

import numpy as np
from astropy.io import fits
from galfitools.galin.MaskDs9 import GetAxis
from galfitools.galin.std import ds9satbox


def makeMask(
    sexfile: str, image: str, maskfile: str, scale: float, satfileout: str
) -> None:
    """Creates a mask file from a catalog of SExtractor

    It creates a mask file for GALFIT using information from a
    SExtractor catalog. It includes masking of saturated regions.

    Parameters
    ----------
    sexfile : str
            name of the Sextractor catalog

    image: str,
            name of the image file

    maskfile: str,
            name of the mask file
    scale: float,
            Scale factor by which the ellipse will be enlarged or diminished.
    satfileout: str
            DS9 region file where the saturation regions will be indicated

    Returns
    -------
    None

    """

    sexarsort = "sexsort.cat"

    satscale = 1
    satoffset = 0

    print("Creating masks....\n")

    (NCol, NRow) = GetAxis(image)

    Total = CatArSort(sexfile, scale, sexarsort, NCol, NRow)

    print("Creating sat region files....\n")

    if not os.path.exists(satfileout):
        print("Saturation file not found. Creating one")
        satfileout = "ds9sat.reg"
        ds9satbox(
            satfileout, sexfile, satscale, satoffset
        )  # crea archivo  Saturacion reg

    # segmentation mask

    MakeImage(maskfile, NCol, NRow)

    MakeMask(maskfile, sexarsort, scale, 0, satfileout)  # offset set to 0
    MakeSatBox(maskfile, satfileout, Total + 1, NCol, NRow)  # make sat region


'''
def ds9satbox(satfileout, sexcat, satscale, satoffset):
    """Creates a DS9 file which contains regions with bad
        saturated regions

    Parameters
    ----------
    satfileout : str, DS9 region output file
    sexcat : str, SExtractor catalog
    satscale : float,
             Scale factor by which the saturation region will be
             enlarged or diminished.
    satoffset: float,
            constant to be added to the size of saturated region

    Returns
    -------
    None


    # repeated
    """

    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag
    check = 0

    f_out = open(satfileout, "w")

    (
        N,
        Alpha,
        Delta,
        X,
        Y,
        Mg,
        Kr,
        Fluxr,
        Isoa,
        Ai,
        E,
        Theta,
        Bkgd,
        Idx,
        Flg,
    ) = np.genfromtxt(sexcat, delimiter="", unpack=True)

    line = "image \n"
    f_out.write(line)

    for idx, item in enumerate(N):

        bi = Ai[idx] * (1 - E[idx])

        Theta[idx] = Theta[idx] * np.pi / 180  # rads!!!

        Rkronx = satscale * 2 * Ai[idx] * Kr[idx] + satoffset
        Rkrony = satscale * 2 * bi * Kr[idx] + satoffset

        if Rkronx == 0:
            Rkronx = 1

        if Rkrony == 0:
            Rkrony = 1

        check = CheckFlag(
            Flg[idx], flagsat, maxflag
        )  # check if object has saturated regions

        if check:

            line = "box({0},{1},{2},{3},0) # color=red move=0 \n".format(
                X[idx], Y[idx], Rkronx, Rkrony
            )
            f_out.write(line)

            line2a = "point({0},{1}) ".format(X[idx], Y[idx])
            line2b = '# point=boxcircle font="times 10 bold" text={{ {0} }} \n'.format(
                N[idx]
            )
            line2 = line2a + line2b
            f_out.write(line2)

    f_out.close()
'''


def MakeMask(maskimage, catfile, scale, offset, regfile):
    """Creates ellipse masks for every object of the SExtractor catalog

    Parameters
    ----------
    maskimage: str
            name of the mask file. This file should already exists
    catfile: str,
            SExtractor catalog
    scale: float,
            Scale factor by which the ellipse will be enlarged or diminished.
    offset: float
            constant to be added to the ellipse size
    regfile: str
            DS9 region file containing the saturated region.

    Returns
    -------
    None

    # repeated
    """

    checkflag = 0
    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    (
        n,
        alpha,
        delta,
        xx,
        yy,
        mg,
        kr,
        fluxrad,
        ia,
        ai,
        e,
        theta,
        bkgd,
        idx,
        flg,
        sxmin,
        sxmax,
        symin,
        symax,
        sxsmin,
        sxsmax,
        sysmin,
        sysmax,
    ) = np.genfromtxt(catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    # print("Creating Masks for sky \n")

    Rkron = scale * ai * kr + offset

    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1

    i = 0  # index indicated where the data is located

    hdu = fits.open(maskimage)
    img = hdu[i].data

    print("Creating ellipse masks for every object \n")

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = CheckFlag(flg[idx], flagsat, maxflag)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
        regflag = CheckSatReg(xx[idx], yy[idx], regfile, Rkron[idx], theta[idx], e[idx])

        if (checkflag is False) and (regflag is False):

            img = MakeKron(
                img,
                n[idx],
                xx[idx],
                yy[idx],
                Rkron[idx],
                theta[idx],
                e[idx],
                sxsmin[idx],
                sxsmax[idx],
                sysmin[idx],
                sysmax[idx],
            )

        # elif(checkflag == True or regflag == True):

    print("ignoring objects where one or more pixels are saturated \n")

    hdu[i].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
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

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    q = 1 - ell
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1 : ymax, xmin - 1 : xmax]

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

    mask = dist < dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


def MakeSatBox(maskimage, region, val, ncol, nrow):
    """Creates saturated regions in the mask

    Parameters
    ----------
    maskimage : str, name of the mask file
    region : str, DS9 box region
    val : int, value to be filled inside the saturated region
    ncol : int, number of columns of the image
    nrow : int, number of rows of the image

    Returns
    -------
    bool
    # repeated
    """

    # k Check

    # 	fileflag=1
    i = 0  # index where data is located
    hdu = fits.open(maskimage)
    img = hdu[i].data

    with open(region) as f_in:

        next(f_in)

        # All lines including the blank ones
        lines = (line.rstrip() for line in f_in)
        lines = (line.split("#", 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            words = line.split(" ")
            if words[0] != "image" and words[0] != "physical" and words[0] != "global":

                (box, info) = line.split("(")

                if box == "box":

                    (xpos, ypos, xlong, ylong, trash) = info.split(",")

                    xpos = float(xpos)
                    ypos = float(ypos)
                    xlong = float(xlong)
                    ylong = float(ylong)

                    xlo = xpos - xlong / 2
                    xhi = xpos + xlong / 2

                    ylo = ypos - ylong / 2
                    yhi = ypos + ylong / 2

                    xlo = int(xlo)
                    xhi = int(xhi)

                    ylo = int(ylo)
                    yhi = int(yhi)

                    if xlo < 1:

                        xlo = 1

                    if xhi > ncol:

                        xhi = ncol

                    if ylo < 1:

                        ylo = 1

                    if yhi > nrow:

                        yhi = nrow

                    img[ylo - 1 : yhi, xlo - 1 : xhi] = val

    hdu[i].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeImage(newfits, sizex, sizey):
    """Creates a new blank Image"""
    # repeated
    # k Check
    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True


def CatArSort(SexCat, scale, SexArSort, NCol, NRow):
    """Sorts the SExtractor catalog by area, from largest to smallest.

    Parameters
    ----------
    SexCat : str, name of the SExtractor catalog
    scale : float,
    SexArSort : new output SExtractor catalog
    NCol : number of columns of the image
    NRow : number of rows of the image


    Returns
    -------
    number of objects

    # repeated
    """

    # sort the sextractor
    # catalog by magnitude,
    # get sizes for objects
    # and write it in a new file

    print("Sorting and getting sizes for objects \n")

    (
        n,
        alpha,
        delta,
        xx,
        yy,
        mg,
        kr,
        fluxrad,
        ia,
        ai,
        e,
        theta,
        bkgd,
        idx,
        flg,
    ) = np.genfromtxt(SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    #    ai = ai.astype(float)
    #    kr = kr.astype(float)

    #    scale = scale.astype(float)

    Rkron = scale * ai * kr

    Rwsky = scale * ai * kr + 10 + 20

    #   considering to use only  KronScale instead of SkyScale
    #    Rwsky = parvar.KronScale * ai * kr + parvar.Offset + parvar.SkyWidth

    Bim = (1 - e) * Rkron

    Area = np.pi * Rkron * Bim * (-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, NCol, NRow)

    (sxsmin, sxsmax, sysmin, sysmax) = GetSize(xx, yy, Rwsky, theta, e, NCol, NRow)

    f_out = open(SexArSort, "w")

    index = Area.argsort()
    for i in index:

        line1 = "{} {} {} {} {} {} {} {} {} {} ".format(
            n[i],
            alpha[i],
            delta[i],
            xx[i],
            yy[i],
            mg[i],
            kr[i],
            fluxrad[i],
            ia[i],
            ai[i],
        )
        line2 = "{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
            e[i],
            theta[i],
            bkgd[i],
            idx[i],
            flg[i],
            int(np.round(sxmin[i])),
            int(np.round(sxmax[i])),
            int(np.round(symin[i])),
            int(np.round(symax[i])),
            int(np.round(sxsmin[i])),
            int(np.round(sxsmax[i])),
            int(np.round(sysmin[i])),
            int(np.round(sysmax[i])),
        )

        line = line1 + line2

        f_out.write(line)

    f_out.close()

    return len(n)


def GetSize(x, y, R, theta, ell, ncol, nrow):
    """This function retrieves the minimum and
    maximum pixel coordinates that enclose the ellipse.

    Parameters
    ----------
    x, y : int, int, position of the ellipse's center
    R : float, ellipse major axis
    theta : float, angular position of the ellipse measured from X-axis
    ell : float, ellipticity of the ellipse
    ncol : number of columns of the image
    nrow : number of rows of the image


    Returns
    -------

    xmin, xmax, ymin, ymax : int, int, int, int
            box delimitation of the ellipse

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
        if isinstance(xmin, np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol

    if mask.any():
        if isinstance(xmax, np.ndarray):
            xmax[mask] = ncol
        else:
            xmax = ncol

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin, np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax, np.ndarray):
            ymax[mask] = nrow
        else:
            ymax = nrow

    return (xmin, xmax, ymin, ymax)


def CheckFlag(val, check, maxx):
    """Check if SExtractor flag contains the check flag

    This function is useful to check if
    a object sextractor flag contains saturated region

    Parameters
    ----------
    val: SExtractor flag of the object
    check: SExtractor flag value to check
    maxx: maximum value of the flag

    Returns
    -------
    bool, returns True if found

    # repeated

    """

    flag = False
    mod = 1

    while mod != 0:

        res = int(val / maxx)

        if maxx == check and res == 1:

            flag = True

        mod = val % maxx

        val = mod
        maxx = maxx / 2

    return flag


def CheckSatReg(x, y, filein, R, theta, ell):
    """Check if object is inside of saturated region.


    Parameters
    ----------
    (x, y): int, int, center's coordinates in pixels
    filein : str, file containing the saturated regions
    R: float, major axis
    theta: float, angular position of the ellipse
    ell: float, ellipticity of the object

    Returns
    -------
    bool : it returns True if at least one pixel is inside


    # repeated

    """

    q = 1 - ell

    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    flag = False

    with open(filein) as f_in:

        lines = (line.rstrip() for line in f_in)  # All lines including the blank ones
        lines = (line.split("#", 1)[0] for line in lines)  # remove comments
        lines = (
            line.rstrip() for line in lines
        )  # remove lines containing only comments
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            words = line.split(" ")
            if words[0] != "image" and words[0] != "physical" and words[0] != "global":

                (box, info) = line.split("(")

                if box == "box":

                    (xpos, ypos, xlong, ylong, trash) = info.split(",")

                    xpos = float(xpos)
                    ypos = float(ypos)
                    xlong = float(xlong)
                    ylong = float(ylong)

                    xlo = xpos - xlong / 2
                    xhi = xpos + xlong / 2

                    ylo = ypos - ylong / 2
                    yhi = ypos + ylong / 2

                    dx = xpos - x
                    dy = ypos - y

                    landa = np.arctan2(dy, dx)

                    if landa < 0:
                        landa = landa + 2 * np.pi

                    landa = landa - theta

                    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

                    xell = (
                        x
                        + R * np.cos(angle) * np.cos(theta)
                        - bim * np.sin(angle) * np.sin(theta)
                    )
                    yell = (
                        y
                        + R * np.cos(angle) * np.sin(theta)
                        + bim * np.sin(angle) * np.cos(theta)
                    )

                    if (xell > xlo and xell < xhi) and (yell > ylo and yell < yhi):

                        flag = True
                        break

    return flag
