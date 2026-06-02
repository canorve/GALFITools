#!/usr/bin/env python3

import os.path
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import GetSize
from galfitools.galin.std import GetPmax
from galfitools.galin.std import GetExpTime
from galfitools.galin.std import Ds9ell2Kronell
from galfitools.galin.std import parse_ds9_ellipse
from matplotlib.path import Path


from galfitools.galout.getN import GetN

import matplotlib.pyplot as plt


def getNDs9(
    ImageFile,
    RegFile,
    maskfile,
    zeropoint,
    plate,
    sky,
    step=1,
    output="cowds9.png",
    dpival=200,
    plot=False,
):
    """computes the magnitude inside a DS9 region file to contruct
        the Curve of Growth

    Computes the magnitude inside the region defined by ellipse at
    different radius to compute the Curve of Growth. The ellipse
    radius is used as a limit for the curve of growth.

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
    plate: float
               plate scale
    sky : float
         sky background in counts

    step: float
         increase in radius for the magnitude integration


    output: plot output

    Returns
    -------
    mag : float
          magnitude measured from DS9 region
    sb: float
          surface brightness measured from DS9 region
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
    # reading reg file
    for raw_line in lines:

        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True

        if flag is True:

            x3 = p[4]
            x4 = x3[:-2]

            (obje, xe, ye, lxe, lye, anglee) = parse_ds9_ellipse(p)

            v0.append(obje)
            v1.append(xe)
            v2.append(ye)
            v3.append(lxe)
            v4.append(lye)
            v5.append(anglee)

            flag = False

    obj = np.array(v0)
    xpos = np.array(v1)
    ypos = np.array(v2)
    rx = np.array(v3)
    ry = np.array(v4)
    angle = np.array(v5)

    Pol = np.array(Pol)

    radFlux = 0
    Ntot = 0

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

            ellFlux, rad, Nell = FluxEllipStep(
                Image,
                xpos[idx],
                ypos[idx],
                rx[idx],
                ry[idx],
                angle[idx],
                ncol,
                nrow,
                step=step,
            )

            radFlux = radFlux + ellFlux
            Ntot = Ntot + Nell

    totFlux = radFlux[-1]

    fractionsa = np.arange(0.1, 0.5, 0.05)
    fractionsb = np.arange(0.55, 1, 0.05)
    fractions = np.concatenate(
        (fractionsa, fractionsb)
    )  # list containing the fraction of light

    target_fluxes = fractions * totFlux

    r_frac = np.interp(target_fluxes, radFlux, rad)

    print(f"R50: {r_frac[2]}, R80:{r_frac[3]}, R90:{r_frac[4]} ")

    mag = -2.5 * np.log10(radFlux / exptime) + zeropoint

    totmag = mag[-1]

    # estimating sersic indexs

    EffRad = r_frac[2]
    R = np.delete(r_frac, 2)
    F = np.delete(fractions, 2)
    ns = GetN().GalNs(EffRad, R, F)

    # plot:
    if plot:

        plt.clf()
        plt.plot(F, ns)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Fraction of light ")
        plt.ylabel("Sersic index")
        plt.savefig(output, dpi=dpival)
    return totmag, exptime, np.mean(ns), np.std(ns)


def FluxEllipStep(Image, xpos, ypos, rx, ry, angle, ncol, nrow, step=1):
    """Gets the flux from an DS9 region ellipse at diffent radius in an image"""

    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos, ypos, rx, ry, angle)

    Flux = np.array([])
    N = 0

    mask = np.array([])  # mask is already used in Image
    R = np.arange(1, Rkron + step, step)

    for r in R:

        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, r, theta, e, ncol, nrow)

        xpeak, ypeak = GetPmax(Image, mask, xmin, xmax, ymin, ymax)
        # flux, N = FluxKron(Image, xx, yy, r, theta, e, xmin, xmax, ymin, ymax)
        flux, N = FluxKron(Image, xpeak, ypeak, r, theta, e, xmin, xmax, ymin, ymax)
        Flux = np.append(Flux, flux)

    return Flux, R, N


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
    n_pix = np.count_nonzero(mask)

    return flux, n_pix


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

    n_pix = np.count_nonzero(mask)

    return flux, n_pix


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

    n_pix = np.count_nonzero(mask)

    return flux, n_pix


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
