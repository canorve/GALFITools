#! /usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy
from galfitools.galin.galfit import (
    GalComps,
    Galfit,
    GalHead,
    GetRadAng,
    SelectGal,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import GetMe

from scipy.interpolate import UnivariateSpline
from scipy.optimize import bisect, newton
from scipy.special import gamma, gammainc, gammaincinv


def getMeRad(
    galfitFile: str, dis: int, rad: float, angle: float, num_comp: int
) -> float:
    """Obtains the surface brightness at a given radius from a set of Sersic functions.

    Given a model composed of multiple Sersic functions,
    it returns the surface brightness at given radius

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    rad: float
         Radius at which the surface brightness is computed
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    num_comp: int
            Number of component from which the center of all
            components will be determined.

    Returns
    -------

    totmag : float
             total magnitude
    meanmerad : float
             mean of surface brightness at rad
             in mag/arcsec**2
    merad : float
             surface brightness at rad in mag/arcsec**2
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    """

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    galcomps = SelectGal(galcomps, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal = galcomps.Active == 1
    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = galcomps.PosAng[maskgal][-1]

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute Re")
        print("exiting..")
        sys.exit(1)

    totmag = getTotMag(head, comps)

    meanme = GetMe().MeanMe(totmag, rad * head.scale)

    meanmerad = getMeanRad(head, comps, rad * head.scale)  # it does not work
    merad = GetMe().Me(head, comps, rad * head.scale, theta)

    return totmag, meanmerad, merad, N, theta


def getTotMag(galhead: GalHead, comps: GalComps) -> float:
    """Returns the total magnitude of a set of sersics"""

    maskgal = comps.Active == 1

    comps.Flux = 10 ** ((galhead.mgzpt - comps.Mag) / 2.5)

    totFlux = comps.Flux[maskgal].sum()

    totmag = -2.5 * np.log10(totFlux) + galhead.mgzpt

    return totmag


def getMeanRad(galhead: GalHead, comps: GalComps, rad: float) -> float:
    """
    mean surface brightness at radius rad
    Note: As contrary to surface brightness at rad it
    does not assume angle direction
    """

    maskgal = comps.Active == 1

    comps.Flux = 10 ** ((galhead.mgzpt - comps.Mag) / 2.5)

    comps.Rad = comps.Rad * galhead.scale  # converting to arcsec
    FluxR = GetFtotser(
        rad,
        comps.Flux[maskgal],
        comps.Rad[maskgal],
        comps.Exp[maskgal],
        comps.AxRat[maskgal],
        comps.PosAng[maskgal],
    )
    ####

    meanIerad = FluxR / (np.pi * rad**2)

    meanmerad = -2.5 * np.log10(meanIerad) + galhead.mgzpt

    return meanmerad


def GetFtotser(
    R: float,
    flux: list,
    rad: list,
    n: list,
    q: list,
    pa: list,
) -> float:
    """partial sersic flux computed from zero up to R for a set of sersics"""

    ftotR = Fser(R, flux, rad, n, q, pa)

    return ftotR.sum()


def Fser(
    R: float,
    Flux: list,
    Re: list,
    n: list,
    q: list,
    pa: list,
) -> float:
    """partial sersic flux computed from zero up to R for a single sersic"""

    k = gammaincinv(2 * n, 0.5)

    X = k * (R / Re) ** (1 / n)

    Fr = Flux * gammainc(2 * n, X)

    return Fr
