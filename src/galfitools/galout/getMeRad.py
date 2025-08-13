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
    meanme : float
             mean of surface brightness at effective radius
             in mag/arcsec**2
    merad : float
             surface brightness at EffRad  in mag/arcsec**2
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
    if angle:
        theta = angle
    else:
        theta = galcomps.PosAng[maskgal][-1]

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    N = numComps(comps, "all")

    if N == 0:
        print("not enough number of components to compute Re")
        print("exiting..")
        sys.exit(1)

    totmag = getTotMag(head, comps)

    meanme = GetMe().MeanMe(totmag, rad * head.scale)
    me = GetMe().Me(head, comps, rad * head.scale, theta)

    return totmag, meanme, me, N, theta


def getTotMag(galhead: GalHead, comps: GalComps) -> float:

    maskgal = comps.Active == 1

    comps.Flux = 10 ** ((galhead.mgzpt - comps.Mag) / 2.5)

    totFlux = comps.Flux[maskgal].sum()

    totmag = -2.5 * np.log10(totFlux) + galhead.mgzpt

    return totmag
