#! /usr/bin/env python

import copy
import sys

import matplotlib.pyplot as plt
import numpy as np
from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import GetMe, GetReff
from scipy.optimize import bisect
from scipy.special import gamma, gammainc, gammaincinv


def getBT(
    galfitFile: str,
    dis: float,
    num_comp: int,
) -> float:
    """gets the Bulge to Total luminosity ratio


    Computes the Bulge to Total luminosity ratio from a
    a GALFIT file which contains a surface brightness model
    of two (or three) surface brightness components.

    It assumes that the first component in the GALFIT file
    is the bulge and the second the disk. If the model is
    composed of bulge, bar and disk, it takes the bulge as
    the first component, bar as second and the disk as the
    third component. Therefore sort the surface brightness
    components accordingly.

    If a bar is modeled it takes its luminosity as part
    of the bulge for the computation of the bulge/Total ratio

    If more than three component is found for the galaxy
    the programs ends and it does not compute the B/T.


    Parameters
    ----------
    galfitFile : str
                name of the GALFIT file
    dis : float
         maximum distance among components
    num_comp: int,
            number of component where the center will be take the center

    Returns
    -------

    bulge_total: float
    totmag : float
            total magnitude


    """

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    galcomps = SelectGal(galcomps, dis, num_comp)

    maskgal = galcomps.Active == 1

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute bulge_total ratio")
        print("exiting..")
        sys.exit(1)

    if N > 3:  # pragma: no cover
        print("maximum number of components = 3")
        print("exiting..")
        sys.exit(1)

    maskgal = comps.Active == 1

    comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

    totFlux = comps.Flux[maskgal].sum()

    totmag = -2.5 * np.log10(totFlux) + head.mgzpt

    mag_bulge = comps.Mag[maskgal][0]

    if N == 1:  # pragma: no cover
        bulge_total = 1

    if N == 2:  # pragma: no cover

        mag_disk = comps.Mag[maskgal][1]
        bulge_total = bulge_to_total(mag_bulge, mag_disk)

    if N == 3:  # pragma: no cover

        mag_bar = comps.Mag[maskgal][1]
        mag_disk = comps.Mag[maskgal][2]
        bulge_total = bulge_to_total(mag_bulge, mag_disk, mag_bar)

    return bulge_total, totmag, N


def bulge_to_total(mag_bulge, mag_disk, mag_bar=99):
    """
    Compute the bulge-to-total luminosity ratio (B/T) from bulge and disk magnitudes.

    Parameters
    ----------
    mag_bulge : float
        Apparent magnitude of the bulge.
    mag_disk : float
        Apparent magnitude of the disk.

    Returns
    -------
    float
        Bulge-to-total luminosity ratio.
    """
    L_bulge = np.power(10, -0.4 * mag_bulge)
    L_disk = np.power(10, -0.4 * mag_disk)
    L_bar = np.power(10, -0.4 * mag_bar)

    if mag_bar != 99:  # pragma: no cover
        L_bulge = L_bulge + L_bar

    bulge_total = L_bulge / (L_bulge + L_disk)

    return bulge_total


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
