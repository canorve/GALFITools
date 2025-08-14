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

    if N == 0:
        print("not enough number of components to compute bulge_total ratio")
        print("exiting..")
        sys.exit(1)

    if N > 3:
        print("maximum number of components = 3")
        print("exiting..")
        sys.exit(1)

    maskgal = comps.Active == 1

    comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

    totFlux = comps.Flux[maskgal].sum()

    totmag = -2.5 * np.log10(totFlux) + head.mgzpt

    mag_bulge = comps.Mag[maskgal][0]

    if N == 1:
        bulge_total = 1

    if N == 2:

        mag_disk = comps.Mag[maskgal][1]
        bulge_total = bulge_to_total(mag_bulge, mag_disk)

    if N == 3:

        mag_bar = comps.Mag[maskgal][1]
        mag_disk = comps.Mag[maskgal][2]
        bulge_total = bulge_to_total(mag_bulge, mag_disk, mag_bar)

    return bulge_total, totmag


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

    if mag_bar != 99:
        L_bulge = L_bulge + L_bar

    bulge_total = L_bulge / (L_bulge + L_disk)

    return bulge_total


class GetN:
    """Computes the Sersic index from photometric parameters



    Methods
    -------
    GalNs : computes the Sersic index from effective radius and
            other fraction of light radius.

    MeMeanMe : computes the Sersic index from surface brightness
                at Re and mean surface brightness at Re.


    """

    def GalNs(self, EffRad: float, EffRadfrac: list, F: list):

        sers = np.array([])

        for idx, f in enumerate(F):

            n = self.ReRfrac(EffRad, EffRadfrac[idx], f)

            sers = np.append(sers, n)

        return sers

    def MeMeanMe(self, me: float, meanme: float) -> float:

        a = 0.1
        b = 12

        Sersic = self.solveSer(a, b, me, meanme)

        return Sersic

    def solveSer(self, a: float, b: float, me: float, meanme: float) -> float:
        "return the sersic index. It uses Bisection"

        try:
            N = bisect(self.funMeMeanMe, a, b, args=(me, meanme))
        except Exception:
            print("unable to solve equation in the given range. Setting n to 99")
            N = 99

        return N

    def funMeMeanMe(self, n: float, me: float, meanme: float) -> float:

        # k = gammaincinv(2 * n, 0.5)

        fn = self.Fn(n)

        result = me - meanme - 2.5 * np.log10(fn)

        return result

    def Fn(self, n: float) -> float:

        k = gammaincinv(2 * n, 0.5)

        fn = ((n * np.exp(k)) / (k ** (2 * n))) * gamma(2 * n)

        return fn

    def ReRfrac(self, Re: float, Rfrac: float, frac: float) -> float:

        a = 0.1
        b = 12

        Sersic = self.solveSerRe(a, b, Re, Rfrac, frac)

        return Sersic

    def solveSerRe(
        self, a: float, b: float, Re: float, Rfrac: float, frac: float
    ) -> float:
        "return the sersic index. It uses Bisection"

        try:
            N = bisect(self.funReRfrac, a, b, args=(Re, Rfrac, frac))
        except Exception:
            print("unable to solve equation in the given range. Setting n to 99")
            N = 99

        return N

    def funReRfrac(self, n: float, Re: float, Rfrac: float, frac: float) -> float:

        k = gammaincinv(2 * n, 0.5)

        x = k * (Rfrac / Re) ** (1 / n)

        result = frac * gamma(2 * n) - gamma(2 * n) * gammainc(2 * n, x)

        return result

    def MeM0(self, me: float, m0: float) -> float:
        """Uses me and m0 to compute n. It is not very realiable
        and takes longer than the other two methods"""
        a = 0
        b = 45

        K = self.solveKm0(a, b, me, m0)

        a = 0.2
        b = 40

        Sersic = self.solveSerK(a, b, K)

        return Sersic

    def solveKm0(self, a: float, b: float, me: float, m0: float) -> float:
        "return the sersic index. It uses Bisection"

        try:
            K = bisect(self.funMeM0, a, b, args=(me, m0))
        except Exception:
            print("solution not found for the given range")
            K = 0

        return K

    def funMeM0(self, K: float, me: float, m0: float) -> float:

        result = me - m0 - 2.5 * K / np.log(10)

        return result

    def solveSerK(self, a: float, b: float, k: float) -> float:

        try:
            sersic = bisect(self.funK, a, b, args=(k))
        except Exception:
            print("unable to solve equation in the given range. Setting n to 99")

        return sersic

    def funK(self, n: float, k: float) -> float:

        result = gammaincinv(2 * n, 0.5) - k

        return result


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
