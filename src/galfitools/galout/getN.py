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


def getN(
    galfitFile: str,
    dis: float,
    frac: float,
    angle: float,
    num_comp: int,
    plot: bool,
    const=0,
) -> float:
    """gets the Sersic index

    Assuming the galaxy is physically composed of a single
    Sersic profile, this function computes the Sersic index
    using two methods: (1) by comparing the mean surface brightness
    at the effective radius to the surface brightness at the
    effective radius, and (2) by using two light fraction
    radii (effective radius and any other radius). Both methods
    employ bisection.

    Although it may seem contradictory to derive the Sersic index
    from a fit involving multiple Sersic components, this approach
    is useful when the galaxy is modeled using MGE (Multi-Gaussian
    Expansion) and the Sersic index needs to be estimated.



    Parameters
    ----------
    galfitFile : str
                name of the GALFIT file
    dis : float
         maximum distance among components
    frac: float,
            fraction of light.
    angle: float,
        Angle of the major axis of the galaxy. If None, it uses
        the angle of the last component
    num_comp: int,
            number of component where the center will be take the center
    plot: bool,
         If True, makes a plot of Sersic index vs. fraction of light radius
         used for computation.
    const : float, optional
         Substract constant from plot (for visualization purposes only)


    Returns
    -------

    sersic : float
    float
        Sersic index mean of the fraction of light radius
    float
        Sersic index standard deviation mean of the fraction of
        light radius method
    totmag : float
            total magnitude
    N : int
        total number of components
    theta : Angular position indicating the direction along
            which the effective radius and surface brightness
            at the effective radius are computed.


    Warnings
    --------

    The computation of the Sersic index assumes
    that all Sersic functions share the same center.
    Therefore, use the `dis` tolerance variable with
    caution. If the distance between components
    is significant, the result will be less accurate.

    Notes
    -----
    The fraction light method, gets the Sersic index using
    two fraction of light radius, for example it uses the
    Re and R90 (50% of light radius and 90% of light radius)
    Hence this method uses different fraction of light radius
    and Re, and it returns the mean and standard deviation of
    all the estimations of the Sersic index.

    The fraction light method calculates the Sersic index
    using two light fraction radii, such as the effective
    radius (Re) and another radius enclosing some fraction
    of the light. This method uses various light fraction
    radii and Re, and returns the mean and standard
    deviation of all Sersic index estimates.

    """

    eff = 0.5

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

    EffRad, totmag = GetReff().GetReSer(head, comps, eff, theta)

    # at the moment, this is ignored
    EffRadfrac, totmag = GetReff().GetReSer(head, comps, frac, theta)

    comps2 = copy.deepcopy(comps)
    meanme = GetMe().MeanMe(totmag, EffRad * head.scale)
    me = GetMe().Me(
        head, comps2, EffRad * head.scale, theta
    )  # substituing effective radius per 0 gives m0

    sersic = GetN().MeMeanMe(me, meanme)

    Fa = np.arange(0.1, 0.5, 0.05)
    Fb = np.arange(0.55, 1, 0.05)
    F = np.concatenate((Fa, Fb))  # list containing the fraction of light

    # computes the radius of those fraction of light F
    R = GetReff().GetRfracSer(head, comps, F, theta)

    ns = GetN().GalNs(EffRad, R, F)

    if plot:

        plt.plot(F, ns - const)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Fraction of light")
        plt.ylabel("Sersic index")
        plt.savefig("Serind.png")

    return sersic, np.mean(ns), np.std(ns), totmag, N, theta


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
