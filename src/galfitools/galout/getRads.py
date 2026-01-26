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


def getBreak(
    galfitFile: str,
    dis: float,
    inicomp: int,
    quick: bool,
    random: int,
    num_comp: int,
    angle: float,
    plot: bool,
    ranx: list,
) -> float:
    """Obtains the break radius from a set of Sersic functions.


    Given a model composed of multiple Sersic functions,
    it returns the radius corresponding to the maximum
    of the second derivative.

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    inicomp: int
            Number of component where it'll obtain the initial parameter
            to search break radius or to generated random initial radius

    quick: bool
           If True, it chooses inicomp as a initial parameter

    random: int

          Number of random radii to use as initial parameters for the
          global maximum search. Random radii will be generated within
          the range from 0 to the effective radius of the component
          specified by the inicomp parameter.

    num_comp: int
            Number of component from which the center of all
            components will be determined.

    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components

    plot: bool
          if True, it will makes plot of the second derivative

    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].


    Returns
    -------
    rbreak : float
            break radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    See Also
    --------
    getBreak2 : a more efficient way to compute Break radius


    Notes
    -----
    The maximum of the second derivative involves to find
    the global maximum among many local maximums. Hence the
    choose of the initial parameters is fundamental.


    """

    galfit = Galfit(galfitFile)

    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal = comps.Active == 1

    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute Re")
        print("exiting..")
        sys.exit(1)

    #########################
    #  computing the slope
    #########################

    if plot:  # pragma: no cover

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin, xmax, 0.1)

        gam = GetBreak().GalBreak(R, comps, theta)

        plt.plot(R, gam)

        plt.xlabel("Rad")
        plt.ylabel("Second derivative")

        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("Break.png")

    if quick:  # pragma: no cover

        rbreak = GetBreak().FindBreak(comps, theta, inicomp)

    else:

        if random:  # pragma: no cover

            radius = np.random.random(random) * comps.Rad[maskgal][inicomp]

        else:

            radius = comps.Rad[maskgal]

        rbreak = MulFindBreak(comps, theta, radius)

    return rbreak, N, theta


def MulFindBreak(comps: GalComps, theta: float, radius: list):
    """Given a list of radii, it finds the break radius by
    evaluating each radius in the list as an initial parameter.

    """

    maskgal = comps.Active == 1

    radsbreak = GetBreak().MulFindBreak(comps, theta, radius)

    betas = np.array([])

    for r in radsbreak:

        beta = GetBreak().BreakSer(
            r,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        betas = np.append(betas, beta)

    # in this version, to find the maximum is to find min(betas)
    idx = np.where(betas == min(betas))

    return radsbreak[idx][0]


class GetBreak:
    """Class to obtain the break radius from a set of Sersic functions.
       Class called by getBreak function.


    Methods
    -------
    GalBreak : Evaluates the Break function

    FindBreak : Return the break radii of a set of Sersic functions

    MulFindBreak : Returns the break radius by evaluating it at various
                  initial parameters derived from a list of effective
                  radii.
    funGalBreakSer : evaluates the break function at a given radius

    BreakSer : Break function to a determined R

    """

    """
    def FullSlopeSer(
        self, R: float, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        # not used here

        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta)

        return SlptotR.sum()

    """

    def GalBreak(
        self, R: list, comps: GalComps, theta: float
    ) -> float:  # pragma: no cover

        maskgal = comps.Active == 1  # using active components only

        kappa = np.array([])

        for r in R:

            beta = self.BreakSer(
                r,
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            )

            kappa = np.append(kappa, beta)

        return kappa

    """
    def funGalSlopeSer(
        self,
        R: float,
        Ie: list,
        Re: list,
        n: list,
        q: list,
        pa: list,
        theta: float,
        slope: float,
    ) -> float:
        # not used here

        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta) - slope

        return fun

    """

    def FindBreak(
        self, comps: GalComps, theta: float, initial_comp: int
    ) -> float:  # pragma: no cover

        maskgal = comps.Active == 1  # using active components only

        init = comps.Rad[maskgal][
            initial_comp
        ]  # component as initial value, hope it works

        breakrad = scipy.optimize.fmin(
            self.funGalBreakSer,
            init,
            args=(
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            ),
        )

        return breakrad[0]

    def MulFindBreak(self, comps: GalComps, theta: float, radius: list) -> float:

        brads = np.array([])

        maskgal = comps.Active == 1  # using active components only

        for idx, item in enumerate(radius):

            init = item

            breakrad = scipy.optimize.fmin(
                self.funGalBreakSer,
                init,
                args=(
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                ),
            )

            line = "Optimized radius: {:.2f} \n ".format(breakrad[0])
            print(line)

            brads = np.append(brads, breakrad[0])

        return brads

    def funGalBreakSer(self, R, Ie, Re, n, q, pa, theta):

        beta = self.BreakSer(R, Ie, Re, n, q, pa, theta)

        return beta

    """
    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"
        # not used here

        maskgal = comps.Active == 1  # using active components only

        a = 0.1
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        try:
            Radslp = bisect(
                self.funGalSlopeSer,
                a,
                b,
                args=(
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                    slope,
                ),
            )

        except Exception:
            print("unable to solve equation in the given range. Setting Radslp to 0")
            Radslp = 0

        return Radslp

    """

    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varx = k * ((R / Re) ** (1 / n) - 1)

        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varxpr = (k / n) * ((R / Re) ** (1 / n - 1)) * (1 / Re)

        return varxpr

    def var_Xprim2(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varxpr2 = (k / n) * ((R / Re) ** (1 / n - 2)) * (1 / Re**2) * (1 / n - 1)

        return varxpr2

    def var_S(self, R: float, Ie: list, Re: list, n: list, X: list):

        S = Ie * np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list):

        Sprim = -Ie * np.exp(-X) * Xprim

        return Sprim.sum()

    def var_Sprim2(
        self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list, Xprim2: list
    ):

        Sprim2 = Ie * np.exp(-X) * Xprim**2 - Ie * np.exp(-X) * Xprim2

        return Sprim2.sum()

    def BreakSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)
        Xprim2 = self.var_Xprim2(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)
        Sprim2 = self.var_Sprim2(Rcor, Ie, Re, n, X, Xprim, Xprim2)

        Break = (Sprim / S) * (R / np.log10(np.e)) + (
            (Sprim2 * S - Sprim**2) / S**2
        ) * (R**2 / np.log10(np.e))

        return Break

    def SlopeSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:  # pragma: no cover
        """slope from sersic function to a determined R"""

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)

        Slp = (Sprim / S) * R

        return Slp


def getBreak2(
    galfitFile: str, dis: float, angle: float, num_comp: int, plot: bool, ranx: list
) -> float:
    """Obtains the break radius from a set of Sersic functions.
        This is an alternative method to getBreak


    Given a model composed of multiple Sersic functions,
    it returns the radius corresponding to the maximum
    of the second derivative. This method is more efficient
    to find the global maximum than getBreak

    Parameters
    ----------
    galfitFile : str
            name of the GALFIT file
    dis : float
        maximum distance among components
    num_comp : int
            Number of component from which the center of all
            components will be determined.
    angle : float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    plot : bool
          if True, it will makes plot of the second derivative
    ranx : list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].


    Returns
    -------
    rbreak : float
            break radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    See Also
    --------
    getBreak : an alternative to find Break radius


    Notes
    -----
    The maximum of the second derivative involves to find
    the global maximum among many local maximums. Hence the
    choose of the appropiate range with the `ranx` variable


    """

    galfit = Galfit(galfitFile)

    # head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    maskgal = comps.Active == 1
    if angle:
        theta = angle
    else:
        # taking the last component position angle for the whole galaxy
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute break radius")
        print("exiting..")
        sys.exit(1)

    #########################
    #  computing the slope
    #########################
    if ranx:
        (xmin, xmax) = ranx[0], ranx[1]
    else:
        xmin = 0.1
        xmax = 100

    R = np.arange(xmin, xmax, 0.1)

    lrad = np.log10(R)
    slp = GetSlope().GalSlope(R, comps, theta)

    if plot:  # pragma: no cover
        plt.close()
        plt.plot(R, slp)

        plt.xlabel("Rad")
        plt.ylabel("Slope")

        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("FirstDerivative.png")

    ###
    yspl = UnivariateSpline(lrad, slp, s=0, k=4)

    # yspl2d = yspl.derivative(n=2)
    yspl1d = yspl.derivative(n=1)

    if plot:  # pragma: no cover
        plt.close()
        plt.plot(10 ** (lrad), yspl1d(lrad))
        plt.xlabel("Rad (pixels)")
        plt.ylabel("Second derivative")

        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("SecondDerivative.png")

    idx = np.where(yspl1d(lrad) == max(yspl1d(lrad)))[0][0]

    brad = lrad[idx]

    rbreak = 10**brad

    return rbreak, N, theta


def getKappa2(
    galfitFile: str, dis: int, angle: float, num_comp: int, plot: bool, ranx: list
) -> float:
    """Obtains the kappa radius from a set of Sersic functions.
        This is an alternative method to getKappa function.


    Given a model composed of multiple Sersic functions,
    it returns the radius corresponding to the maximum
    curvature. This method is more efficient
    to find the global maximum than getKappa

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    num_comp: int
            Number of component from which the center of all
            components will be determined.
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    plot: bool
          if True, it will makes plot of the second derivative
    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].


    Returns
    -------
    rbreak : float
            break radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    See Also
    --------
    getKappa: an alternative to find Kappa radius


    Notes
    -----
    The maximum of the curvature involves to find
    the global maximum among many local maximums. Hence the
    choose of the appropiate range with the `ranx` variable


    """

    galfit = Galfit(galfitFile)

    # head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    maskgal = comps.Active == 1
    if angle:
        theta = angle
    else:
        # taking the last component position angle for the whole galaxy
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute kappa radius")
        print("exiting..")
        sys.exit(1)

    #########################
    #  computing the slope
    #########################
    if ranx:
        (xmin, xmax) = ranx[0], ranx[1]
    else:
        xmin = 0.1
        xmax = 100

    R = np.arange(xmin, xmax, 0.1)

    lrad = np.log10(R)
    slp = GetSlope().GalSlope(R, comps, theta)

    ###
    yspl = UnivariateSpline(lrad, slp, s=0, k=4)

    yspl1d = yspl.derivative(n=1)  # this is the second derivative of the Sersics
    # yspl2d = yspl.derivative(n=2)  # this is the Third derivative of the Sersics

    kap1d = np.abs(yspl1d(lrad)) / (1 + slp**2) ** (3 / 2)

    if plot:  # pragma: no cover
        plt.close()
        plt.plot(10 ** (lrad), kap1d)
        plt.xlabel("Rad (pixels)")
        plt.ylabel("Curvature radius")

        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("kappaRadius.png")

    idx = np.where(kap1d == max(kap1d))[0][0]

    krad = lrad[idx]

    rkappa = 10**krad

    return rkappa, N, theta


############################
############################
############################


def getFWHM(galfitFile: str, dis: int, angle: float, num_comp: int):
    """Obtains  gets the FWHM from a set of Sersics functions.

    Given a model composed of multiple Sersic functions,
    it returns the radius with Full Width Half Maximum (FWHM)

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    num_comp: int
            Number of component from which the center of all
            components will be determined.

    Returns
    -------

    fwhm : float
           FWHM radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    """

    galfit = Galfit(galfitFile)

    # head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal = comps.Active == 1

    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute FWHM")
        print("exiting..")
        sys.exit(1)

    #########################
    #  computing the slope
    #########################

    fwhm = GetFWHM().FindFWHM(comps, theta)

    return fwhm, N, theta


class GetFWHM:
    """Class to obtain the FWHM for the whole model
       Class called by getFWHM function.


    Methods
    -------
    FindFWHM : return the fwhm of a set of Sersic functions.
                It uses Bisection

    FullFWHMSer : Surface brightness I(R) from sersic function evaluated at R

    GalFWHM : Evaluates the surface britghtness from a list of R

    funGalFWHMSer : Evaluates Surface brightness at R

    FWHMSer :  Surface brightness I(R) from sersic function evaluated at R


    """

    def FullFWHMSer(
        self, R: float, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:  # pragma: no cover

        SlptotR = self.FWHMSer(R, Re, n, q, pa, theta)

        return SlptotR.sum()

    def GalFWHM(
        self, R: list, comps: GalComps, theta: float
    ) -> float:  # pragma: no cover

        maskgal = comps.Active == 1  # using active components only

        gam = np.array([])

        for r in R:

            slp = self.FWHM(
                r,
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            )
            gam = np.append(gam, slp)

        return gam

    def funGalFWHMSer(
        self,
        R: float,
        S0: float,
        Ie: list,
        Re: list,
        n: list,
        q: list,
        pa: list,
        theta: float,
    ) -> float:

        fun = self.FWHMSer(R, Ie, Re, n, q, pa, theta) - S0 / 2

        return fun

    def FindFWHM(self, comps: GalComps, theta: float) -> float:

        maskgal = comps.Active == 1  # using active components only

        # find the surface brightness at S(R=0)
        X = self.var_X(0, comps.Rad[maskgal], comps.Exp[maskgal])
        S0 = self.var_S(0, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], X)

        a = 0.001
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        try:
            Radslp = bisect(
                self.funGalFWHMSer,
                a,
                b,
                args=(
                    S0,
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                ),
            )

        except Exception:  # pragma: no cover
            print("solution not found in given range")
            Radslp = 0

        return 2 * Radslp

    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varx = k * ((R / Re) ** (1 / n) - 1)

        return varx

    def var_Xprim(self, R: float, Re: list, n: list):  # pragma: no cover

        k = gammaincinv(2 * n, 0.5)

        varxpr = (k / n) * ((R / Re) ** (1 / n - 1)) * (1 / Re)

        return varxpr

    def var_S(self, R: float, Ie: list, Re: list, n: list, X: list):

        S = Ie * np.exp(-X)

        return S.sum()

    def var_Sprim(
        self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list
    ):  # pragma: no cover

        Sprim = Ie * np.exp(-X) * Xprim

        return Sprim.sum()

    def FWHMSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """I(R) from sersic function to a determined R"""

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        # Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        # Sprim = self.var_Sprim(Rcor, Ie,  Re, n, X, Xprim)

        Ifwhm = S

        return Ifwhm


############################
############################
############################


def getKappa(
    galfitFile: str,
    dis: int,
    inicomp: int,
    quick: bool,
    random: int,
    angle: float,
    num_comp: int,
    plot: bool,
    ranx: list,
) -> float:
    """Obtains the Kappa radius from a set of Sersic functions.

    Given a model composed of multiple Sersic functions,
    it returns the radius corresponding to the maximum
    curvature.

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    inicomp: int
            Number of component where it'll obtain the initial parameter
            to search kappa radius or to generated random initial radius

    quick: bool
           If True, it chooses inicomp as a initial parameter

    random: int

          Number of random radii to use as initial parameters for the
          global maximum search. Random radii will be generated within
          the range from 0 to the effective radius of the component
          specified by the inicomp parameter.

    num_comp: int
            Number of component from which the center of all
            components will be determined.

    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components

    plot: bool
          if True, it will makes plot of the second derivative

    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].


    Returns
    -------
    rkappa : float
            kappa radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    See Also
    --------
    getKappa2 : a more efficient way to compute Kappa radius


    Notes
    -----
    The maximum of the curvature involves to find
    the global maximum among many local maximums. Hence the
    choose of the initial parameters is fundamental.


    """

    galfit = Galfit(galfitFile)

    # head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal = comps.Active == 1
    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute Re")
        print("exiting..")
        sys.exit(1)

    line = "Using a theta value of : {:.2f} degrees\n".format(theta)
    print(line)

    #########################
    #  computing the slope
    #########################

    if plot:  # pragma: no cover

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin, xmax, 0.1)

        gam = GetKappa().GalKappa(R, comps, theta)

        plt.plot(R, gam)
        plt.xlabel("Radius")
        plt.ylabel("Curvature")
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("Kappa.png")

    if quick:  # pragma: no cover

        rkappa = GetKappa().FindKappa(comps, theta, inicomp)

    else:  # pragma: no cover

        if random:

            radius = np.random.random(random) * comps.Rad[maskgal][inicomp]

        else:

            radius = comps.Rad[maskgal]

        rkappa = MultiFindKappa(comps, theta, radius)

    return rkappa, N, theta


def MultiFindKappa(comps, theta, radius):
    """Given a list of radii, it finds the kappa radius by
    evaluating each radius in the list as an initial parameter.

    """

    maskgal = comps.Active == 1

    radskappa = GetKappa().MulFindKappa(comps, theta, radius)

    kappas = np.array([])

    print("finding global minimum:")

    for r in radskappa:

        kappa = GetKappa().funGalKappaSer(
            r,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        kappas = np.append(kappas, kappa)

    # radmask = radskappa < radius[-1]  # hope it works

    idx = np.where(kappas == min(kappas))

    return radskappa[idx][0]


class GetKappa:
    """Class to obtain the kappa radius from a set of Sersic functions.
       Class called by getKappa function.


    Methods
    -------
    GalKappa : Evaluates the kappa function

    FindKappa : Return the kappa radii of a set of Sersic functions

    MulFindKappa: Returns the kappa radius by evaluating it at various
                  initial parameters derived from a list of effective
                  radii.
    funGalKappaSer : evaluates the kappa function at a given radius

    BetaSer : kappa function to a determined R


    """

    """
    def FullSlopeSer(
        self, R: float, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        # not used here

        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta)

        return SlptotR.sum()

    """

    def GalKappa(
        self, R: list, comps: GalComps, theta: float
    ) -> float:  # pragma: no cover

        maskgal = comps.Active == 1  # using active components only

        kappa = np.array([])

        for r in R:

            beta = self.BetaSer(
                r,
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            )

            gam = self.SlopeSer(
                r,
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            )

            kap = np.abs(beta) / (1 + gam**2) ** (3 / 2)

            kappa = np.append(kappa, kap)

        return kappa

    def funGalSlopeSer(
        self,
        R: float,
        Ie: list,
        Re: list,
        n: list,
        q: list,
        pa: list,
        theta: float,
        slope: float,
    ) -> float:  # pragma: no cover

        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta) - slope

        return fun

    def FindKappa(
        self, comps: GalComps, theta: float, initial_comp: int
    ) -> float:  # pragma: no cover
        "return the break radius of a set of Sersic functions"

        maskgal = comps.Active == 1  # using active components only

        init = comps.Rad[maskgal][
            initial_comp
        ]  # component as initial value, hope it works

        kapparad = scipy.optimize.fmin(
            self.funGalKappaSer,
            init,
            args=(
                comps.Ie[maskgal],
                comps.Rad[maskgal],
                comps.Exp[maskgal],
                comps.AxRat[maskgal],
                comps.PosAng[maskgal],
                theta,
            ),
        )

        return kapparad[0]

    def funGalKappaSer(self, R, Ie, Re, n, q, pa, theta):

        beta = self.BetaSer(R, Ie, Re, n, q, pa, theta)

        gam = self.SlopeSer(R, Ie, Re, n, q, pa, theta)

        kappa = np.abs(beta) / (1 + gam**2) ** (3 / 2)

        return -kappa

    def MulFindKappa(self, comps: GalComps, theta: float, radius: list) -> float:
        "return the kappa radius evaluated at different effective radius"

        krads = np.array([])

        maskgal = comps.Active == 1  # using only active components

        for idx, item in enumerate(radius):

            init = item

            kapparad = scipy.optimize.fmin(
                self.funGalKappaSer,
                init,
                args=(
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                ),
            )

            line = "Optimized radius: {:.2f} \n ".format(kapparad[0])
            print(line)

            krads = np.append(krads, kapparad[0])

        return krads

    """
    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"
        # not used here

        maskgal = comps.Active == 1  # using active components only

        a = 0.1
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        try:
            Radslp = bisect(
                self.funGalSlopeSer,
                a,
                b,
                args=(
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                    slope,
                ),
            )

        except Exception:
            print("solution not found in given range")
            Radslp = 0

        return Radslp

    """

    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varx = k * ((R / Re) ** (1 / n) - 1)

        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varxpr = (k / n) * ((R / Re) ** (1 / n - 1)) * (1 / Re)

        return varxpr

    def var_Xprim2(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varxpr2 = (k / n) * ((R / Re) ** (1 / n - 2)) * (1 / Re**2) * (1 / n - 1)

        return varxpr2

    def var_S(self, R: float, Ie: list, Re: list, n: list, X: list):

        S = Ie * np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list):

        Sprim = -Ie * np.exp(-X) * Xprim

        return Sprim.sum()

    def var_Sprim2(
        self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list, Xprim2: list
    ):

        Sprim2 = Ie * np.exp(-X) * Xprim**2 - Ie * np.exp(-X) * Xprim2

        return Sprim2.sum()

    def BetaSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """Kappa from sersic function to a determined R"""

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)
        Xprim2 = self.var_Xprim2(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)
        Sprim2 = self.var_Sprim2(Rcor, Ie, Re, n, X, Xprim, Xprim2)

        Beta = (Sprim / S) * (R / np.log10(np.e)) + (
            (Sprim2 * S - Sprim**2) / S**2
        ) * (R**2 / np.log10(np.e))

        return Beta

    def SlopeSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """slope from sersic function to a determined R"""

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)

        Slp = (Sprim / S) * R

        return Slp


############################
############################
############################


def getReComp(
    galfitFile: str, dis: int, eff: float, angle: float, num_comp: int
) -> float:
    """Obtains the effective radius from a set of Sersic functions.

    Given a model composed of multiple Sersic functions,
    it returns the effective radius or any other radius
    that represents a percent of total light.

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    eff: float
         fraction of total light. The function will find
         the radius at this value. eff must be between 0 and 1
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    num_comp: int
            Number of component from which the center of all
            components will be determined.

    Returns
    -------

    EffRad : float
            radius found at `eff` of total light in pixels
    totmag : float
             total magnitude
    meanme : float
             mean of surface brightness at effective radius
             in mag/arcsec**2
    me : float
             surface brightness at EffRad  in mag/arcsec**2
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    Notes
    -----
    The function call a class to find EffRad using the
    Bisection method.


    """

    assert (eff > 0) and (eff <= 1), "effrad must be a value between 0 and 1"

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

    EffRad, totmag = GetReff().GetReSer(head, comps, eff, theta)

    meanme = GetMe().MeanMe(totmag, EffRad * head.scale)
    me = GetMe().Me(head, comps, EffRad * head.scale, theta)

    # EffRad_arc = EffRad*head.scale

    return EffRad, totmag, meanme, me, N, theta


class GetMe:
    """Class to obtain the surface brightness at effective radius
       and the mean surface brightness at effective radius
       in units of mag/''

    Methods
    -------
    MeanMe : Mean Surface brightness I(R) at effective radius

    Me : Surface brightness I(R) at R

    Itotser : Evaluates the Sersic surface brightness flux
             from a list of R

    Iser : Sersic surface brightness flux at R


    """

    def MeanMe(self, magtot: float, effrad: float) -> float:

        meanme = magtot + 2.5 * np.log10(2 * np.pi * effrad**2)

        return meanme

    def Me(self, head, comps, EffRad, theta):

        comps.Rad = comps.Rad * head.scale  # converting to arcsec
        comps.Flux = 10 ** ((-comps.Mag) / 2.5)

        k = gammaincinv(2 * comps.Exp, 0.5)

        denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
        denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
        # denom3 = (gamma(2*comps.Exp))*(comps.AxRat)
        denom3 = gamma(2 * comps.Exp)

        denom = denom1 * denom2 * denom3

        comps.Ie = comps.Flux / denom

        maskgal = comps.Active == 1
        Itotr = self.Itotser(
            EffRad,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        me = -2.5 * np.log10(Itotr)

        return me

    def Itotser(
        self, R: float, Ie: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        ItotR = self.Iser(R, Ie, rad, n, q, pa, theta)

        return ItotR.sum()

    def Iser(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        k = gammaincinv(2 * n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta)

        Ir = Ie * np.exp(-k * ((Rcor / Re) ** (1 / n) - 1))

        return Ir


class GetReff:
    """Class to obtain the effective radius of the composed
        surface brightness model

    class is called by GetReComp function. Main method is
    GetReSer.


    Attributes
    ----------
    Attributes of the main method GetReser.
    galhead: GalHead data class defined in galin/galfit.py
             data class that stores header information of
             Galfit file.
    comps: GalComps data class defined in galin/galfit.py
            data class that stores the components parameters
            of galfit file.
    eff: float
         fraction of total light. The method will find
         the radius at this value. eff must be between 0 and 1
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees


    Methods
    -------
    GetReSer : computes the effective radius from the
                composed Sersic model using the bisection method


    """

    def GetRfracSer(self, head, comps, F, theta):

        rads = np.array([])

        for f in F:

            r, totmag = GetReff().GetReSer(head, comps, f, theta)

            rads = np.append(rads, r)

        return rads

    def GetReSer(
        self, galhead: GalHead, comps: GalComps, eff: float, theta: float
    ) -> float:

        maskgal = comps.Active == 1

        comps.Flux = 10 ** ((galhead.mgzpt - comps.Mag) / 2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5 * np.log10(totFlux) + galhead.mgzpt

        a = 0.1
        b = comps.Rad[maskgal][-1] * 1000  # hope it doesn't crash

        Reff = self.solveSerRe(
            a,
            b,
            comps.Flux[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            totFlux,
            eff,
            theta,
        )

        return Reff, totmag

    def solveSerRe(
        self,
        a: float,
        b: float,
        flux: list,
        rad: list,
        n: list,
        q: list,
        pa: list,
        totFlux: float,
        eff: float,
        theta: float,
    ) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"

        try:
            Re = bisect(
                self.funReSer, a, b, args=(flux, rad, n, q, pa, totFlux, eff, theta)
            )

        except Exception:  # pragma: no cover
            print("solution not found in given range")
            Re = 0

        return Re

    def funReSer(
        self,
        R: float,
        flux: list,
        rad: list,
        n: list,
        q: list,
        pa: list,
        totFlux: float,
        eff: float,
        theta: float,
    ) -> float:

        fun = self.Ftotser(R, flux, rad, n, q, pa, theta) - totFlux * eff

        return fun

    def Ftotser(
        self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """partial sersic flux computed from zero up to R for a set of sersics"""

        ftotR = self.Fser(R, flux, rad, n, q, pa, theta)

        return ftotR.sum()

    def Fser(
        self, R: float, Flux: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """partial sersic flux computed from zero up to R for a single sersic"""

        k = gammaincinv(2 * n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta)

        X = k * (Rcor / Re) ** (1 / n)

        Fr = Flux * gammainc(2 * n, X)

        return Fr


def getSlope(
    galfitFile: str,
    dis: int,
    slope: float,
    angle: float,
    num_comp: int,
    plot: bool,
    ranx: list,
) -> float:
    """Obtains the slope radius from a set of Sersic functions.

    Given a model composed of multiple Sersic functions,
    this function returns the radius corresponding to the
    derivative value specified by the slope variable. It calls
    to class GetSlope to computed using the bisection method


    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    dis: float
        maximum distance among components
    slope : float
        Value of the slope at which the radius will be determined.
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    num_comp: int
            Number of component from which the center of all
            components will be determined.
    plot: bool
          if True, it will makes plot of the second derivative

    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].

    Returns
    -------
    rgam: float
          the slope radius in pixels
    N : int
        total number of components
    theta : float
           Angular position indicating the direction along
           which break radius is computed. In degrees



    """

    galfit = Galfit(galfitFile)

    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    k = gammaincinv(2 * comps.Exp, 0.5)

    denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
    denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
    denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)

    denom = denom1 * denom2 * denom3

    comps.Ie = comps.Flux / denom

    comps = SelectGal(comps, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal = comps.Active == 1

    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute Re")
        print("exiting..")
        sys.exit(1)

    #########################
    # computing the slope
    #########################

    if plot:  # pragma: no cover

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin, xmax, 0.1)

        gam = GetSlope().GalSlope(R, comps, theta)

        plt.plot(R, gam)
        plt.xlabel("Radius")
        plt.ylabel("Slope")
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("slope.png")

    rgam = GetSlope().FindSlope(comps, theta, slope)

    return rgam, N, theta


class GetSlope:
    """Class to obtain the slope radius from a set of Sersic functions.

    This class is called by GetSlope Function to compute
    the radius at the given slope value.
    The main method is FindSlope.

    Methods
    -------
    GalSlope: Evaluates the first derivative at the specified radius (R).

    FindSlope: Return the slope radius of a set of Sersic functions.
               It uses bisection to find the radius at the specified slope


    """

    def FullSlopeSer(
        self, R: float, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:  # pragma: no cover

        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta)

        return SlptotR.sum()

    def GalSlope(self, R: list, comps: GalComps, theta: float) -> float:

        maskgal = comps.Active == 1  # using active components only

        gam = np.array([])

        if comps.NameComp[maskgal][1] == "ferrer":

            masksersic = comps.NameComp[maskgal] == "sersic"
            maskferrer = comps.NameComp[maskgal] == "ferrer"

            for r in R:

                slp1 = self.SlopeFerrer(
                    r,
                    comps.Exp2[maskgal][maskferrer],  # change to beta
                    comps.Rad[maskgal][maskferrer],
                    comps.Exp[maskgal][maskferrer],
                    comps.AxRat[maskgal][maskferrer],
                    comps.PosAng[maskgal][maskferrer],
                    theta,
                )

                slp2 = self.SlopeSer(
                    r,
                    comps.Ie[maskgal][masksersic],  # change to beta
                    comps.Rad[maskgal][masksersic],
                    comps.Exp[maskgal][masksersic],
                    comps.AxRat[maskgal][masksersic],
                    comps.PosAng[maskgal][masksersic],
                    theta,
                )

                gam = np.append(gam, slp1 + slp2)

        else:

            for r in R:
                slp = self.SlopeSer(
                    r,
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                )
                gam = np.append(gam, slp)

        return gam

    def funGalSlopeSer(
        self,
        R: float,
        Ie: list,
        Re: list,
        n: list,
        q: list,
        pa: list,
        theta: float,
        slope: float,
    ) -> float:

        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta) - slope

        return fun

    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"

        maskgal = comps.Active == 1  # using active components only

        a = 0.1
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        try:
            Radslp = bisect(
                self.funGalSlopeSer,
                a,
                b,
                args=(
                    comps.Ie[maskgal],
                    comps.Rad[maskgal],
                    comps.Exp[maskgal],
                    comps.AxRat[maskgal],
                    comps.PosAng[maskgal],
                    theta,
                    slope,
                ),
            )

        except Exception:  # pragma: no cover
            print("solution not found in given range")
            Radslp = 0

        return Radslp

    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varx = k * ((R / Re) ** (1 / n) - 1)

        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2 * n, 0.5)

        varxpr = (k / n) * ((R / Re) ** (1 / n - 1)) * (1 / Re)

        return varxpr

    def var_S(self, R: float, Ie: list, Re: list, n: list, X: list):

        S = Ie * np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list):

        Sprim = Ie * np.exp(-X) * Xprim

        return Sprim.sum()

    def SlopeSer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """slope from sersic function to a determined R"""

        Rcor = GetRadAng(R, q, pa, theta)

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie, Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)

        Slp = (Sprim / S) * R

        return Slp

    def SlopeFerrer(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """slope from Ferrer function to a determined R

        def dlogSigma_dlogr_analytic2(r, Sigma0=1.0, rout=10.0, beta=0.5, alpha=2.0):
        d log10(Σ) / d log10(r) = (r/Σ) dΣ/dr = d ln(Σ) / d ln(r)

        for Σ(r)=Σ0(1-x)^α, x=(r/rout)^(2-β), n=2-β:
        ln Σ = ln Σ0 + α ln(1-x)
        d ln Σ / d ln r = α * [1/(1-x)] * d(1-x)/d ln r
                       = α * [1/(1-x)] * (-(d x / d ln r))
        and  dx/d ln r = n x,
        => d log Σ / d log r = - α * n * x / (1-x)

        Note: parameters here are not Sersic anymore, but
              these use the same positional arguments in
              sersic file. For instance Re here really means
              Outer truncation radius of the Sersic Functions
              Ie does not have an equivalent, but it is used as beta here.
        """

        beta = Ie
        alpha = n
        Rout = Re

        Rcor = GetRadAng(R, q, pa, theta)

        X = Rcor / Rout  # x = r / rout

        gam = 2.0 - beta

        dev = alpha * (beta - 2) * X ** (gam) / (1 - X**gam)

        return dev


def getBulgeRad(galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx):
    """Determines the radius at which the surface brightness of
    two models are equal.


    Given two surface brightness models, this function identifies
    the radius where their surface brightness values are equal. This
    is useful for determining the radius where the surface brightness
    of the bulge and disk are equal, or for finding the radius between
    the core and coreless regions in elliptical galaxies.


    Parameters
    ----------
    galfitFile1 : str
            name of the GALFIT file with model 1
    galfitFile2 : str
            name of the GALFIT file with model 2

    dis: float
        maximum distance among components
    num_comp: int
            Number of component from which the center of all
            components will be determined.
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the last components
    plot: bool
          if True, it will makes plot of the second derivative
    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].

    Returns
    -------
    rbulge: float
            radius where the I1 and I2 are equal in pixels
    N1 : int
        total number of components for model 1
    N2 : int
        total number of components for model 2

    theta : float
           Angular position indicating the direction along
           which this radius is computed. In degrees


    """

    galfit1 = Galfit(galfitFile1)
    head1 = galfit1.ReadHead()
    galcomps1 = galfit1.ReadComps()

    galfit2 = Galfit(galfitFile2)
    galcomps2 = galfit2.ReadComps()

    galcomps1 = SelectGal(galcomps1, dis, num_comp)

    galcomps2 = SelectGal(galcomps2, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    maskgal2 = galcomps2.Active == 1

    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = galcomps2.PosAng[maskgal2][-1]

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps1 = conver2Sersic(galcomps1)
    comps2 = conver2Sersic(galcomps2)

    N1 = numComps(comps1, "all")
    N2 = numComps(comps2, "all")

    if (N1 == 0) or (N2 == 0):  # pragma: no cover
        print("not enough number of components of one of the two models to proceed")
        print("exiting..")
        sys.exit(1)

    if plot:  # pragma: no cover

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin, xmax, 0.1)

        Idiff = getDiff(head1, comps1, comps2, R, theta)

        plt.clf()
        plt.plot(R, Idiff)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Rad")
        plt.ylabel("I1 - I2")
        plt.savefig("Idiff.png")

        I1, I2 = getIs(head1, comps1, comps2, R, theta)

        plt.clf()
        plt.plot(R, I1)
        plt.plot(R, I2)
        plt.xscale("log")
        plt.xlabel("Rad")
        plt.ylabel("I")
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("BulgeRad.png")

    # computing bulge radius
    try:
        rbulge = newton(getDiffx, 0, args=(head1, comps1, comps2, theta))
    except Exception:  # pragma: no cover
        print("solution not found for initial parameter")
        rbulge = 0

    return rbulge, N1, N2, theta


def getDiff(head1, comps1, comps2, R, theta):  # pragma: no cover
    """Calculates the surface brightness difference
    between two models at various radii I1 - I2.

    """
    Idiff = np.array([])

    for r in R:

        Ir1 = GetIr().Ir(head1, comps1, r, theta)
        Ir2 = GetIr().Ir(head1, comps2, r, theta)

        Ird = Ir1 - Ir2
        Idiff = np.append(Idiff, Ird)

    return Idiff


def getIs(head1, comps1, comps2, R, theta):  # pragma: no cover
    """Calculates the surface brightness of two
    models at various radii.

    """

    I1 = np.array([])
    I2 = np.array([])

    for r in R:

        Ir1 = GetIr().Ir(head1, comps1, r, theta)
        Ir2 = GetIr().Ir(head1, comps2, r, theta)

        I1 = np.append(I1, Ir1)
        I2 = np.append(I2, Ir2)

    return I1, I2


def getDiffx(r, head1, comps1, comps2, theta):
    """Calculates the surface brightness difference
    between two models at a specified radii I1 - I2.

    """

    Ir1 = GetIr().Ir(head1, comps1, r, theta)
    Ir2 = GetIr().Ir(head1, comps2, r, theta)

    Irdx = Ir1 - Ir2

    return Irdx


class GetIr:
    """
    Class called by getDiff, getDiffx, getIs
    to obtain the surface brightness from  a set
    of Sersic functions at different radii

    the main method is Ir

    Methods
    -------
    Ir : obtains the surface brightness at R

    Itotser : method called by Ir to obtain the
              total surface brightness of the
              set components at R
    Iser : method called by Itotser to obtain
           the surface brightnness
           per component at R


    """

    def Ir(self, head, comps, R, theta):

        # comps.Rad = comps.Rad*head.scale
        comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

        k = gammaincinv(2 * comps.Exp, 0.5)

        denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
        denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
        # denom3 = (gamma(2*comps.Exp))*(comps.AxRat)
        denom3 = gamma(2 * comps.Exp)

        denom = denom1 * denom2 * denom3

        comps.Ie = comps.Flux / denom

        maskgal = comps.Active == 1

        Itotr = self.Itotser(
            R,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        return Itotr

    def Itotser(
        self, R: float, Ie: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        ItotR = self.Iser(R, Ie, rad, n, q, pa, theta)

        return ItotR.sum()

    def Iser(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """sersic flux to a determined R"""

        k = gammaincinv(2 * n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta)

        Ir = Ie * np.exp(-k * ((Rcor / Re) ** (1 / n) - 1))

        return Ir
