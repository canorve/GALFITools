#! /usr/bin/env python

import sys

import matplotlib.pyplot as plt
import numpy as np
from galfitools.galin.galfit import (
    GalComps,
    Galfit,
    GalHead,
    GetRadAng,
    SelectGal,
    conver2Sersic,
    numComps,
)

from galfitools.galout.getRads import GetReff
from scipy.optimize import bisect
from scipy.special import gamma, gammainc, gammaincinv


def getCOW(
    galfitFile: str,
    dis: int,
    angle: float,
    frac: float,
    num_comp: int,
    plotname: str,
    dpival: int,
    galfitF2: str,
    maxdiff: bool,
) -> float:
    """plots the curve-of-growth

    Makes a plot of the Curver-of-Growth from the galfit.XX file.
    The curve model is made from the Sersic functions.

    Parameters
    ----------
    galfitFile: str
                name of the GALFIT file
    dis: int
        maximum distance among components
    angle: float
            angular position of the major axis galaxy.
            by default (if it is set to None)
            it will take the last component
            of the model in the GALFIT file.
    frac: float
          fraction of light radius. This is the upper limit of
          X-Axis.
    num_comp: int,
             Number of component where it'll obtain the center (X,Y)
    plotname: str
            name of the output plot fileh
    dpival: int
            value of the dpi (dots per inch) for the plot
    galfitF2: str
            Second GALFIT file to add to the plot (optional).
    maxdiff: bool
            plot the maximum difference as a vertical line between
            model 1 and 2 (galfitF2)

    Returns
    -------

    totmag : float
            total magnitude
    N : int
        total number of components of the galaxy
    theta : float
           Angular position indicating the direction of
           the galaxy's curve of growth.

    Warnings
    --------
    use `dis` parameter with precaution. The equations assume
    that the components of the same galaxy share the same center.
    greater values of dis will produce wrong computations of COW


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
    # print('number of model components: ',N)

    if N == 0:  # pragma: no cover
        print("not enough number of components to plot COW")
        print("exiting..")
        sys.exit(1)

    # line = "Using a theta value of : {:.2f} degrees \n".format(theta)
    # print(line)

    EffRadfrac, totmag = GetReff().GetReSer(
        head, comps, frac, theta
    )  # obtains the radius at 95% of total flux

    # computing the Sersic indexes for different radius

    R = np.arange(0.1, EffRadfrac, 0.1)

    cows, totmag = GetCOW().GetCOWrad(head, comps, theta, R)

    if galfitF2:  # pragma: no cover

        # it uses the same num_comp as the first file
        head2, comps2, theta2 = readGalfitF2(galfitF2, dis, angle, num_comp)

        cows2, totmag2 = GetCOW().GetCOWrad(
            head2, comps2, theta, R
        )  # same theta as galfitF1

        diff = np.abs(cows - cows2)

        idx = np.where(diff == np.max(diff))

        xline = R[idx][0]

        ymin = cows[idx]
        ymax = cows2[idx]

    plt.plot(R, cows, color="blue", label="model 1")
    plt.xlabel("Rad")
    plt.ylabel("curve of growth")
    plt.grid(True)
    plt.minorticks_on()

    xmin = 0.1
    xmax = np.max(R)

    plt.gca().invert_yaxis()

    plt.xlabel("Radius (pixels)")
    plt.title("Curve of Growth ")
    plt.ylabel("magnitude (< R) ")

    if galfitF2:  # pragma: no cover

        plt.plot(R, cows2, color="green", label="model 2")

        plt.xlabel("Rad")
        plt.ylabel("curve of growth")

        if maxdiff:
            plt.vlines(xline, ymin, ymax, color="red", label="max diff")

    plt.hlines(totmag, xmin, xmax, color="black", label="total magnitude")

    plt.legend(loc="lower right")
    plt.savefig(plotname, dpi=dpival)

    return totmag, N, theta


def readGalfitF2(galfitF2, dis, angle, num_comp):
    """
    Reads the GALFIT file to obtain header and
    components information

    """

    galfit = Galfit(galfitF2)

    head2 = galfit.ReadHead()
    galcomps2 = galfit.ReadComps()

    galcomps2 = SelectGal(galcomps2, dis, num_comp)

    # taking the last component position angle for the whole galaxy
    maskgal = galcomps2.Active == 1
    if angle:
        theta2 = angle
    else:
        theta2 = galcomps2.PosAng[maskgal][-1]

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps2 = conver2Sersic(galcomps2)

    N = numComps(comps2, "all")

    if N == 0:
        print("not enough number of components to plot COW")
        print("exiting..")
        sys.exit(1)

    return head2, comps2, theta2


class GetCOW:
    """Obtains the Curve-of-Growth at a given radius

    The main method is GetCOWrad which calls to the other methods
    This is the one to be used. See the parameters of this method
    below:

    Parameters
    ----------
    head :  GalHead data class defined in galin/galfit.py
    comps : GalComps data class defined in galin/galfit.py
    theta : float
           Angular position indicating the direction of
           the galaxy's curve of growth.
    R : array
        Array indicating the radius at each point where the
        curve of growth will be computed.

    Methods
    -------
    GetCOWrad : given a array of R, obtains the curve of growth

    GetCOWFluxR : selects the components of the galaxy and gets the
                total flux of the COW at a given R

    funCOWSer : obtains the function of the COW

    COWtotser : obtains the total Sersic flux of all components
                to a determined R

    COWser : obtains the Sersic flux of one component to a determined R


    """

    def GetCOWrad(self, head, comps, theta, R):

        cowfluxs = np.array([])

        for r in R:  # check if functions works without the for

            cowflux = self.GetCOWFluxR(head, comps, theta, r)

            cowfluxs = np.append(cowfluxs, cowflux)

        maskgal = comps.Active == 1

        comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5 * np.log10(totFlux) + head.mgzpt

        cows = -2.5 * np.log10(cowfluxs) + head.mgzpt

        return cows, totmag

    def GetCOWFluxR(
        self, galhead: GalHead, comps: GalComps, theta: float, rad: float
    ) -> float:

        maskgal = comps.Active == 1

        comps.Flux = 10 ** ((galhead.mgzpt - comps.Mag) / 2.5)

        # rad refers to the radius where flux will be computed
        # comps.Rad refers to the effective radius of every component

        cowR = self.funCOWSer(
            rad,
            comps.Flux[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        return cowR

    def funCOWSer(
        self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        fun = self.COWtotser(R, flux, rad, n, q, pa, theta)  # - totFlux*eff

        return fun

    def COWtotser(
        self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        ftotR = self.COWser(R, flux, rad, n, q, pa, theta)

        return ftotR.sum()

    def COWser(
        self, R: float, Flux: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """sersic flux to a determined R"""

        k = gammaincinv(2 * n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta)

        X = k * (Rcor / Re) ** (1 / n)

        Fr = Flux * gammainc(2 * n, X)

        return Fr


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
