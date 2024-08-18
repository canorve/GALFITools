#! /usr/bin/env python

import sys

import numpy as np
from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    conver2Sersic,
    numComps,
)
from scipy.special import gamma, gammainc, gammaincinv


def getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad):
    """gets the missing light from the difference between two GALFIT models

    From two surface brightness models of the same galaxy, this function
    computes the magnitude of the difference between the two models. The
    calculation extends from the center to the radius provided by the user

    Parameters
    ----------
    galfitFile1 : str
                Galfit File containing the coreless
                surface brightness model
    galfitFile2 : str
                Galfit File containing the core surface
                brightness model
    dis : float
        Maximum distance among components
    num_comp : int
            number of component to select center of galaxy. It will
            apply for both GALFIT files
    rad : float
          upper limit of radius to integrate the missing light in pixels

    Returns
    -------
    magmiss : float
            magnitude of the missing light

    N1 : number of components of model 1
    N2 : number of components of model 2

    Warnings
    --------
    It only works for Sersic functions and their derivations
    such as de Vaucouleurs, exponential, gaussian.

    """

    galfit1 = Galfit(galfitFile1)
    head1 = galfit1.ReadHead()
    galcomps1 = galfit1.ReadComps()

    galfit2 = Galfit(galfitFile2)
    head2 = galfit2.ReadHead()
    galcomps2 = galfit2.ReadComps()

    galcomps1 = SelectGal(galcomps1, dis, num_comp)
    galcomps2 = SelectGal(galcomps2, dis, num_comp)

    # taking the last component position angle for the whole galaxy

    # maskgal1 = galcomps1.Active == 1
    # maskgal2 = galcomps2.Active == 1

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps1 = conver2Sersic(galcomps1)
    comps2 = conver2Sersic(galcomps2)

    N1 = numComps(comps1, "all")
    N2 = numComps(comps2, "all")

    if (N1 == 0) or (N2 == 0):
        print("not enough number of components of one of the two models to proceed")
        print("exiting..")
        sys.exit(1)

    Flux1, mag1 = GetMag().GetFlux(head1, comps1, rad)
    Flux2, mag2 = GetMag().GetFlux(head2, comps2, rad)

    magmiss = getMiss(head1, mag1, mag2)

    return magmiss, N1, N2


class GetMag:
    """Obtains the total flux and magnitude up
    to a given radius, applicable only for Sersic functions.


    """

    def GetFlux(self, head, comps, gamRad):

        # comps.Rad = comps.Rad*head.scale
        comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

        k = gammaincinv(2 * comps.Exp, 0.5)

        denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
        denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
        # axis ratio must be considered if mag total will be used:
        denom3 = (gamma(2 * comps.Exp)) * (comps.AxRat)
        # denom3 = (gamma(2*comps.Exp))

        denom = denom1 * denom2 * denom3

        comps.Ie = comps.Flux / denom

        maskgal = comps.Active == 1

        Ftotr = self.Ftotser(
            gamRad,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
        )

        mag = -2.5 * np.log10(Ftotr) + head.mgzpt

        return Ftotr, mag

    def Ftotser(
        self, R: float, Ie: list, rad: list, n: list, q: list, pa: list
    ) -> float:

        FtotR = self.Fser(R, Ie, rad, n, q, pa)

        return FtotR.sum()

    def Fser(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list) -> float:
        """sersic flux to a determined R"""

        k = gammaincinv(2 * n, 0.5)

        # Rcor = GetRadAng(R, q, pa)
        x = k * (R / Re) ** (1 / n)

        Fr = (
            Ie
            * (Re**2)
            * (2 * np.pi * q * n)
            * np.exp(k)
            * (k ** (-2 * n))
            * gammainc(2 * n, x)
            * gamma(2 * n)
        )

        return Fr


def getMiss(head, mag1, mag2):
    """obtains the magnitude difference

    Parameters
    ----------
    head : GalHead data class defined in galfit.py
    mag1 : Magnitude 1
    mag2 : Magnitude 2

    Returns
    -------
    magMiss : float
            magnitude difference

    """

    Flux1 = 10 ** ((head.mgzpt - mag1) / 2.5)
    Flux2 = 10 ** ((head.mgzpt - mag2) / 2.5)

    if Flux1 > Flux2:

        Fluxmiss = Flux1 - Flux2
    else:
        print("Warning: Flux1 is less than Flux2 for the selected radius. Inverting")

        Fluxmiss = Flux2 - Flux1

    magMiss = -2.5 * np.log10(Fluxmiss) + head.mgzpt

    return magMiss


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
