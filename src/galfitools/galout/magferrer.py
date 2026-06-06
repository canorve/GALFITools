#! /usr/bin/env python3

import numpy as np
import argparse

from galfitools.shell.prt import printWelcome
from scipy.special import beta, betainc
from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
)


def ferrer_luminosity(R, Sigma0, r_out, alpha, beta_par):
    """
    Luminosity inside radius R for a Ferrer-like profile.

    Parameters
    ----------
    R : float
        Radius where the luminosity is evaluated.
    Sigma0 : float
        Central surface brightness in linear flux units.
    r_out : float
        Outer truncation radius.
    alpha : float
        Outer truncation shape parameter.
    beta_par : float
        Inner slope parameter.

    Returns
    -------
    L : float
        Luminosity inside radius R.
    """
    gamma = 2.0 - beta_par

    if gamma <= 0:
        raise ValueError("The parameter 2 - beta must be positive.")

    R_eff = min(R, r_out)
    x = (R_eff / r_out) ** gamma

    a = 2.0 / gamma
    b = alpha + 1.0

    # scipy.special.betainc gives the regularized incomplete beta function.
    incomplete_beta = betainc(a, b, x) * beta(a, b)

    L = (2.0 * np.pi * Sigma0 * r_out**2 / gamma) * incomplete_beta

    return L


def ferrers_magnitude(
    R: float,
    mu0_mag_arcsec2: float,
    r_out: float,
    alpha: float,
    beta_par: float,
    pixscale: float,
) -> float:
    """
    Compute the magnitude of the GALFIT ferrers model
    inside radius R.

    Parameters
    ----------
    R: float
         radius in pixels.
    mu0_mag_arcsec2 : float
        Central surface brightness in mag/arcsec^2.
    r_out: float
        outer truncation radius in pixels.
    alpha: float
         outer truncation shapness
    beta_par: float
        central slope.
    pixscale: float
        Pixel scale in arcsec/pixel.


    Returns
    -------
    float
        magnitude inside radius R

    """

    r_out_arcsec = r_out * pixscale
    R_arcsec = R * pixscale

    Sigma0 = 10 ** (-mu0_mag_arcsec2 / 0.4)

    lum = ferrer_luminosity(R_arcsec, Sigma0, r_out_arcsec, alpha, beta_par)

    mag = -2.5 * np.log10(lum)

    return mag


def ferrer_total_luminosity(Sigma0, r_out, alpha, beta_par):
    """
    Total luminosity integrated from 0 to r_out.
    """
    gamma = 2.0 - beta_par

    if gamma <= 0:
        raise ValueError("The parameter 2 - beta must be positive.")

    a = 2.0 / gamma
    b = alpha + 1.0

    Ltot = (2.0 * np.pi * Sigma0 * r_out**2 / gamma) * beta(a, b)

    return Ltot


def ferrers_total_magnitude(
    mu0_mag_arcsec2: float,
    r_out: float,
    alpha: float,
    beta_par: float,
    pixscale: float,
) -> float:
    """
    Compute the total magnitude of the GALFIT ferrers model.

    Parameters
    ----------
    mu0_mag_arcsec2 : float
        Central surface brightness in mag/arcsec^2.
    r_out: float
        outer truncation radius. in pixels
    alpha: float
         outer truncation shapness
    beta_par: float
        central slope.
    pixscale: float
        Pixel scale in arcsec/pixel.


    Returns
    -------
    float
        Total magnitude of the ferrers model.

    """

    r_out_arcsec = r_out * pixscale

    Sigma0 = 10 ** ((-mu0_mag_arcsec2) / 2.5)

    lumtot = ferrer_total_luminosity(Sigma0, r_out_arcsec, alpha, beta_par)

    mag_total = -2.5 * np.log10(lumtot)

    return mag_total


def magFerrers(galfile, numferrer=2):
    """Computes the magnitude of the Ferrers function.

    It computes magnitude
    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    numferrer: int
            component number (position in input file) of the ferrer function
            default = 2

    Returns
    -------
    magedge: magnitude of the EdgeDisk function

    See Also
    --------
    magEdge: Computes the Ferrers  function.


    """

    galfit = Galfit(galfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # pixscale=.262
    pixscale = galhead.scale

    num_comp = 1  # it assumes the galaxy is the first component. Not related to numexp
    dis = 3

    comps = SelectGal(galcomps, dis, num_comp)

    maskgal = comps.Active == 1

    numferrer = numferrer - 1  # changing to index value

    mu0_mag_arcsec2 = comps.Mag[maskgal][numferrer]

    r_out = comps.Rad[maskgal][numferrer]
    alpha = comps.Exp[maskgal][numferrer]
    beta_par = comps.Exp2[maskgal][numferrer]

    magferrers = ferrers_total_magnitude(
        mu0_mag_arcsec2, r_out, alpha, beta_par, pixscale
    )

    return magferrers


def mainmagFerrers(argv=None) -> int:
    printWelcome()

    parser = argparse.ArgumentParser(
        description="computes the magnitud of the Ferrers function"
    )
    parser.add_argument("galfile", help="GALFIT input file")

    parser.add_argument(
        "-nf",
        "--numferrer",
        type=int,
        default=2,
        help="component number of the ferrer function in GALFIT file. Default=2",
    )

    args = parser.parse_args(argv)

    mag = magFerrers(args.galfile, args.numferrer)

    print("magnitude of Ferrers function: ")
    print(f"  {mag:.2f}")

    return 0


if __name__ == "__main__":
    mainmagFerrers()
