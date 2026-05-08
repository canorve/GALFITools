#! /usr/bin/env python3

import numpy as np
import argparse

from galfitools.shell.prt import printWelcome
from galfitools.galin.galfit import (
    GalComps,
    Galfit,
    GalHead,
    GetRadAng,
    SelectGal,
    conver2Sersic,
    conver2Edge,
    numComps,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)


def edgedisk_total_magnitude(
    mu0_mag_arcsec2: float,
    hs_pix: float,
    rs_pix: float,
    pixscale: float,
) -> float:
    """
    Compute the total magnitude of a GALFIT edgedisk model.

    Parameters
    ----------
    mu0_mag_arcsec2 : float
        Central surface brightness in mag/arcsec^2.
    hs_pix : float
        Vertical scale height in pixels.
    rs_pix : float
        Radial scale length in pixels.
    pixscale : float
        Pixel scale in arcsec/pixel.

    Returns
    -------
    float
        Total magnitude of the edgedisk model.

    Raises
    ------
    ValueError
        If any geometric parameter is non-positive.
    """
    if hs_pix <= 0:
        raise ValueError("hs_pix must be > 0.")
    if rs_pix <= 0:
        raise ValueError("rs_pix must be > 0.")
    if pixscale <= 0:
        raise ValueError("pixscale must be > 0.")

    area_factor = 2.0 * np.pi * rs_pix * hs_pix * pixscale**2
    mag_total = mu0_mag_arcsec2 - 2.5 * np.log10(area_factor)

    return mag_total


def magEdge(galfile, numedge=2):
    """Computes the magnitude of the EdgeDisk function.

    It converts a exponential function (or Sersic with n = 1)
    to a edgedisk function

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    numexp: int
            component number (position in input file) of the exponential function
            default = 2

    Returns
    -------
    magedge: magnitude of the EdgeDisk function

    See Also
    --------
    exp2edge: Convert a Exponential function to EdgeDisk function.


    """

    galfit = Galfit(galfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # pixscale=.262
    pixscale = galhead.scale

    num_comp = (
        1  # it assumes the edge-on galaxy is the first component. Not related to numexp
    )
    dis = 3

    comps = SelectGal(galcomps, dis, num_comp)

    # z
    # mu0_mag_arcsec2=  18.4556
    # hs_pix=  3.4134
    # rs_pix= 6.8061

    maskgal = comps.Active == 1

    mu0_mag_arcsec2 = comps.Mag[maskgal][numedge]

    hs_pix = comps.Rad[maskgal][numedge]
    rs_pix = comps.Exp[maskgal][numedge]

    magedge = edgedisk_total_magnitude(mu0_mag_arcsec2, hs_pix, rs_pix, pixscale)

    return magedge


def mainmagEdge(argv=None) -> int:
    printWelcome()

    parser = argparse.ArgumentParser(
        description="computes the magnitud of a edgedisk function"
    )
    parser.add_argument("galfile", help="GALFIT input file")

    parser.add_argument(
        "-ne",
        "--numedge",
        type=int,
        default=2,
        help="component number of the edgedisk function in GALFIT file. Default=2",
    )

    args = parser.parse_args(argv)

    mag = magEdge(args.galfile, args.numedge)

    print("magnitude of the edgedisk: ")
    print(f"  {mag:.2f}")

    return 0


if __name__ == "__main__":
    mainmagEdge()
