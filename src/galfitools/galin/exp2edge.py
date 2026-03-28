#! /usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy
import os
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


def Exp2Edge(
    galfitFile: str,
    fileout: str,
    numexp: int,
) -> GalComps:
    """Convert a Exponential function to EdgeDisk function.

    It converts a exponential function (or Sersic with n = 1)
    to a edgedisk function

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    fileout: str
            name of the GALFIT output file
    numexp: int
            component number (position in input file) of the exponential function

    Returns
    -------
    galcomps: GalComps
              data class defined in galfitools.galin.galfit

    See Also
    --------
    conver2Sersic : converts gaussian, exponential de Vaucoulers
                    to Sersic models

    conver2Ferrer: convert gaussian bar to a Ferrer bar


    """

    galfit = Galfit(galfitFile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    num_comp = 1  # it assumes the barred galaxy is the first component
    dis = 3

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    comps = SelectGal(comps, dis, num_comp)

    comps = conver2Edge(comps, galhead.scale, numexp)  # converts only the 2 component

    # printing output file
    fout = open(fileout, "w")

    filename = galhead.outimage
    root, extension = os.path.splitext(filename)
    newname = root + "-edge.fits"
    galhead.outimage = newname

    galPrintHeader(fout, galhead)

    index = 0

    for index, item in enumerate(comps.N):
        galPrintComp(fout, index + 1, index, comps)

    galPrintSky(fout, index + 1, galsky)

    fout.close()

    return comps
