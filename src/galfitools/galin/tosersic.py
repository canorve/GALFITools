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


def toSersic(
    galfitFile: str,
    fileout: str,
    nfree: bool,
) -> GalComps:
    """Convert to Sersic  components.

    It converts a exponential, de Vaucouleurs, Gaussian functions
    to a Sersic function

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    fileout: str
            name of the GALFIT output file
    nfree: Boolean
          leaves sersic index as free to fit


    Returns
    -------
    galcomps: GalComps
              data class defined in galfitools.galin.galfit

    Note
    ----
    the components to convert to Sersic must be
    derived from Sersic: Gaussian, de Vaucoulers, exponential.
    This does not apply to Ferrer, Nuker, King, etc.

    See Also
    --------
    conver2Sersic : converts gaussian, exponential de Vaucoulers
                    to Sersic models

    conver2Ferrer: convert gaussian bar to a Ferrer bar

    conver2Edge: convert to edgedisk

    """

    galfit = Galfit(galfitFile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    # num_comp = 1
    # dis = 3

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    # printing output file
    fout = open(fileout, "w")

    galPrintHeader(fout, galhead)
    index = 0
    for index, item in enumerate(comps.N):

        if nfree:
            comps.ExpFree[index] = 1
        else:
            comps.ExpFree[index] = 0

        galPrintComp(fout, index + 1, index, comps)

    galPrintSky(fout, index + 1, galsky)

    fout.close()

    return comps
