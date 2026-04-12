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


def magCorr(
    galfitFile: str,
    fileout: str,
    A=0.0,
    K=0.0,
) -> GalComps:
    """Corrects magnitudes for Extinction and K correction.

    It reads a galfit file and corrects all the magnitudes
    for Extinction and K correction

    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    fileout: str
            name of the GALFIT output file
    A: float
       magnitud correction for Extinction
    K: float
       K Correction

    Returns
    -------
    galcomps: GalComps
              data class defined in galfitools.galin.galfit

    Note
    ----
    user must provide the corrections and applies to all components
    in the GALFIT file.

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

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    # comps = conver2Sersic(galcomps)

    # makes corrections:
    galcomps.Mag = galcomps.Mag - A - K

    # printing output file
    fout = open(fileout, "w")

    galPrintHeader(fout, galhead)
    index = 0
    for index, item in enumerate(galcomps.N):

        galPrintComp(fout, index + 1, index, galcomps)

    galPrintSky(fout, index + 1, galsky)

    fout.close()

    return galcomps
