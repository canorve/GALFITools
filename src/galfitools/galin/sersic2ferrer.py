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
    conver2Ferrer,
    numComps,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)


def Sersic2Ferrer(
    galfitFile: str,
    alpha: bool,
    beta: bool,
    fileout: str,
) -> GalComps:
    """Convert a Sersic function to Ferrer funciton.


    It converts a galfit (prefered fitted) file which
    contains of a Bulge (Sersic), Bar (gaussian/sersic n = 0.5)
    and disk (exponential/sersic n =1) surface brightness model
    (in that order) into a file containing  Bulge (Sersic no
    changes), Bar (Ferrer alpha = 1, beta = 0) and disk
    (exponential/sersic n =1 no changes).


    Parameters
    ----------
    galfitFile: str
            name of the GALFIT file
    alpha: bool
           If True, it leaves the Ferrer alpha parameter as free

    beta: bool
           If True, it leaves the Ferrer beta parameter as free

    fileout: str
            name of the GALFIT output file

    Returns
    -------
    galcomps: GalComps
              data class defined in galfitools.galin.galfit

    See Also
    --------
    conver2Sersic : converts gaussian, exponential de Vaucoulers
                    to Sersic models

    Notes
    -----
    To make the convertion from Sersic (n=0.5) to Ferrer
    it takes the first two terms of an Taylor expantion
    of the Sersic to make the approximation to Ferrer.

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

    # taking the last component position angle for the whole galaxy

    maskgal = comps.Active == 1

    if alpha:
        comps.ExpFree = 1
    else:
        comps.ExpFree = 0

    if beta:
        comps.Exp2Free = 1
    else:
        comps.Exp2Free = 0

    comps[1] = conver2Ferrer(comps[1])  # converts only the 2 component

    # printing output file
    fout = open(fileout, "w")

    galPrintHeader(fout, head)

    index = 0

    for index, item in enumerate(comps.N):

        galPrintComp(fout, index + 1, index, comps)

    galPrintSky(fout, index + 1, galsky)

    fout.close()

    return comps
