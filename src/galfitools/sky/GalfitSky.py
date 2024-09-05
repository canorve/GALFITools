#!/usr/bin/env python3

import os
import os.path
import subprocess as sp

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.MakeMask import ds9satbox
from galfitools.mge.mge2galfit import PrintHeader

# change it and use skyRem instead


def galfitSky(imgname, maskfile, mgzpt, scale, X, Y, sky) -> None:
    """Computes the sky using GALFIT

    It uses the FITS image and a mask to compute the
    to make an initial parameter file to compute sky background
    by GALFIT

    Parameters
    ----------
    imgname : str
            name of the FITS image

    maskfile : str
            name of the mask image
    mgzpt: float
            magnitud zeropoint
    scale : float
            plate scale
    X, Y : float, float x,y position of center
    sky : float
        initial parameter of sky

    Notes
    -----
    The function returns a GALFIT initial parameter file
    named sky.txt. The user must run GALFIT with this file
    to obtain the sky background


    """
    flagpos = False

    (NCol, NRow) = GetAxis(imgname)

    convbox = 100

    #    sky=0.55
    fit = 1
    Z = 0

    #    scale = 0.0455  # arcsec/pixel

    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    #    imgname = "ngc4342.fits"
    root_ext = os.path.splitext(imgname)
    TNAM = root_ext[0]

    ##################
    #################

    # Here we use an accurate four gaussians MGE PSF for
    # the HST/WFPC2/F814W filter, taken from Table 3 of
    # Cappellari et al. (2002, ApJ, 578, 787)

    #    sigmapsf = [0.494, 1.44, 4.71, 13.4]      # In PC1 pixels
    #    normpsf = [0.294, 0.559, 0.0813, 0.0657]  # total(normpsf)=1

    parfile = "sky.txt"
    outname = TNAM

    rmsname = "none"
    psfname = "none"

    consfile = "none"

    T1 = "{}".format(imgname)
    T2 = outname + "-sky.fits"
    T3 = "{}".format(rmsname)

    if flagpos is False:
        xlo = 1
        ylo = 1

        (xhi, yhi) = GetAxis(imgname)

    plate = 1

    fout1 = open(parfile, "w")

    PrintHeader(
        fout1,
        T1,
        T2,
        T3,
        psfname,
        1,
        maskfile,
        consfile,
        xlo,
        xhi,
        ylo,
        yhi,
        convbox,
        convbox,
        mgzpt,
        plate,
        plate,
        "regular",
        0,
        0,
    )

    index = 0

    index += 1

    PrintSky(fout1, index, sky, Z, fit)

    fout1.close()

    print("To compute sky run: galfit sky.txt")


####################################################

# GALFIT functions


def PrintSky(hdl, ncomp, sky, Z, fit):
    """Print GALFIT sky function to filehandle

    # repeated
    """

    # k Check

    line00 = "# Object number: {}   \n".format(ncomp)
    line01 = " 0)      sky       #    Object type  \n"
    line02 = " 1) {}     {}  # sky background  \n".format(sky, fit)
    line03 = " 2) 0.000      0    # dsky/dx (sky gradient in x) \n"
    line04 = " 3) 0.000      0     # dsky/dy (sky gradient in y)\n"
    line05 = " Z) {}        # Skip this model in output image?  \n".format(Z)
    line06 = "\n"
    line07a = "================================="
    line07b = "===============================================\n"
    line07 = line07a + line07b

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)

    return True


#############################################################################
#    End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
