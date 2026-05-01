#!/usr/bin/env python3

import os
import os.path
import subprocess as sp
import sys

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import MakeImage


def skyRem(imageFile, mean, sig, nsig, output="mask.fits", bor_flag=False, borValue=0):
    """Creates a mask image for GALFIT by subtracting the sky background.

    Parameters
    ----------
    imageFile : str
                FITS image where the data will be taken
    mean : float
            mean of the sky  background
    sig : float
          standard deviation of sky
    nsig : float
          number of sky standard deviations to be removed from image
    output: str
                output name of the new mask image
    bor_flag : False, optional
               if True, it will mask the border of the image. This is
               for those region where the image matrix is larger than the
               data matrix, e.g.  Hubble images
    borValue : float
               value of the border in the image

    Returns
    -------
    None

    """

    bor_val = 100

    if not os.path.exists(imageFile):
        print("{} image filename does not exist!".format(imageFile))
        sys.exit()

    (ncol, nrow) = GetAxis(imageFile)

    if not os.path.exists(output):
        print("{} mask image does not exist, creating one ... ".format(output))
    else:
        print("overwriting {} mask image ".format(output))

    MakeImage(output, ncol, nrow)

    i = 0  # index of data
    # original file
    hdu = fits.open(imageFile)
    dataImage = hdu[i].data

    # output file
    hdu2 = fits.open(output)
    maskImage = hdu2[0].data

    topsky = mean + sig * nsig

    mask = dataImage >= topsky

    if mask.any():
        maskImage[mask] = dataImage[mask] - topsky

    # masking the border in case:

    if bor_flag:
        print("masking the border")
        bor_mask = dataImage == borValue

        if bor_mask.any():
            maskImage[bor_mask] = bor_val

    # writing mask file:

    hdu2[0].data = maskImage

    hdu2.writeto(output, overwrite=True)

    hdu.close()
    hdu2.close()


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
