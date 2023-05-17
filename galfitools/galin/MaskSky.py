#!/usr/bin/env python3

import numpy as np
import sys
import os
import stat
import subprocess as sp
import os.path
from astropy.io import fits
import scipy
import scipy.special
import matplotlib.pyplot as plt

import argparse


# programa que hace parches (mascaras) en una seccion de una imagen FITS

def main(): 

    bor_flag = False

    parser = argparse.ArgumentParser(description="creates a mask image for GALFIT using original image and sky mean and sigma")


    parser.add_argument("ImageFile", help="original data image ")
    parser.add_argument("MaskFile", help="Name of the new Mask file")
    parser.add_argument("-sm","--skymean",default=0,type=float, help="mean of the sky background")
    parser.add_argument("-ss","--skysigma",default=0,type=float, help="sigma of the sky background")
    parser.add_argument("-ns","--numbersig",default=1,type=float, help="number of times that the sigma of the sky will be multiplied to remove the sky background")

    parser.add_argument("-b","--border", action="store_true", help="Mask the borders when their value is zero")

    parser.add_argument("-bv","--borValue",default=0,type=float, help="value of the border if it is different from zero")
    args = parser.parse_args()

    image = args.ImageFile 
    mask = args.MaskFile
    sky_mean = args.skymean
    sky_sig = args.skysigma
    nsig = args.numbersig

    bor_flag = args.border
    borValue = args.borValue


    SkyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)

#################################################################
#################################################################

def SkyRem(imageFile,maskFile,mean,sig, nsig,borValue,bor_flag=False):

    bor_val = 100

    if not os.path.exists(imageFile):
        print ('{} image filename does not exist!'.format(imageFile))
        sys.exit()


    (ncol, nrow) = GetAxis(imageFile)



    if not os.path.exists(maskFile):
        print ('{} mask image does not exist, creating one ... '.format(maskFile)) 
    else:
        print ('overwriting {} mask image '.format(maskFile)) 

    MakeImage(maskFile, ncol, nrow)




      # original file
    hdu=fits.open(imageFile)
    dataImage = hdu[0].data


    # output file
    hdu2=fits.open(maskFile)
    maskImage = hdu2[0].data


    topsky = mean + sig*nsig


    mask = dataImage >= topsky 



    if mask.any():
        maskImage[mask] = dataImage[mask] - topsky 


    #masking the border in case:

    if bor_flag:
        print("masking the border")
        bor_mask = dataImage == borValue 

        if bor_mask.any():
            maskImage[bor_mask] = bor_val 
 


    #writing mask file:

    hdu2[0].data=maskImage

    hdu2.writeto(maskFile,overwrite=True)


    hdu.close()
    hdu2.close()




def GetAxis(Image):
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]  # for hubble images
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow




def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex),dtype=np.float64)
    hdu.writeto(newfits, overwrite=True)

    return True




  
#############################################################################
######################### End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
if __name__ == '__main__':
  main()


