#!/usr/bin/env python3


import numpy as np
from astropy.io import fits
import sys
import os

def imarith(ImageFile: str, output: str, image2: str, add: float, mul: float, div: float, sub: float)-> None: 

    if not os.path.exists(ImageFile):
        print ('{} image filename does not exist!'.format(ImageFile))
        sys.exit()

    if image2:
        if not os.path.exists(image2):
            print ('{} image filename does not exist!'.format(image2))
            sys.exit()

        hdu2 = fits.open(image2)
        dataImage2 = hdu2[0].data




    # original file
    hdu = fits.open(ImageFile)
    dataImage = hdu[0].data


    if add:
        if image2:
            newdata = dataImage + dataImage2
        else:
            newdata = dataImage + add

    if mul:
        if image2:
            newdata = dataImage*dataImage2
        else:
            newdata = dataImage*mul

    if div:
        if image2:
            newdata = dataImage/dataImage2
        else:
            newdata = dataImage/div

    if sub:
        if image2:
            newdata = dataImage - dataImage2
        else:
            newdata = dataImage - sub



    hdu[0].data = newdata 
    hdu.writeto(output, overwrite=True)


    hdu.close()

    if image2:
        hdu2.close()


