#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits

import argparse

from galfitools.galin.MaskDs9 import GetAxis 

def makeSim(image, GAIN, skymean, skystd, newimage)-> None: 


    sizex, sizey = GetAxis(image) 

    hdu = fits.open(image)

    img = hdu[0].data

    eimg = GAIN*img

    noisyimg = np.random.poisson(eimg)

    pimg = noisyimg/GAIN

    hdu[0].data = pimg  

    hdu.writeto("poissonimg.fits", overwrite=True)

    sky = np.random.normal(skymean, skystd, (sizex,sizey))
   

    hdusky = fits.PrimaryHDU()
    hdusky.data = sky

    hdusky.writeto("skynoise.fits", overwrite = True)



    #sumar las dos partes:

    imgsim = pimg + sky


    hdu[0].data = imgsim

    hdu.writeto(newimage, overwrite = True)


    hdu.close()



#end of program
if __name__ == '__main__':
    main()




