#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits

import argparse

def main():


    parser = argparse.ArgumentParser(description="simulates a observed galaxy from a GALFIT model")

    parser.add_argument("image", help="the GALFIT galaxy model")
    parser.add_argument("newimage", help="the name of the new galaxy image")


    parser.add_argument("-s","--sky", type=float, help="the sky background value. default = 0", default=0)
    parser.add_argument("-std","--std", type=float, help="the sky standard deviation. default = 1", default=1)

    parser.add_argument("-g","--gain", type=float, help="the gain value of the image. default =  1", default=1)



    args = parser.parse_args()

    image = args.image 
    GAIN = args.gain 

    skymean = args.sky
    skystd = args.std 

    newimage = args.newimage 




    MakeSim(image, GAIN, skymean, skystd, newimage)




def MakeSim(image, GAIN, skymean, skystd, newimage)-> None: 


    sizex,sizey =GetAxis(image) 

    hdu = fits.open(image)

    img=hdu[0].data

    eimg = GAIN*img

    noisyimg= np.random.poisson(eimg)

    pimg = noisyimg / GAIN

    hdu[0].data = pimg  

    hdu.writeto("poissonimg.fits",overwrite=True)

    sky = np.random.normal(skymean,skystd,(sizex,sizey))
   

    hdusky = fits.PrimaryHDU()
    hdusky.data = sky

    hdusky.writeto("skynoise.fits",overwrite=True)



    #sumar las dos partes:

    imgsim = pimg + sky


    hdu[0].data = imgsim

    hdu.writeto(newimage,overwrite=True)


    hdu.close()


def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow


#end of program
if __name__ == '__main__':
    main()




