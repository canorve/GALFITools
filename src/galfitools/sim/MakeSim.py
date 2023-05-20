#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits


#check modidy argparse 
def main():

    if len(sys.argv[1:]) != 5:
        print ('Missing arguments')
        print ("Usage:\n %s [ImageFile] [GAIN] [skymean] [skystd] [newimage]" % sys.argv[0])
        print ("Example:\n %s mod.fits 10 0.11 3.2 sim.fits " % sys.argv[0])
        sys.exit()

    image = sys.argv[1]
    GAIN = float(sys.argv[2])

    skymean= float(sys.argv[3])
    skystd =  float(sys.argv[4])

    newimage = sys.argv[5]

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




