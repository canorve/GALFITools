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

#introducir sky box y sky ring
#use maskds9 para obtener los pixeles de las regiones

def sky(imgname, maskimage, filereg) -> None:



    xlo,xhi,ylo,yhi = GetRegionDs9(filereg)
    ################################################
    ################################################


    hdu = fits.open(imgname)
    imgdat = hdu[0].data.astype(float)
    hdu.close()

    hdu = fits.open(maskimage)
    maskdat = hdu[0].data
    hdu.close()

    imgdat  = imgdat[ylo - 1:yhi, xlo - 1:xhi]
    maskdat = maskdat[ylo - 1:yhi, xlo - 1:xhi]


    maskdat=np.array(maskdat,dtype=bool)



    mask = maskdat == False

    img=imgdat[mask]  



    #print("mean sky: {:.3f} ".format(img.mean()))
    #print("std sky: {:.3f} ".format(img.std()))
    #print("rms sky: {:.3f} ".format(rms(img)))


    #print("Excluding the top and bottom 20%:") 
 
    flatimg = img.flatten()  
    flatimg.sort()

    tot = len(flatimg)

    top = round(.8*tot)
    bot = round(.2*tot)

    img2 = flatimg[bot:top]

    mean = img2.mean()
    sig = img2.std()


    return mean, sig


    #print("mean sky: {:.3f} ".format(mean))
    #print("std sky: {:.3f} ".format(sig))


    #img3=img2.copy()

    #mask2  = np.abs(img3 - mean) <= 3* sig 

    #msky = img3[mask2].mean()
    #ssky = img3[mask2].std()

    #return msky, ssky 
    



def rms(array):
   return np.sqrt(np.mean(array ** 2))


def GetRegionDs9(filein):
    "Get the size (xmin,xmax,ymin,ymax) from ds9 region file "

    xhi=0
    xlo=0
    ylo=0
    yhi=0

    with open(filein) as f_in:

        lines = (line.rstrip() for line in f_in) # All lines including the blank ones
        lines = (line.split('#', 1)[0] for line in lines) # remove comments
        lines = (line.rstrip() for line in lines)   # remove lines containing only comments
        lines = (line for line in lines if line) # Non-blank lines

        for line in lines:

            (chunks)=line.split(' ')
            if (chunks[0] != "image" and chunks[0] != "physical" and chunks[0] != "global"  ):

                (box,info)=line.split('(')

                if(box == "box"):

                    (xpos,ypos,xlong,ylong,trash)=info.split(',')

                    xpos=float(xpos)
                    ypos=float(ypos)
                    xlong=float(xlong)
                    ylong=float(ylong)


                    xlo = xpos - xlong/2
                    xhi = xpos + xlong/2

                    ylo = ypos - ylong/2
                    yhi = ypos + ylong/2



    xlo=int(round(xlo))
    xhi=int(round(xhi))
    ylo=int(round(ylo))
    yhi=int(round(yhi))


    return xlo,xhi,ylo,yhi




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
