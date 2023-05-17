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

# computes sky 

def main():


    if len(sys.argv[1:]) != 2 and len(sys.argv[1:]) != 6 and len(sys.argv[1:]) != 3:
        print ('Missing arguments')
        print ("Usage:\n %s [ImageFile] [MaskImage] [optional Ds9RegFile] [optional: Xmin Xmax Ymin Ymax]" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits ds9.reg" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits 330 450 200 700" % sys.argv[0])

        sys.exit()


    flagseg=False
    flagreg=False
    xmin=xmax=ymin=ymax=0
    if len(sys.argv[1:]) == 6:
        flagseg=True
        xmin=int(sys.argv[3])
        xmax=int(sys.argv[4])
        ymin=int(sys.argv[5])
        ymax=int(sys.argv[6])

    if len(sys.argv[1:]) == 3:
        flagreg=True
        filereg=sys.argv[3]
        xlo,xhi,ylo,yhi=GetRegionDs9(filereg)

   
    imgname= sys.argv[1]
    maskimage= sys.argv[2]

    ################################################
    ################################################


    hdu = fits.open(imgname)
    imgdat = hdu[0].data.astype(float)
    hdu.close()

    hdu = fits.open(maskimage)
    maskdat = hdu[0].data
    hdu.close()

    if flagseg or flagreg:
        imgdat  = imgdat[ylo - 1:yhi, xlo - 1:xhi]
        maskdat = maskdat[ylo - 1:yhi, xlo - 1:xhi]


    maskdat=np.array(maskdat,dtype=bool)



    mask = maskdat == False

    img=imgdat[mask]  



    print("mean sky: {:.3f} ".format(img.mean()))
    print("std sky: {:.3f} ".format(img.std()))
    #print("rms sky: {:.3f} ".format(rms(img)))


    print("Excluding the top and bottom 20%:") 
 
    flatimg=img.flatten()  
    flatimg.sort()

    tot=len(flatimg)

    top=round(.8*tot)
    bot=round(.2*tot)

    img2=flatimg[bot:top]

    mean=img2.mean()
    sig=img2.std()

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))

    print("Sky within 5 sigma:") 

    img3=img2.copy()

    mask2  = np.abs(img3 - mean) <= 5* sig 

    print("mean sky: {:.3f} ".format(img3[mask2].mean()))
    print("std sky: {:.3f} ".format(img3[mask2].std()))
    



####################################
#####################################
######################################

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
