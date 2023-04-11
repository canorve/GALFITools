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

from matplotlib.path import Path

from shapely.geometry import Point, Polygon

### This program is the python version of polyfit.f from
## GALFIT's webpage


def main():


    parser = argparse.ArgumentParser(description="creates (or modify) a mask image for GALFIT using a polygon shape from a Ds9 region file ")

    parser.add_argument("MaskFile", help="the Mask image file to modify or create")
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument("-f","--fill", type=int, help="the value in counts to fill into the Ds9 regions. Default = 0 (remove)",default=0)

    parser.add_argument("-i","--image", type=str, help="image to obtain the size (this img will not be modified")


    args = parser.parse_args()

    MaskFile = args.MaskFile 
    RegFile = args.RegFile 
    fill = args.fill
    image = args.image



    if not os.path.exists(MaskFile):

        print ('%s: image filename does not exist!' %(sys.argv[1]))
        print ('Creating a new image file ')

        hdu=fits.PrimaryHDU()

        if image:
            (ncol, nrow) = GetAxis(image)
        else:
            nx = input("enter numbers of pixels in X ") 
            ny = input("enter numbers of pixels in Y ") 

            nx = np.int64(nx)
            ny = np.int64(ny)

        Image = np.zeros([ny,nx])
        hdu.data=Image
        hdu.writeto(MaskFile,overwrite=True) 





    f1 = open(RegFile,'r')
    lines = f1.readlines()
    f1.close()


    polv = lines[3].split(")")
    pol,ver = polv[0].split("(")

    if pol == "polygon":

        points = ver.split(",")

        N=len(points)

        i=0

        x=np.array([])
        y=np.array([])

        x=x.astype(float)
        y=y.astype(float)

        tupVerts=[]

        #make tuple for vertices
        while i < N:

            tupVerts = tupVerts + [(round(float(points[i])-1),round(float(points[i+1])-1)) ]
            i=i+2

        print("vertices: ", tupVerts)

    else:

        print("can't find a polygon in ds9 region file")
        sys.exit()


    poly = Polygon(tupVerts)

    (ncol, nrow) = GetAxis(MaskFile)

    hdu = fits.open(MaskFile)

    image = hdu[0].data

    x, y = np.meshgrid(np.arange(ncol), np.arange(nrow)) # make a canvas with coordinates
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T 
    p = Path(tupVerts) # make a polygon

    grid = p.contains_points(points)
    mask = grid.reshape(nrow, ncol) # now you have a mask with points inside a polygon


    image[mask] = fill

    #imgdata = ~mask*image
    imgdata = image

    #ypos, xpos = np.mgrid[0: nrow, 0: ncol]

    #writing mask file

    #hdu[0].data = imgdata.astype(int)
    hdu[0].data = imgdata
    hdu.writeto(MaskFile,overwrite=True) 
    hdu.close()




def GetAxis(Image):
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]  # for hubble images
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow

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
