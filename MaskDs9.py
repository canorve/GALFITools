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


    parser = argparse.ArgumentParser(description="creates (or modify) a mask image for GALFIT from a Ds9 region file ")

    parser.add_argument("MaskFile", help="the Mask image file to modify or create")
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument("-f","--fill", type=int, help="the value in counts to fill into the Ds9 regions. Default = 0 (remove)",default=0)

    parser.add_argument("-i","--image", type=str, help="image to obtain the size  ")





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


    (ncol, nrow) = GetAxis(MaskFile)


    hdu=fits.open(MaskFile)

    Image = hdu[0].data



    if not os.path.exists(RegFile):
        print ('%s: reg filename does not exist!' %(sys.argv[2]))
        sys.exit()

    v0=[]
    v1=[]
    v2=[]
    v3=[]
    v4=[]
    v5=[]

    f1 = open(RegFile,'r')

    lines = f1.readlines()

    f1.close()

    flag=False
    #reading reg file
    for line in lines:


        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True

        if b1[0] == "box":

            x0 = "box"
            x2 = x1[4:]
            flag = True

        if (flag == True):
            x3 = p[4]
            x4 = x3[:-2]


            v0.append(x0)
            v1.append(float(x2) - 1)
            v2.append(float(p[1]) - 1)
            v3.append(float(p[2]))
            v4.append(float(p[3]))
            v5.append(float(x4))

            flag=False


    obj  = np.array(v0)
    xpos = np.array(v1) 
    ypos = np.array(v2)
    rx = np.array(v3)
    ry = np.array(v4)
    angle = np.array(v5)



    for idx, item in enumerate(obj):

        #converting ellipse from DS9 ell to kron ellipse

        if obj[idx] == "ellipse":

            
            Image = MakeEllip(Image,fill,xpos[idx],ypos[idx],rx[idx],ry[idx],angle[idx],ncol,nrow)


        if obj[idx] == "box":


            Image = MakeBox(Image,fill,xpos[idx],ypos[idx],rx[idx],ry[idx],angle[idx],ncol,nrow)



    #writing mask file

    hdu.data=Image
    hdu.writeto(MaskFile,overwrite=True) 
    hdu.close()


def MakeEllip(Image,fill,xpos,ypos,rx,ry,angle,ncol,nrow):
    "Make an ellipse in an image"

    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos,ypos,rx,ry,angle)
    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)
    Image = MakeKron(Image, fill, xx, yy, Rkron, theta, e, xmin, xmax, ymin, ymax)

    return Image


def MakeBox(Image,fill,xpos,ypos,rx,ry,angle,ncol,nrow):
    "Make a box in an image"

    xlo = xpos - rx/2.
    xhi = xpos + rx/2.

    ylo = ypos - ry/2.
    yhi = ypos + ry/2.

    if xlo  < 1:
        xlo = 1
    if ylo  < 1:
        ylo = 1

    if xhi  > ncol:
        xhi = ncol 
 
    if yhi  > nrow:
        yhi = nrow
 
    xlo=int(np.round(xlo))
    xhi=int(np.round(xhi))

    ylo=int(np.round(ylo))
    yhi=int(np.round(yhi))

    Image[ylo-1:yhi,xlo-1:xhi]=fill

    return Image



def GetAxis(Image):
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]  # for hubble images
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow



def Ds9ell2Kronell(xpos,ypos,rx,ry,angle):


    if rx >= ry:

        q = ry/rx
        e = 1 - q
        Rkron = rx
        theta = angle
        xx = xpos
        yy = ypos
    else:
        q = rx/ry
        e = 1 - q
        Rkron = ry
        theta = angle + 90
        xx = xpos
        yy = ypos
 
     
    return xx, yy, Rkron, theta, e 



def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

    # Check

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1 : ymax + 1, xmin - 1: xmax +1]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2*np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist <= dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


def GetSize(x, y, R, theta, ell, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    # k Check
    q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    constx =  np.sqrt((R**2)*(np.cos(theta))**2 + (bim**2)*(np.sin(theta))**2)
    consty =  np.sqrt((R**2)*(np.sin(theta))**2 + (bim**2)*(np.cos(theta))**2)



    xmin = x - constx                       
    xmax = x + constx
    ymin = y - consty
    ymax = y + consty
                    

    mask = xmin < 1
    if mask.any():
        xmin = 1

    mask = xmax > ncol
    if mask.any():
        xmax = ncol

    mask = ymin < 1
    if mask.any():
        ymin = 1

    mask = ymax > nrow
    if mask.any():
        ymax = nrow

    return (round(xmin), round(xmax), round(ymin), round(ymax))



  
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


