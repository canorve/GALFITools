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

# computes sky using GALFIT
#change it and use skyRem instead

def galfitSky(imgname, maskfile, mgzpt, scale, X, Y, sky)-> None:
    

    flagpos=False




    (NCol, NRow) = GetAxis(imgname)



    convbox=100

    #    sky=0.55
    fit=1
    Z=0

    xpos = X
    ypos = Y
    anglegass=35.8


    #    scale = 0.0455  # arcsec/pixel

    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    #    imgname = "ngc4342.fits"
    root_ext = os.path.splitext(imgname)
    TNAM = root_ext[0]
 

    exptime = GetExpTime(imgname)

    ##################
    #################

    # Here we use an accurate four gaussians MGE PSF for
    # the HST/WFPC2/F814W filter, taken from Table 3 of
    # Cappellari et al. (2002, ApJ, 578, 787)

    #    sigmapsf = [0.494, 1.44, 4.71, 13.4]      # In PC1 pixels
    #    normpsf = [0.294, 0.559, 0.0813, 0.0657]  # total(normpsf)=1

    parfile="sky.txt"
    outname=TNAM

    rmsname="none"
    psfname="none"

    consfile="none"

    T1 = "{}".format(imgname)
    T2 = outname + "-sky.fits"
    T3 = "{}".format(rmsname)

    if flagpos == False:
        xlo=1
        ylo=1

        (xhi,yhi)=GetAxis(imgname)

    plate=1

    fout1 = open(parfile, "w")

    PrintHeader(fout1, T1, T2, T3, psfname, 1, maskfile, consfile, xlo, xhi, ylo,
                yhi, convbox, convbox, mgzpt, plate, plate, "regular", 0, 0)

    index = 0

    index+=1



    PrintSky(fout1, index, sky, Z, fit)

    fout1.close()


    print("To compute sky run: galfit sky.txt")





####################################
#####################################
######################################

def ds9satbox (satfileout,output,satscale,satoffset):
    "Creates a file for ds9 which selects bad saturated regions"

    scaleflag=1
    offsetflag=1
    regfileflag=1
    magflag=1
    clasflag=1

    flagsat=4      ## flag value when object is saturated (or close to)
    maxflag=128    ## max value for flag
    check=0
    regflag = 0    ## flag for saturaded regions




    f_out = open(satfileout, "w")

    N,Alpha,Delta,X,Y,Mg,Kr,Fluxr,Isoa,Ai,E,Theta,Bkgd,Idx,Flg=np.genfromtxt(output,delimiter="",unpack=True)


    line="image \n"
    f_out.write(line)



    for idx, item in enumerate(N):


        bi=Ai[idx]*(1-E[idx])

        Theta[idx] = Theta[idx] * np.pi /180  #rads!!!


        Rkronx = satscale * 2 * Ai[idx] * Kr[idx]  + satoffset
        Rkrony = satscale * 2 * bi * Kr[idx]  + satoffset



        if Rkronx == 0:
            Rkronx = 1

        if Rkrony == 0:
            Rkrony = 1

        check=CheckFlag(Flg[idx],flagsat)  ## check if object has saturated regions

        if (check):

            line="box({0},{1},{2},{3},0) # color=red move=0 \n".format(X[idx],Y[idx],Rkronx,Rkrony)
            f_out.write(line)

            line2="point({0},{1}) # point=boxcircle font=\"times 10 bold\" text={{ {2} }} \n".format(X[idx],Y[idx],N[idx])
            f_out.write(line2)


    f_out.close()





def MakeMask(maskimage, catfile, scale, offset, regfile):
    "Create a mask image using ellipses for every Object of catfile. Now includes offset"
# k Check

    checkflag = 0
    flagsat = 4  # flag value when object is saturated (or close to)
    maxflag = 128  # max value for flag

    regflag = 0  # flag for saturaded regions

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg, sxmin, sxmax, symin, symax, sxsmin, sxsmax, sysmin, sysmax = np.genfromtxt(
        catfile, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

    print("Creating Masks for sky \n")

    Rkron = scale * ai * kr + offset

    mask = Rkron < 1
    if mask.any():
        Rkron[mask] = 1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    for idx, val in enumerate(n):

        # check if object doesn't has saturaded regions
        checkflag = CheckFlag(flg[idx], flagsat)
        # check if object is inside of a saturaded box region indicated by
        # user in ds9
        regflag = CheckSatReg(xx[idx], yy[idx], regfile, Rkron[idx], theta[idx], e[idx])

        if (checkflag == False) and (regflag == False):

            print ("Creating ellipse mask for object {}  \n".format(n[idx]))
            img = MakeKron(img, n[idx], xx[idx], yy[idx], Rkron[idx], theta[idx], e[
                idx], sxsmin[idx], sxsmax[idx], sysmin[idx], sysmax[idx])

        elif(checkflag == True or regflag == True):

            print ("Skipping object {}, one or more pixels are saturated \n".format(n[idx]))

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True


def MakeKron(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine create a Kron ellipse within a box defined by: xmin, xmax, ymin, ymax"

# Check

    xmin = np.int(xmin)
    xmax = np.int(xmax)
    ymin = np.int(ymin)
    ymax = np.int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist < dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat


def MakeSatBox(maskimage, region, val, ncol, nrow):
    "Create a mask for saturated regions"
    "Regions must be in DS9 box regions format"

# k Check

#	fileflag=1

    hdu = fits.open(maskimage)
    img = hdu[0].data

    with open(region) as f_in:

        next(f_in)
#        next(f_in)
#        next(f_in)

        # All lines including the blank ones
        lines = (line.rstrip() for line in f_in)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        for line in lines:

            (box, info) = line.split('(')

            if(box == "box"):

                (xpos, ypos, xlong, ylong, trash) = info.split(',')

                xpos = float(xpos)
                ypos = float(ypos)
                xlong = float(xlong)
                ylong = float(ylong)

                xlo = (xpos - xlong / 2)
                xhi = (xpos + xlong / 2)

                ylo = (ypos - ylong / 2)
                yhi = (ypos + ylong / 2)

                xlo = int(xlo)
                xhi = int(xhi)

                ylo = int(ylo)
                yhi = int(yhi)

                if (xlo < 1):

                    xlo = 1

                if (xhi > ncol):

                    xhi = ncol

                if (ylo < 1):

                    ylo = 1

                if (yhi > nrow):

                    yhi = nrow

                img[ylo - 1:yhi, xlo - 1:xhi] = val

    hdu[0].data = img
    hdu.writeto(maskimage, overwrite=True)
    hdu.close()

    return True




def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"
# k Check
    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True


def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow


def CatArSort(SexCat,scale,SexArSort,NCol,NRow):
    # k Check

    # sort the sextractor
    # catalog by magnitude,
    # get sizes for objects
    # and write it in a new file

    # The sextractor catalog must contain the following parameters:
    #   1 NUMBER                 Running object number
    #   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
    #   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
    #   4 X_IMAGE                Object position along x                                    [pixel]
    #   5 Y_IMAGE                Object position along y                                    [pixel]
    #   6 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   7 KRON_RADIUS            Kron apertures in units of A or B
    #   8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
    #   9 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
    #  10 A_IMAGE                Profile RMS along major axis                               [pixel]
    #  11 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
    #  12 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
    #  13 BACKGROUND             Background at centroid position                            [count]
    #  14 CLASS_STAR             S/G classifier output
    #  15 FLAGS                  Extraction flags


    print("Sorting and getting sizes for objects \n")

    n, alpha, delta, xx, yy, mg, kr, fluxrad, ia, ai, e, theta, bkgd, idx, flg = np.genfromtxt(
        SexCat, delimiter="", unpack=True)

    n = n.astype(int)
    flg = flg.astype(int)

#    ai = ai.astype(float)
#    kr = kr.astype(float)

#    scale = scale.astype(float)


    Rkron = scale * ai * kr

    Rwsky = scale * ai * kr + 10  + 20


#   considering to use only  KronScale instead of SkyScale
#    Rwsky = parvar.KronScale * ai * kr + parvar.Offset + parvar.SkyWidth

    Bim = (1 - e) * Rkron

    Area = np.pi * Rkron * Bim *(-1)

    (sxmin, sxmax, symin, symax) = GetSize(xx, yy, Rkron, theta, e, NCol, NRow)

    (sxsmin, sxsmax, sysmin, sysmax) = GetSize(xx, yy, Rwsky, theta, e, NCol,NRow)

#    print(sxmin.size, sxmax.size, symin, symax)

    f_out = open(SexArSort, "w")

    index = Area.argsort()
    for i in index:

        line = "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(n[i], alpha[i], delta[i], xx[i], yy[i], mg[i], kr[i], fluxrad[i], ia[i], ai[i], e[i], theta[i], bkgd[i], idx[i], flg[i], np.int(
            np.round(sxmin[i])), np.int(np.round(sxmax[i])), np.int(np.round(symin[i])), np.int(np.round(symax[i])), np.int(np.round(sxsmin[i])), np.int(np.round(sxsmax[i])), np.int(np.round(sysmin[i])), np.int(np.round(sysmax[i])))

        f_out.write(line)

    f_out.close()

    return len(n)


def GetSize(x, y, R, theta, ell, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    # k Check
    q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

# getting size

    xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        xmin[mask] = 1

    mask = xmax > ncol

    if mask.any():
        xmax[mask] = ncol

    mask = ymin < 1
    if mask.any():
        ymin[mask] = 1

    mask = ymax > nrow
    if mask.any():
        ymax[mask] = nrow

    return (xmin, xmax, ymin, ymax)

#def CheckFlag(val,check,max):
 #  "Check for flag contained in $val, returns 1 if found "

  # flag = False
  # mod = 1

   #while (mod != 0):


    #   res = int(val/max)

    #3  if (max == check and res == 1 ):

    #       flag=True


     #  mod = val % max

      # val = mod
       #max = max/2


#   return flag


def CheckSatReg(x,y,filein,R,theta,ell):
   "Check if object is inside of saturated region. returns True if at least one pixel is inside"
## check if object is inside of
## saturaded region as indicated by ds9 box region
## returns 1 if object center is in saturaded region


   q = (1 - ell)

   bim = q * R

   theta = theta * np.pi /180  ## Rads!!!

   flag = False
   fileflag =1



   with open(filein) as f_in:

       lines = (line.rstrip() for line in f_in) # All lines including the blank ones
       lines = (line.split('#', 1)[0] for line in lines) # remove comments
       lines = (line.rstrip() for line in lines)   # remove lines containing only comments
       lines = (line for line in lines if line) # Non-blank lines

       for line in lines:


           if (line != "image"):

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

                   dx = xpos - x
                   dy = ypos - y

                   landa=np.arctan2( dy,dx )

                   if landa < 0:
                       landa=landa + 2 * np.pi


                   landa = landa - theta

                   angle = np.arctan2(np.sin(landa)/bim, np.cos(landa)/R)

                   xell =  x + R * np.cos(angle)* np.cos(theta)  - bim * np.sin(angle) * np.sin(theta)
                   yell =  y + R * np.cos(angle)* np.sin(theta)  + bim * np.sin(angle) * np.cos(theta)


                   if ( (xell > xlo and xell < xhi) and (yell > ylo and yell < yhi)  ):

                       flag=True
                       break


   return flag













#################333333333333333333333333333333333
####################################################
####################################################
#####################################################





### GALFIT functions

def PrintHeader(hdl, A, B, C, D, E, F, G, xlo, xhi, ylo, yhi, convx, convy, J, platedx, platedy, O, P, S):
    "print GALFIT header in a file"

    # k Check
    # print to filehandle
    # the header for GALFIT

    lineZ = "==================================================================================================\n"
    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}                                   # Input Data image (FITS file)                            \n".format(A)
    lineB = "B) {}                                   # Output data image block                                 \n".format(B)
    lineC = "C) {}                                   # Sigma image name (made from data if blank or \"none\")  \n".format(C)
    lineD = "D) {}                                   # Input PSF image and (optional) diffusion kernel         \n".format(D)
    lineE = "E) {}                                   # PSF fine sampling factor relative to data               \n".format(E)
    lineF = "F) {}                                   # Bad pixel mask (FITS image or ASCII coord list)         \n".format(F)
    lineG = "G) {}                                   # File with parameter constraints (ASCII file)            \n".format(G)
    lineH = "H) {} {} {} {}                          # Image region to fit (xmin xmax ymin ymax)               \n".format(xlo, xhi, ylo, yhi)
    lineI = "I) {} {}                                # Size of the convolution box (x y)                       \n".format(convx, convy)
    lineJ = "J) {}                                   # Magnitude photometric zeropoint                         \n".format(J)
    lineK = "K) {} {}                                # Plate scale (dx dy). \[arcsec per pixel\]               \n".format(platedx, platedy)
    lineO = "O) {}                                   # Display type (regular, curses, both)                    \n".format(O)
    lineP = "P) {}                                   # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps      \n".format(P)
    lineS = "S) {}                                   # Modify/create objects interactively?                    \n".format(S)
    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS                                                     \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:                                      \n"
    line3 = "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,             \n"
    line4 = "#       ferrer, powsersic, sky, and isophote.                                    \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified:                    \n"
    line7 = "#      C0 (diskyness/boxyness),                                                  \n"
    line8 = "#      Fn (n=integer, Azimuthal Fourier Modes),                                  \n"
    line9 = "#      R0-R10 (PA rotation, for creating spiral structures).                     \n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number                                                               \n"
    line12 = "# column 2:                                                                                 \n"
    line13 = "#          -- Parameter 0:    the allowed functions are: sersic, nuker, expdisk             \n"
    line14 = "#                             edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky    \n"
    line15 = "#          -- Parameter 1-10: value of the initial parameters                               \n"
    line16 = "#          -- Parameter C0:   For diskiness/boxiness                                        \n"
    line17 = "#                             <0 = disky                                                    \n"
    line18 = "#                             >0 = boxy                                                     \n"
    line19 = "#          -- Parameter Z:    Outputting image options, the options are:                    \n"
    line20 = "#                             0 = normal, i.e. subtract final model from the data to create \n"
    line21 = "#                             the residual image                                            \n"
    line22 = "#                             1 = Leave in the model -- do not subtract from the data       \n"
    line23 = "#                                                                                           \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0)                                       \n"
    line25 = "# column 4: comment                                                                         \n"
    line26 = " \n"

    line27 = "==================================================================================================\n"

    hdl.write(lineZ)
    hdl.write(lineX)
    hdl.write(lineA)
    hdl.write(lineB)
    hdl.write(lineC)
    hdl.write(lineD)
    hdl.write(lineE)
    hdl.write(lineF)
    hdl.write(lineG)
    hdl.write(lineH)
    hdl.write(lineI)
    hdl.write(lineJ)
    hdl.write(lineK)
    hdl.write(lineO)
    hdl.write(lineP)
    hdl.write(lineS)
    hdl.write(lineY)

    hdl.write(line0)
    hdl.write(line1)
    hdl.write(line2)
    hdl.write(line3)
    hdl.write(line4)
    hdl.write(line5)
    hdl.write(line6)
    hdl.write(line7)
    hdl.write(line8)
    hdl.write(line9)
    hdl.write(line10)

    hdl.write(line11)
    hdl.write(line12)
    hdl.write(line13)
    hdl.write(line14)
    hdl.write(line15)
    hdl.write(line16)
    hdl.write(line17)
    hdl.write(line18)
    hdl.write(line19)
    hdl.write(line20)
    hdl.write(line21)
    hdl.write(line22)
    hdl.write(line23)
    hdl.write(line24)
    hdl.write(line25)
    hdl.write(line26)
    hdl.write(line27)

    return True


def PrintSky(hdl, ncomp, sky, Z, fit):
    "Print GALFIT sky function to filehandle"

    # k Check

    line00 = "# Object number: {}                                                             \n".format(ncomp)
    line01 = " 0)      sky            #    Object type                                        \n"
    line02 = " 1) {}         {}       # sky background        [ADU counts]                    \n".format(sky, fit)
    line03 = " 2) 0.000      0        # dsky/dx (sky gradient in x)                           \n"
    line04 = " 3) 0.000      0        # dsky/dy (sky gradient in y)                           \n"
    line05 = " Z) {}                  # Skip this model in output image?  (yes=1, no=0)       \n".format(Z)
    line06 = "\n"
    line07 = "================================================================================\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)

    return True


def PrintSersic(hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit):
    "print GALFIT Sersic function to filehandle"
    # k Check

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}                                                             \n".format(
            ncomp)
    line01 = " 0)     sersic               #  Object type                                     \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}            #  position x, y     [pixel]                       \n".format(
        xpos, ypos, fit, fit)
    line03 = " 3) {:.2f}       {}              #  total magnitude                                 \n".format(
        magser, fit)
    line04 = " 4) {:.2f}       {}              #  R_e         [Pixels]                            \n".format(
            reser, fit)
    line05 = " 5) {}       {}              #  Sersic exponent (deVauc=4, expdisk=1)           \n".format(
        nser, fit)
    line06 = " 9) {:.2f}       {}              #  axis ratio (b/a)                                \n".format(
        axratser, fit)
    line07 = "10) {:.2f}       {}              #  position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        angleser, fit)
    line08 = " Z) {}                       #  Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line09 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)
    hdl.write(line09)

    return True

def PrintGauss(hdl, ncomp, xpos, ypos, magass, fwhm, axratgass, anglegass, Z, fit):
    "print GALFIT GAUSS function to filehandle"

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}                                                             \n".format(
            ncomp)
    line01 = " 0)     gaussian               #  Object type                                     \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}            #  position x, y     [pixel]                       \n".format(
        xpos, ypos, fit, fit)
    line03 = " 3) {:.2f}       {}              #  total magnitude                                 \n".format(
        magass, fit)
    line04 = " 4) {:.2f}       {}              #  FWHM         [Pixels]                            \n".format(
            fwhm, fit)
    line05 = " 9) {:.2f}       {}              #  axis ratio (b/a)                                \n".format(
        axratgass, fit)
    line06 = "10) {:.2f}       {}              #  position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        anglegass, fit)
    line07 = " Z) {}                       #  Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line08 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)

    return True

def PrintExp(hdl, ncomp, xpos, ypos, magexp, rsexp, axratexp, angleexp, Z, fit):
    "print GALFIT exponential function to filehandle"

    # k Check

    # print to filehandle
    # a exponential function given
    # by the parameters

    line00 = "# Object number: $ncomp                                                        \n"
    line01 = " 0)     expdisk              # Object type                                     \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}           # position x, y     [pixel]                       \n".format(
        xpos, ypos, fit, fit)
    line03 = " 3) {:.2f}        {}             # total magnitude                                 \n".format(
        magexp, fit)
    line04 = " 4) {:.2f}        {}             #      Rs  [Pixels]                               \n".format(
        rsexp, fit)
    line05 = " 9) {:.2f}        {}             # axis ratio (b/a)                                \n".format(
        axratexp, fit)
    line06 = "10) {:.2f}        {}             # position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        angleexp, fit)
    line07 = " Z) {}                       # Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line08 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)

    return True


def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"
# k Check


    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.open(Image)
    dat = hdu[0].data[ylo - 1:yhi, xlo - 1:xhi]
    hdu[0].data = dat
    hdu.writeto(Imageout, overwrite=True)
    hdu.close()

def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow

def GetExpTime(Image):
    # k Check
    "Get exposition time from the image"

    hdu = fits.open(Image)
    exptime = hdu[0].header["EXPTIME"]
    hdu.close()
    return exptime


def CheckFlag(val, check):
    "Check for flag contained in val, returns True if found "

    flag = False
    mod = 1
    maxx=128


    while (mod != 0):

        res = int(val / maxx)

        if (maxx == check and res == 1):

            flag = True

        mod = val % maxx

        val = mod
        maxx = maxx / 2

    return flag


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
    mainGalfitSky()
