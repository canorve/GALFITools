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

import mgefit
from mgefit.find_galaxy import find_galaxy
#from mgefit.mge_fit_1d import mge_fit_1d
from mgefit.sectors_photometry import sectors_photometry
#from mgefit.sectors_photometry_twist import sectors_photometry_twist


from mgefit.mge_fit_sectors import mge_fit_sectors
#from mgefit.mge_print_contours import mge_print_contours
#from mgefit.mge_fit_sectors_twist import mge_fit_sectors_twist


def main():

    if len(sys.argv[1:]) != 5 and len(sys.argv[1:]) != 6:
      print ('Missing arguments')
      print ("Usage: %s [Image] [X] [Y] [MagZpt] [Box Size] [Sky (Optional)] " % (sys.argv[0]))
#      print ("Usage: %s [Image] [MagZpt] [PSF sigma] [twist yes 0=No] [Sky (Optional)] [Mask (optional)] " % (sys.argv[0]))
      print ("Example: %s  image.fits 1000 1050 25 50 0.5" % (sys.argv[0]))
      print ("Example: %s  image.fits 1000 1050 25 100 " % (sys.argv[0]))
      sys.exit()

    flagsky=False
    flageven=False

    imgname= sys.argv[1]

    maskfile ="none"

    X = np.float(sys.argv[2])
    Y = np.float(sys.argv[3])

    mgzpt= sys.argv[4]
    mgzpt=np.float(mgzpt)

    boxsize= np.float(sys.argv[5])
    boxsize=np.int(np.round(boxsize))

    if boxsize % 2 == 0:
        flageven=True

    if len(sys.argv[1:]) == 6:
        flagsky =True
        sky = np.float(sys.argv[6])


#    exptime= sys.argv[3]
#    exptime=np.float(exptime)

    convbox=100

#    sky=0.55
    fit=1
    Z=0

    xpos=367
    ypos=357
    anglegass=35.8

#    scale = 0.0455  # arcsec/pixel

    scale = 0.68  # arcsec/pixel set to default


    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

#    imgname = "ngc4342.fits"

    if imgname.find(".") != -1:
        (TNAM, trash) = imgname.split(".")
    else:
        TNAM=imgname

    exptime = GetExpTime(imgname)

##################
#################

    # Here we use an accurate four gaussians MGE PSF for
    # the HST/WFPC2/F814W filter, taken from Table 3 of
    # Cappellari et al. (2002, ApJ, 578, 787)

#    sigmapsf = [0.494, 1.44, 4.71, 13.4]      # In PC1 pixels
#    normpsf = [0.294, 0.559, 0.0813, 0.0657]  # total(normpsf)=1

#    ang=90-ang

######################
#    sky=back
#################

    hdu = fits.open(imgname)
    img = hdu[0].data

#    skylev = 0.55   # counts/pixel

    minlevel = 2  # counts/pixel
    ngauss = 12

    plt.clf()
#    f = find_galaxy(img, fraction=0.04, plot=1)

    eps,theta,xpeak,ypeak,back=StarSextractor(imgname,X,Y)
    print("star found at ",xpeak,ypeak)
    print("Ellipticity, Angle = ",eps,theta)

    xpos=xpeak
    ypos=ypeak

    if flagsky == True:
        img = img - sky   # subtract sky
        print("Sky = ",sky)

    else:
        print("Sky = ",back)

        sky=back
        img = img - back   # subtract sky

# I have to switch x and y coordinates, don't ask me why
    xtemp=xpeak
    xpeak=ypeak
    ypeak=xtemp

    plt.clf()

#    s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel, plot=1)
    s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel)


#            s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel, plot=1)

    plt.pause(2)  # Allow plot to appear on the screen

    ###########################

#            m = mge_fit_sectors(s.radius, s.angle, s.counts, eps,
#                                ngauss=ngauss, sigmapsf=sigmapsf, normpsf=normpsf,
#                                scale=scale, plot=1, bulge_disk=0, linear=0)

#    m = mge_fit_sectors(s.radius, s.angle, s.counts, eps,
#                        ngauss=ngauss, scale=scale, plot=1, bulge_disk=0, linear=0)

    m = mge_fit_sectors(s.radius, s.angle, s.counts, eps,
                        ngauss=ngauss, scale=scale, bulge_disk=0, linear=0)



    #    print(len(m.sol.T))
    #    print(len(m.sol))

    (counts,sigma,axisrat)=m.sol
    anglegass = 90 - theta

    hdu.close()

### print GALFIT files
####################
#####################

    parfile="psfgas.txt"
    outname=TNAM

    rmsname="none"
    psfname="none"

    consfile="constraints"

    T1 = "{}".format(imgname)
    T2 = outname + "-gas.fits"
    T3 = "{}".format(rmsname)

#    xlo=1
#    ylo=1

#    (xhi,yhi)=GetAxis(imgname)

## computing PSF size

    if flageven:

        xlo= xpos - boxsize/2
        ylo= ypos - boxsize/2

        xhi= xpos + boxsize/2
        yhi= ypos + boxsize/2

        xlo=np.int(np.round(xlo))
        ylo=np.int(np.round(ylo))
        xhi=np.int(np.round(xhi))
        yhi=np.int(np.round(yhi))

    else:

        xlo= xpeak - boxsize/2 + 0.5
        ylo= ypeak - boxsize/2 + 0.5

        xhi= xpeak + boxsize/2 + 0.5
        yhi= ypeak + boxsize/2 + 0.5

        xlo=np.int(np.round(xlo))
        ylo=np.int(np.round(ylo))
        xhi=np.int(np.round(xhi))
        yhi=np.int(np.round(yhi))


#    scale = 0.0455

    fout1 = open(parfile, "w")

    PrintHeader(fout1, T1, T2, T3, psfname, 1, maskfile, consfile, xlo, xhi, ylo,
                yhi, convbox, convbox, mgzpt, scale, scale, "regular", 0, 0)

    index = 0

#    index+=1

    while index < len(counts):

        TotCounts = counts[index]
        SigPix =  sigma[index]
        qobs =  axisrat[index]

        TotCounts=np.float(TotCounts)
        SigPix=np.float(SigPix)
        qobs=np.float(qobs)


        C0=TotCounts/(2*np.pi*qobs*SigPix**2)

        Ftot = 2*np.pi*SigPix**2*C0*qobs

        mgemag = mgzpt + 0.1  + 2.5*np.log10(exptime)  - 2.5*np.log10(Ftot)

        FWHM = 2.35482 * SigPix

    #        mgesb= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1


        outline = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f}  q: {:.2f} angle: {:.2f} \n".format(mgemag, SigPix, FWHM, qobs, anglegass)
        print(outline)

        PrintGauss(fout1, index, xpos, ypos, mgemag, FWHM, qobs, anglegass, Z, fit)


        index+=1


# sky removed as requested by GALFIT
    PrintSky(fout1, index, sky, 1, fit)
    fout1.close()

    print("Done. Run GALFIT as : 'galfit -o1 psfgas.txt' to create psf model ")
    print("or  Run GALFIT as : 'galfit psfgas.txt' to adjust the mge psf model ")



####################################################
####################################################
#####################################################


def StarSextractor(filimage,X,Y):

    sexfile="default.sex"
    CatSex="sim.cat"

    runcmd="sextractor -c {} {} ".format(sexfile,filimage)
    print(runcmd)

    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  # Run GALFIT

    KronScale=1
###############  Read in sextractor sorted data ################
    if os.stat(CatSex).st_size != 0 :
        Num, RA, Dec, XPos, YPos, Mag, Kron, FluxRad, IsoArea, AIm, E, Theta, Background, Class, Flag = np.genfromtxt(
           CatSex, delimiter="", unpack=True)  # sorted
## checking if they are numpy arrays:


#        index =  Mag.argsort()
#        i=0
        dist=10

        for idx, item in enumerate(Num):
#        for i in index:

            dx= XPos[idx] - X
            dy= YPos[idx] - Y
            dt= np.sqrt(dx**2 + dy**2)


            cflag=CheckFlag(4,Flag[idx])
            if cflag==False and dt < dist:
                dist = dt
                sindex=idx



        Num=Num[sindex]
        RA=RA[sindex]
        Dec=Dec[sindex]
        XPos=XPos[sindex]
        YPos=YPos[sindex]
        Mag=Mag[sindex]
        Kron=Kron[sindex]
        FluxRad=FluxRad[sindex]
        IsoArea=IsoArea[sindex]
        AIm=AIm[sindex]
        E=E[sindex]
        Theta=Theta[sindex]
        Background=Background[sindex]
        Class=Class[sindex]
        Flag=Flag[sindex]

##
        Angle = np.abs(Theta)
        AR = 1 - E
        RKron = KronScale * AIm * Kron
        Sky = Background


############

    XPos=np.round(XPos)
    YPos=np.round(YPos)

    XPos=np.int(XPos)
    YPos=np.int(YPos)


    print("Sextractor done")
    return E,Angle,XPos,YPos,Background






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
    main()
