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
from scipy.special import gammaincinv

import mgefit
from mgefit.find_galaxy import find_galaxy
from mgefit.sectors_photometry import sectors_photometry
from mgefit.sectors_photometry_twist import sectors_photometry_twist


from mgefit.mge_fit_sectors import mge_fit_sectors
from mgefit.mge_fit_sectors_twist import mge_fit_sectors_twist

from mgefit.mge_fit_sectors_regularized import mge_fit_sectors_regularized 

import argparse



def mge2gal(args) -> None:


    image = args.image
    regfile = args.Ds9regFile 
    magzpt = args.magzpt
    twist = args.twist
    regu = args.regu
    center = args.center
    maskfile = args.mask

    psf = args.psf
    sky = args.sky
    scale = args.plate
    gauss = args.gauss

    psfile = args.psfile
    sigfile = args.sigfile

    freeser = args.freeser
    freesky = args.freesky

    regu = args.regu


    numgauss = args.numgauss 

    initgauss = 12



    mgeoutfile="mgegas.txt"

    #default convolution box sizes:
    convbox=100
    convboxy=100
    

    if regu:
        try:
            from mgefit.mge_fit_sectors_twist_regularized import mge_fit_sectors_twist_regularized 
        except: 
            print("mge_fit_sectors_twist_regularized is not installed. Regu desactivated")
            regu = False

    ##########################################
    ############################################


    if regu:
        print('regularization mode activated')


    root_ext = os.path.splitext(image)
    timg= root_ext[0]


    namepng=timg+".png"

    sectorspng="sectors.png"

    magzpt=float(magzpt)



    fit=1
    serfit=0
    skyfit=0

    Z=0


    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    #    image = "ngc4342.fits"

    root_ext = os.path.splitext(image)
    TNAM= root_ext[0]




    exptime = GetExpTime(image)



    ######################
    #    Mask file
    #################


    if maskfile:

        errmsg="file {} does not exist".format(maskfile)
        assert os.path.isfile(maskfile), errmsg
        hdu = fits.open(maskfile)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()

    else:
        mask=np.array([])

    ######################
    #    sky=back
    ######################

    hdu = fits.open(image)
    img = hdu[0].data

    img=img.astype(float)

    minlevel = 0  

    plt.clf()

    (ncol, nrow) = GetAxis(image)
    
    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)
    xx, yy, Rkron, theta, eps = Ds9ell2Kronell(xpos,ypos,rx,ry,angle)
    if center:
        print('center of ds9 ellipse region will be used')
        xpeak, ypeak = xpos, ypos
    else:        
        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, eps, ncol, nrow)
        xpeak, ypeak = GetPmax(img, mask, xmin, xmax, ymin, ymax)
 
    print("galaxy found at ", xpeak + 1, ypeak + 1)
    print("Ellipticity, Angle = ", eps, theta)

    print("Sky = ", sky)

    img = img - sky   # subtract sky

    # I have to switch x and y coordinates, don't ask me why
    xtemp=xpeak
    xpeak=ypeak
    ypeak=xtemp

    plt.clf()

    if twist:
        print("twist mode mge fit sectors is used")
                        # Perform galaxy photometry
        if maskfile:
            s = sectors_photometry_twist(img, theta, xpeak, ypeak, 
                                         minlevel=minlevel, badpixels=maskb ,plot=1)
        else:
            s = sectors_photometry_twist(img, theta, xpeak, ypeak, minlevel=minlevel, 
                                         plot=1)

        plt.savefig(sectorspng)

        plt.pause(1)  # Allow plot to appear on the screen

        plt.clf()

        tolpsf = 0.001

        if np.abs(psf) > tolpsf:
            if regu:
                m = mge_fit_sectors_twist_regularized(s.radius, s.angle, s.counts, eps, ngauss=initgauss,
                                            sigmapsf=psf, scale=scale, plot=1)
            else:
                m = mge_fit_sectors_twist(s.radius, s.angle, s.counts, eps, ngauss=initgauss,
                                            sigmapsf=psf, scale=scale, plot=1)
 
        else:
            print("No convolution")
            if regu:
                m = mge_fit_sectors_twist_regularized(s.radius, s.angle, s.counts, eps, ngauss=initgauss,
                                             scale=scale, plot=1)
            else:
                m = mge_fit_sectors_twist(s.radius, s.angle, s.counts, eps, ngauss=initgauss,
                                             scale=scale, plot=1)
 


        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(namepng)


    elif twist == False:

        if maskfile:
            s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel, badpixels=maskb ,plot=1)
        else:

            s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel, plot=1)

#            s = sectors_photometry(img, eps, theta, xpeak, ypeak, minlevel=minlevel, plot=1)

        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(sectorspng)

    ###########################

        plt.clf()



        tolpsf = 0.001
        if np.abs(psf) > tolpsf:

            if regu:
 
                m = mge_fit_sectors_regularized(s.radius, s.angle, s.counts, eps,
                                    ngauss=initgauss, sigmapsf=psf,
                                    scale=scale, plot=1, bulge_disk=0, linear=0)
            else:
                m = mge_fit_sectors(s.radius, s.angle, s.counts, eps,
                                    ngauss=initgauss, sigmapsf=psf,
                                    scale=scale, plot=1, bulge_disk=0, linear=0)
 

        else:
            print("No convolution")

            if regu:
                m = mge_fit_sectors_regularized(s.radius, s.angle, s.counts, eps,
                                    ngauss=initgauss, scale=scale, plot=1, 
                                    bulge_disk=0, linear=0)
            else:
                m = mge_fit_sectors(s.radius, s.angle, s.counts, eps,
                                    ngauss=initgauss, scale=scale, plot=1, 
                                    bulge_disk=0, linear=0)
 

        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(namepng)


    if twist:

        (counts, sigma, axisrat, pa) = m.sol

        theta2 = 270 - theta
        alpha1 = pa - theta2
        alphaf = alpha1 - 90

    elif twist == False:
        (counts, sigma, axisrat) = m.sol
        #anglegass = 90 - theta
        anglegass =  theta


    ####################
    ### print GALFIT files
    #####################

    if gauss:
        parfile="mgeGALFIT.txt"
    else:
        parfile="mseGALFIT.txt"

    outname=TNAM

    rmsname = sigfile 
    if psfile:
        psfname = psfile
    else:
        psfname = "psf.fits"


    # switch back
    xtemp=xpeak
    xpeak=ypeak
    ypeak=xtemp


    consfile="constraints.txt"

    T1 = "{}".format(image)
    T2 = outname + "-mge.fits"
    T3 = "{}".format(rmsname)

    xlo=1
    ylo=1

    (xhi, yhi) = (ncol, nrow)


    fout1 = open(parfile, "w")
    fout2 = open(mgeoutfile, "w")


    outline2 = "# Mag Sig(pixels) FWHM(pixels) q angle \n"
    fout2.write(outline2)

    #sersic K for n = 0.5

    K = gammaincinv(1,0.5)


    if psfile:
        (ncol, nrow) = GetAxis(psfname)
        convbox = ncol + 1
        convboxy = nrow + 1



    PrintHeader(fout1, T1, T2, T3, psfname, 1, maskfile, consfile, xlo, xhi, ylo,
                yhi, convbox, convboxy, magzpt, scale, scale, "regular", 0, 0)


    if freeser:
        serfit = 1
    else:
        serfit = 0


    #removing the last components:
    totGauss = len(counts)

    if not(numgauss):
        numgauss = totGauss

    print("total number of gaussians by mge_fit_sectors: ", totGauss)

    if numgauss > totGauss:
        print("number of gaussians is greater than the ones fitted by mge_fit_sectors")
        print("total number of gaussians of mge_fit_sectors will be used instead")
        numgauss = totGauss

    print("total number of gaussians used by GALFIT: ", numgauss)



    index = 0

    while index < numgauss: 


        TotCounts = counts[index]
        SigPix =  sigma[index]
        qobs =  axisrat[index]

        TotCounts = float(TotCounts)
        SigPix = float(SigPix)
        qobs = float(qobs)

        if twist:
            anglegass = alphaf[index] #- 90
            anglegass = float(anglegass)


        if gauss:
            
            C0=TotCounts/(2*np.pi*qobs*SigPix**2)

            Ftot = 2*np.pi*SigPix**2*C0*qobs

            mgemag = magzpt + 0.1  + 2.5*np.log10(exptime)  - 2.5*np.log10(Ftot)

            FWHM = 2.35482*SigPix

            outline = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f}  q: {:.2f} angle: {:.2f} \n".format(mgemag, SigPix, FWHM, qobs, anglegass)
            print(outline)

            outline2 = "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} \n".format(mgemag, SigPix, FWHM, qobs, anglegass)
            fout2.write(outline2)


            PrintGauss(fout1, index+1, xpeak + 1, ypeak + 1, mgemag, FWHM, qobs, anglegass, Z, fit)


        else:

            C0=TotCounts/(2*np.pi*qobs*SigPix**2)

            Ftot = 2*np.pi*SigPix**2*C0*qobs

            mgemag = magzpt + 0.1  + 2.5*np.log10(exptime)  - 2.5*np.log10(Ftot)
            
            FWHM = 2.35482*SigPix

            h  = np.sqrt(2) * SigPix
             
            Re = (K**(0.5))*h

            outline = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f} Re: {:.2f}  q: {:.2f} angle: {:.2f} \n".format(mgemag, SigPix, FWHM, Re, qobs, anglegass)
            print(outline)

            outline2 = "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} \n".format(mgemag, SigPix, FWHM, Re,  qobs, anglegass)
            fout2.write(outline2)


            PrintSersic(fout1, index+1, xpeak + 1, ypeak + 1, mgemag, Re, 0.5, qobs, anglegass, Z, fit, serfit)



        index+=1


    if freesky:
        skyfit = 1
    else:
        skyfit = 0

    PrintSky(fout1, index+1, sky, Z, skyfit)
    fout1.close()
    fout2.close()

    makeConstraints(consfile, len(counts))

    print("Done. Gaussians are stored in {}, and {} for galfit format ".format(mgeoutfile,parfile))

####################################################
####################################################
#####################################################

def ReadMgePsf(psfile):


    counts,mag,sigpix,fwhm,qobs,angle = np.genfromtxt(psfile, skip_header=1, delimiter="", unpack=True)

    tot=counts.sum()
    psfs=sigpix

    normpsf=counts/tot

    return psfs,normpsf


def Sextractor(filimage,X,Y):

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
        dist=100

        for idx, item in enumerate(Num):


#        if(isinstance(Num,np.ndarray)):
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

        XPos=round(XPos)
        YPos=round(YPos)

        XPos=int(XPos)
        YPos=int(YPos)

    else:
        print("Sextractor cat file not found ")

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


def PrintSersic(hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit, serfit):
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
        nser, serfit)
    #line05 = " 5) {}       {}              #  Sersic exponent (deVauc=4, expdisk=1)           \n".format(
    #    nser, fit)
    line06 = " 6)  0.0000       0           #  ----------------                                \n"
    line07 = " 7)  0.0000       0           #  ----------------                                \n"
    line08 = " 8)  0.0000       0           #  ----------------                                \n"
    line09 = " 9) {:.2f}       {}              #  axis ratio (b/a)                                \n".format(
        axratser, fit)
    line10 = "10) {:.2f}       {}              #  position angle (PA)  [Degrees: Up=0, Left=90]   \n".format(
        angleser, fit)
    lineZ = " Z) {}                       #  Skip this model in output image?  (yes=1, no=0) \n".format(
        Z)
    line11 = "\n"

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
    hdl.write(line10)
    hdl.write(lineZ)
    hdl.write(line11)

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

    try:
        hdu = fits.open(Image)
        exptime = hdu[0].header["EXPTIME"]
        hdu.close()
    except: 
        exptime = 1
    return float(exptime)

def makeConstraints(consfile: str, numcomp: int) -> True:
    "creates a contraints file with the number of components"


    fout = open(consfile, "w")


    line00 = "#Component/    parameter   constraint	Comment           \n"
    line01 = "# operation	(see below)   range         \n"

    linempty = "                  \n"				

    comp = np.arange(2,numcomp+1)

    cad = "1"
    for idx, item in enumerate(comp):
        cad = cad + "_" + str(item)

    line02 = "    " + cad + "    " + "x" + "    "  + "offset \n"				
    line03 = "    " + cad + "    " + "y" + "    "  + "offset \n"				


    fout.write(line00)
    fout.write(line01)
    fout.write(linempty)
    fout.write(line02)
    fout.write(line03)

    fout.close()



    return True





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


def GetInfoEllip(regfile):

    if not os.path.exists(regfile):
        print ('%s: reg filename does not exist!' %(regfile))
        sys.exit()

    f1 = open(regfile,'r')

    lines = f1.readlines()

    f1.close()

    flag = False
    found = False

    #reading reg file
    for line in lines:


        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True
            found = True

        if (flag == True):
            x3 = p[4]
            x4 = x3[:-2]


            v0 = x0

            v1 = float(x2)
            v2 = float(p[1])
            v3 = float(p[2])
            v4 = float(p[3])
            v5 = float(x4)

            flag = False

    if found:
        obj  = v0
        xpos = v1
        ypos = v2
        rx = v3
        ry = v4
        angle = v5

        if rx >= ry:
            axratio = ry/rx 
            eps = 1 - axratio 
            theta = angle + 90 
        else:
            axratio = rx/ry 
            eps = 1 - axratio 
            theta = angle 

        return obj, xpos, ypos, rx, ry, angle
        #return eps, theta, xpos, ypos
    else:
        print("ellipse region was not found in file. Exiting.. ")
        sys.exit()

    return 0, 0, 0, 0



##new: 

#Image = MakeEllip(Image,Value,xpos[idx],ypos[idx],rx[idx],ry[idx],angle[idx],ncol,nrow)

#def MakeEllip(Image,Value,xpos,ypos,rx,ry,angle,ncol,nrow):
#    "Make an ellipse in an image"

#    xx, yy, Rkron, theta, e = Ds9ell2Kronell(xpos,ypos,rx,ry,angle)
#    (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, e, ncol, nrow)
#    Image = MakeKron(Image, Value, xx, yy, Rkron, theta, e, xmin, xmax, ymin, ymax)

#    return Image


def Ds9ell2Kronell(xpos,ypos,rx,ry,angle):


    if rx >= ry:

        q = ry/rx
        e = 1 - q
        Rkron = rx
        theta = angle + 90
        xx = xpos
        yy = ypos
    else:
        q = rx/ry
        e = 1 - q
        Rkron = ry
        theta = angle# + 90
        xx = xpos
        yy = ypos
 
     
    return xx, yy, Rkron, theta, e 


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



def GetPmax(image, mask, xmin, xmax, ymin, ymax):


    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    chuckimg = image[ymin - 1:ymax, xmin - 1:xmax]
    if mask.any():
        chuckmsk = mask[ymin - 1:ymax, xmin - 1:xmax]

        invmask = np.logical_not(chuckmsk)

        invmask = invmask*1

        chuckimg = chuckimg*invmask

    maxy, maxx = np.where(chuckimg == np.max(chuckimg))
    
    xpos = maxx[0] + xmin - 1
    ypos = maxy[0] + ymin - 1

    return (xpos, ypos)





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
