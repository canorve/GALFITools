#!/usr/bin/env python3

import os
import os.path
import subprocess as sp

import numpy as np
from astropy.io import fits
from galfitools.galin.MaskDs9 import GetAxis
from galfitools.galin.std import ds9satbox

# change it and use skyRem instead


def galfitSky(imgname, maskfile, mgzpt, scale, X, Y, sky) -> None:
    """Computes the sky using GALFIT

    It uses the FITS image and a mask to compute the
    to make an initial parameter file to compute sky background
    by GALFIT

    Parameters
    ----------
    imgname : str
            name of the FITS image

    maskfile : str
            name of the mask image
    mgzpt: float
            magnitud zeropoint
    scale : float
            plate scale
    X, Y : float, float x,y position of center
    sky : float
        initial parameter of sky

    Notes
    -----
    The function returns a GALFIT initial parameter file
    named sky.txt. The user must run GALFIT with this file
    to obtain the sky background


    """
    flagpos = False

    (NCol, NRow) = GetAxis(imgname)

    convbox = 100

    #    sky=0.55
    fit = 1
    Z = 0

    #    scale = 0.0455  # arcsec/pixel

    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    #    imgname = "ngc4342.fits"
    root_ext = os.path.splitext(imgname)
    TNAM = root_ext[0]

    ##################
    #################

    # Here we use an accurate four gaussians MGE PSF for
    # the HST/WFPC2/F814W filter, taken from Table 3 of
    # Cappellari et al. (2002, ApJ, 578, 787)

    #    sigmapsf = [0.494, 1.44, 4.71, 13.4]      # In PC1 pixels
    #    normpsf = [0.294, 0.559, 0.0813, 0.0657]  # total(normpsf)=1

    parfile = "sky.txt"
    outname = TNAM

    rmsname = "none"
    psfname = "none"

    consfile = "none"

    T1 = "{}".format(imgname)
    T2 = outname + "-sky.fits"
    T3 = "{}".format(rmsname)

    if flagpos is False:
        xlo = 1
        ylo = 1

        (xhi, yhi) = GetAxis(imgname)

    plate = 1

    fout1 = open(parfile, "w")

    PrintHeader(
        fout1,
        T1,
        T2,
        T3,
        psfname,
        1,
        maskfile,
        consfile,
        xlo,
        xhi,
        ylo,
        yhi,
        convbox,
        convbox,
        mgzpt,
        plate,
        plate,
        "regular",
        0,
        0,
    )

    index = 0

    index += 1

    PrintSky(fout1, index, sky, Z, fit)

    fout1.close()

    print("To compute sky run: galfit sky.txt")


####################################################

# GALFIT functions


def PrintHeader(
    hdl,
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    xlo,
    xhi,
    ylo,
    yhi,
    convx,
    convy,
    J,
    platedx,
    platedy,
    Op,
    P,
    S,
):
    """prints GALFIT header in a file

    # repeated

    """

    # k Check
    # print to filehandle
    # the header for GALFIT

    lineZa = "=========================================="
    lineZb = "========================================================\n"
    lineZ = lineZa + lineZb
    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}       # Input Data image (FITS file) \n".format(A)
    lineB = "B) {}   # Output data image block  \n".format(B)
    lineC = 'C) {}   # Sigma image name (made from data if blank or "none")  \n'.format(
        C
    )
    lineD = "D) {}   # Input PSF image and (optional) diffusion kernel\n".format(D)
    lineE = "E) {}   # PSF fine sampling factor relative to data  \n".format(E)
    lineF = "F) {}     # Bad pixel mask (FITS image or ASCII coord list) \n".format(F)
    lineG = "G) {}     # File with parameter constraints (ASCII file)  \n".format(G)
    lineH = "H) {} {} {} {}    # Image region to fit (xmin xmax ymin ymax)\n".format(
        xlo, xhi, ylo, yhi
    )
    lineI = "I) {} {}    # Size of the convolution box (x y)  \n".format(convx, convy)
    lineJ = "J) {}       # Magnitude photometric zeropoint  \n".format(J)
    lineK = "K) {} {}    # Plate scale (dx dy). [arcsec per pixel] \n".format(
        platedx, platedy
    )
    lineO = "O) {}    # Display type (regular, curses, both)   \n".format(Op)
    lineP = "P) {}    # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps\n".format(P)
    lineS = "S) {}    # Modify/create objects interactively?  \n".format(S)
    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS          \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:\n"
    line3 = "#    nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n"
    line4 = "#    ferrer, powsersic, sky, and isophote.   \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified:\n"
    line7 = "#  C0 (diskyness/boxyness),   \n"
    line8 = "#  Fn (n=integer, Azimuthal Fourier Modes),  \n"
    line9 = "#  R0-R10 (PA rotation, for creating spiral structures).  \n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number\n"
    line12 = "# column 2:     \n"
    line13 = "#  -- Parameter 0: the allowed functions are: sersic, nuker, expdisk\n"
    line14 = "#        edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky \n"
    line15 = "#        -- Parameter 1-10: value of the initial parameters \n"
    line16 = "#        -- Parameter C0:   For diskiness/boxiness \n"
    line17 = "#              <0 = disky   \n"
    line18 = "#              >0 = boxy    \n"
    line19 = "#       -- Parameter Z: Outputting image options, the options are:\n"
    line20 = "#  0 = normal, i.e. subtract final model from the data to create \n"
    line21 = "#            the residual image  \n"
    line22 = "#            1 = Leave in the model -- do not subtract from the data\n"
    line23 = "#            \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0)\n"
    line25 = "# column 4: comment  \n"
    line26 = " \n"

    line27a = "========================================="
    line27b = "=========================================================\n"
    line27 = line27a + line27b

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
    """Print GALFIT sky function to filehandle

    # repeated
    """

    # k Check

    line00 = "# Object number: {}   \n".format(ncomp)
    line01 = " 0)      sky       #    Object type  \n"
    line02 = " 1) {}     {}  # sky background  \n".format(sky, fit)
    line03 = " 2) 0.000      0    # dsky/dx (sky gradient in x) \n"
    line04 = " 3) 0.000      0     # dsky/dy (sky gradient in y)\n"
    line05 = " Z) {}        # Skip this model in output image?  \n".format(Z)
    line06 = "\n"
    line07a = "================================="
    line07b = "===============================================\n"
    line07 = line07a + line07b

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)

    return True


'''
def PrintSersic(
    hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit
):
    """print GALFIT Sersic function to filehandle

    # not used
    # repeated

    """
    # k Check

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}  \n".format(ncomp)
    line01 = " 0)     sersic      #  Object type    \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  #  position x, y   \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}  #  total magnitude     \n".format(magser, fit)
    line04 = " 4) {:.2f}  {}   #  R_e  [Pixels]\n".format(reser, fit)
    line05 = " 5) {}  {}    #  Sersic exponent (deVauc=4, expdisk=1)\n".format(
        nser, fit
    )
    line06 = " 9) {:.2f}  {}  #  axis ratio (b/a)  \n".format(axratser, fit)
    line07 = "10) {:.2f}   {}  #  position angle (PA)  \n".format(angleser, fit)
    line08 = " Z) {}       #  Skip this model in output image? \n".format(Z)
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
    """print GALFIT GAUSS function to filehandle

    # not used
    # repeated

    """

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}  \n".format(ncomp)
    line01 = " 0)     gaussian     #  Object type  \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {} #  position x, y \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}  #  total magnitude  \n".format(magass, fit)
    line04 = " 4) {:.2f}   {}  #  FWHM  [Pixels]  \n".format(fwhm, fit)
    line05 = " 9) {:.2f}   {}  #  axis ratio (b/a)  \n".format(axratgass, fit)
    line06 = "10) {:.2f}   {}   #  position angle (PA)  \n".format(anglegass, fit)
    line07 = " Z) {}    #  Skip this model in output image? \n".format(Z)
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
    """print GALFIT exponential function to filehandle

    # repeated
    # not used

    """

    # k Check

    # print to filehandle
    # a exponential function given
    # by the parameters

    line00 = "# Object number: $ncomp    \n"
    line01 = " 0)     expdisk     # Object type       \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  # position x, y  [pixel] \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}   # total magnitude \n".format(magexp, fit)
    line04 = " 4) {:.2f}   {} # Rs  [Pixels]  \n".format(rsexp, fit)
    line05 = " 9) {:.2f}   {} # axis ratio (b/a) \n".format(axratexp, fit)
    line06 = "10) {:.2f}   {}  # position angle (PA)  \n".format(angleexp, fit)
    line07 = " Z) {} # Skip this model in output image? \n".format(Z)
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
    """Get a piece from the FITS image

    # repeated
    # not used

    """
    # k Check

    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        sp.run(
            [runcmd],
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True,
        )

    hdu = fits.open(Image)
    dat = hdu[0].data[ylo - 1 : yhi, xlo - 1 : xhi]
    hdu[0].data = dat
    hdu.writeto(Imageout, overwrite=True)
    hdu.close()


def GetExpTime(Image):
    """Get exposition time from the image header

    # repeated
    # not used


    """

    hdu = fits.open(Image)
    exptime = hdu[0].header["EXPTIME"]
    hdu.close()
    return exptime


def CheckFlag(val, check):
    """Check for flag contained in val, returns True if found

    # repeated
    # used in function not used here

    """

    flag = False
    mod = 1
    maxx = 128

    while mod != 0:

        res = int(val / maxx)

        if maxx == check and res == 1:

            flag = True

        mod = val % maxx

        val = mod
        maxx = maxx / 2

    return flag


'''

#############################################################################
#    End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
