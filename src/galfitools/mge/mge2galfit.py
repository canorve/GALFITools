#!/usr/bin/env python3

import os
import os.path
import subprocess as sp
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from galfitools.galin.galfit import Galfit
from galfitools.galin.MaskDs9 import GetAxis
from mgefit.mge_fit_sectors import mge_fit_sectors
from mgefit.mge_fit_sectors_regularized import mge_fit_sectors_regularized
from mgefit.mge_fit_sectors_twist import mge_fit_sectors_twist
from mgefit.sectors_photometry import sectors_photometry
from mgefit.sectors_photometry_twist import sectors_photometry_twist
from scipy.special import gammaincinv


def mge2gal(
    galfitFile,
    regfile,
    center,
    psf,
    twist,
    gauss,
    freeser,
    freesky,
    numgauss,
    xypos=None,
    ellip=None,
    posang=None,
) -> str:
    """Creates MGE initial parameters for GALFIT

    Creates a Multi-Gaussian Expansion (MGE) model and
    formats it into an initial parameter file for fitting by GALFIT.

    galfitFile : str
                GALFIT file from which the header information
                and sky value will be extracted to create the MGE model.
    regfile : str
             DS9 ellipse region file, which must enclose the galaxy to be fitted.
    center : Bool
            if True it will take the geometric's center of the DS9 ellipse
            as the center, otherwise it will take the pixel with the peak's value
    psf : float
          value of the PSF sigma
    twist : bool
            If True, the twist option for MGE is enabled, allowing the
            angular positions of the Gaussians to differ from one another.
    gauss : bool
            if True, it uses the gaussian model instead of the Sersic model
            with n = 0.5
    freeser : bool
            leaves the sersic index as a free parameter to fit
    freesky : bool
            leaves the sky parameter as a free parameter to fit
    numgauss : int
            maximum number of gaussians allowed to fit
    xypos : list, optional
            provides the (x y) position center of the object to fit
    ellip : float
            ellipticity of the object.
    posang : position angle of object. Measured from Y-axis


    Returns
    -------
    T2 : str
        name of the output FITS


    Notes
    -----
    The creation of initial parameters for MGE are created through
    the routine of Cappellari, MNRAS 333,400-410 (2002)


    """
    # reading options from galfit header
    galfit = Galfit(galfitFile)
    head = galfit.ReadHead()
    galsky = galfit.ReadSky()

    image = head.inputimage
    imageout = head.outimage
    magzpt = head.mgzpt
    maskfile = head.maskimage

    sky = galsky.sky

    scale = head.scale
    psfile = head.psfimage

    sigfile = head.sigimage

    convbox = head.convx
    convboxy = head.convy

    consfile = head.constraints
    xlo = head.xmin
    ylo = head.ymin
    xhi = head.xmax
    yhi = head.ymax

    #################

    initgauss = 12
    regu = None  # regu option now available

    mgeoutfile = "mgegas.txt"

    #  Updating variables for empty files ###########

    if (maskfile == "None") or (maskfile == "none"):
        maskfile = None

    if (sigfile == "None") or (sigfile == "none"):
        sigfile = None

    if (psfile == "None") or (psfile == "none"):
        psfile = None

    if regu:
        try:
            from mgefit.mge_fit_sectors_twist_regularized import (
                mge_fit_sectors_twist_regularized,
            )
        except Exception:
            print(
                "mge_fit_sectors_twist_regularized is not installed. Regu desactivated"
            )
            regu = False

    ##########################################
    ############################################

    if regu:
        print("regularization mode activated")

    root_ext = os.path.splitext(image)
    timg = root_ext[0]

    namepng = timg + ".png"

    sectorspng = "sectors.png"

    magzpt = float(magzpt)

    fit = 1
    serfit = 0
    skyfit = 0

    Z = 0

    ###################################################
    #  Create GALFIT input file header to compute sky #
    ###################################################

    #    image = "ngc4342.fits"

    root_ext = os.path.splitext(imageout)
    TNAM = root_ext[0]

    exptime = GetExpTime(image)

    ######################
    #    Mask file
    #################
    if maskfile:

        errmsg = "file {} does not exist".format(maskfile)
        assert os.path.isfile(maskfile), errmsg
        hdu = fits.open(maskfile)
        mask = hdu[0].data
        maskb = np.array(mask, dtype=bool)
        hdu.close()

    else:
        mask = np.array([])

    ######################
    #    sky=back
    ######################

    hdu = fits.open(image)
    img = hdu[0].data

    img = img.astype(float)

    minlevel = 0

    plt.clf()

    (ncol, nrow) = GetAxis(image)

    eps = 0
    theta = 0

    if regfile:
        if os.path.isfile(regfile):

            obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)

            xx, yy, Rkron, theta, eps = Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle)
            if center:
                print("center of ds9 ellipse region will be used")
                xpeak, ypeak = xpos, ypos
            else:
                (xmin, xmax, ymin, ymax) = GetSize(
                    xx, yy, Rkron, theta, eps, ncol, nrow
                )
                xpeak, ypeak = GetPmax(img, mask, xmin, xmax, ymin, ymax)

    if xypos:
        xpeak = xypos[0] - 1
        ypeak = xypos[1] - 1

    if ellip:
        eps = ellip

    if posang:
        theta = posang

    print("galaxy found at ", xpeak + 1, ypeak + 1)
    print("Ellipticity, Angle = ", eps, theta)

    print("Sky = ", sky)

    img = img - sky  # subtract sky

    # I have to switch x and y coordinates, don't ask me why
    xtemp = xpeak
    xpeak = ypeak
    ypeak = xtemp

    plt.clf()

    if twist:
        print("twist mode mge fit sectors is used")
        # Perform galaxy photometry
        if maskfile:
            s = sectors_photometry_twist(
                img, theta, xpeak, ypeak, minlevel=minlevel, badpixels=maskb, plot=1
            )
        else:
            s = sectors_photometry_twist(
                img, theta, xpeak, ypeak, minlevel=minlevel, plot=1
            )

        plt.savefig(sectorspng)

        plt.pause(1)  # Allow plot to appear on the screen

        plt.clf()

        tolpsf = 0.001

        if np.abs(psf) > tolpsf:
            if regu:
                m = mge_fit_sectors_twist_regularized(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    sigmapsf=psf,
                    scale=scale,
                    plot=1,
                )
            else:
                m = mge_fit_sectors_twist(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    sigmapsf=psf,
                    scale=scale,
                    plot=1,
                )

        else:
            print("No convolution")
            if regu:
                m = mge_fit_sectors_twist_regularized(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    scale=scale,
                    plot=1,
                )
            else:
                m = mge_fit_sectors_twist(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    scale=scale,
                    plot=1,
                )

        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(namepng)

    else:
        # elif twist == False:

        if maskfile:
            s = sectors_photometry(
                img,
                eps,
                theta,
                xpeak,
                ypeak,
                minlevel=minlevel,
                badpixels=maskb,
                plot=1,
            )
        else:

            s = sectors_photometry(
                img, eps, theta, xpeak, ypeak, minlevel=minlevel, plot=1
            )

        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(sectorspng)

        plt.clf()

        tolpsf = 0.001
        if np.abs(psf) > tolpsf:

            if regu:

                m = mge_fit_sectors_regularized(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    sigmapsf=psf,
                    scale=scale,
                    plot=1,
                    bulge_disk=0,
                    linear=0,
                )
            else:
                m = mge_fit_sectors(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    sigmapsf=psf,
                    scale=scale,
                    plot=1,
                    bulge_disk=0,
                    linear=0,
                )

        else:
            print("No convolution")

            if regu:
                m = mge_fit_sectors_regularized(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    scale=scale,
                    plot=1,
                    bulge_disk=0,
                    linear=0,
                )
            else:
                m = mge_fit_sectors(
                    s.radius,
                    s.angle,
                    s.counts,
                    eps,
                    ngauss=initgauss,
                    scale=scale,
                    plot=1,
                    bulge_disk=0,
                    linear=0,
                )

        plt.pause(1)  # Allow plot to appear on the screen

        plt.savefig(namepng)

    if twist:

        (counts, sigma, axisrat, pa) = m.sol

        theta2 = 270 - theta
        alpha1 = pa - theta2
        alphaf = alpha1 - 90

    else:
        # elif twist == False:
        (counts, sigma, axisrat) = m.sol
        # anglegass = 90 - theta
        anglegass = theta

    ####################
    #  print GALFIT files
    #####################

    if gauss:
        parfile = "mgeGALFIT.txt"
    else:
        parfile = "mseGALFIT.txt"

    outname = TNAM

    rmsname = sigfile
    if psfile:
        psfname = psfile
    else:
        psfname = "None"

    # switch back
    xtemp = xpeak
    xpeak = ypeak
    ypeak = xtemp

    T1 = "{}".format(image)
    T2 = outname + "-mge.fits"
    T3 = "{}".format(rmsname)

    # now this is read from galfit file

    # xlo=1
    # ylo=1
    # (xhi, yhi) = (ncol, nrow)

    fout1 = open(parfile, "w")
    fout2 = open(mgeoutfile, "w")

    outline2 = "# Mag Sig(pixels) FWHM(pixels) q angle \n"
    fout2.write(outline2)

    # sersic K for n = 0.5

    K = gammaincinv(1, 0.5)

    if psfile:
        (ncol, nrow) = GetAxis(psfname)
        convbox = ncol + 1
        convboxy = nrow + 1

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
        convboxy,
        magzpt,
        scale,
        scale,
        "regular",
        0,
        0,
    )

    if freeser:
        serfit = 1
    else:
        serfit = 0

    # removing the last components:
    totGauss = len(counts)

    if not (numgauss):
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
        SigPix = sigma[index]
        qobs = axisrat[index]

        TotCounts = float(TotCounts)
        SigPix = float(SigPix)
        qobs = float(qobs)

        if twist:
            anglegass = alphaf[index]  # - 90
            anglegass = float(anglegass)

        if gauss:

            C0 = TotCounts / (2 * np.pi * qobs * SigPix**2)

            Ftot = 2 * np.pi * SigPix**2 * C0 * qobs

            mgemag = magzpt + 0.1 + 2.5 * np.log10(exptime) - 2.5 * np.log10(Ftot)

            FWHM = 2.35482 * SigPix

            outlinea = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f} ".format(
                mgemag, SigPix, FWHM
            )
            outlineb = "q: {:.2f} angle: {:.2f} \n".format(qobs, anglegass)

            outline = outlinea + outlineb
            print(outline)

            outline2 = "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} \n".format(
                mgemag, SigPix, FWHM, qobs, anglegass
            )
            fout2.write(outline2)

            PrintGauss(
                fout1,
                index + 1,
                xpeak + 1,
                ypeak + 1,
                mgemag,
                FWHM,
                qobs,
                anglegass,
                Z,
                fit,
            )

        else:

            C0 = TotCounts / (2 * np.pi * qobs * SigPix**2)

            Ftot = 2 * np.pi * SigPix**2 * C0 * qobs

            mgemag = magzpt + 0.1 + 2.5 * np.log10(exptime) - 2.5 * np.log10(Ftot)

            FWHM = 2.35482 * SigPix

            h = np.sqrt(2) * SigPix

            Re = (K ** (0.5)) * h

            outlinea = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f} ".format(
                mgemag, SigPix, FWHM
            )
            outlineb = "Re: {:.2f}  q: {:.2f} angle: {:.2f} \n".format(
                Re, qobs, anglegass
            )
            outline = outlinea + outlineb

            print(outline)

            outline2 = "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} \n".format(
                mgemag, SigPix, FWHM, Re, qobs, anglegass
            )
            fout2.write(outline2)

            PrintSersic(
                fout1,
                index + 1,
                xpeak + 1,
                ypeak + 1,
                mgemag,
                Re,
                0.5,
                qobs,
                anglegass,
                Z,
                fit,
                serfit,
            )

        index += 1

    if freesky:
        skyfit = 1
    else:
        skyfit = 0

    PrintSky(fout1, index + 1, sky, Z, skyfit)
    fout1.close()
    fout2.close()

    makeConstraints(consfile, len(counts))

    print(
        "Done. Gaussians are stored in {}, and {} for galfit format ".format(
            mgeoutfile, parfile
        )
    )

    # returns output name
    return T2


####################################################
####################################################
#####################################################


def ReadMgePsf(psfile):
    """Reads the file created by mge2gal
    returns the normalized psf.

    Notes
    -----
    Not used in GALFITools and will be removed
    in the future

    """

    counts, mag, sigpix, fwhm, qobs, angle = np.genfromtxt(
        psfile, skip_header=1, delimiter="", unpack=True
    )

    tot = counts.sum()
    psfs = sigpix

    normpsf = counts / tot

    return psfs, normpsf


def Sextractor(filimage, X, Y):
    """
    Runs Sextractor and recover photometry from
    the object indicated by coordinates X, Y

    Notes
    -----
    Not used in the library and will be removed in the future

    """

    sexfile = "default.sex"
    CatSex = "sim.cat"

    runcmd = "sextractor -c {} {} ".format(sexfile, filimage)
    print(runcmd)

    sp.run(
        [runcmd], shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    )  # Run GALFIT

    # KronScale = 1
    #    Read in sextractor sorted data
    if os.stat(CatSex).st_size != 0:
        (
            Num,
            RA,
            Dec,
            XPos,
            YPos,
            Mag,
            Kron,
            FluxRad,
            IsoArea,
            AIm,
            E,
            Theta,
            Background,
            Class,
            Flag,
        ) = np.genfromtxt(
            CatSex, delimiter="", unpack=True
        )  # sorted
        # checking if they are numpy arrays:

        #        index =  Mag.argsort()
        #        i=0
        dist = 100

        for idx, item in enumerate(Num):

            #        if(isinstance(Num,np.ndarray)):
            dx = XPos[idx] - X
            dy = YPos[idx] - Y
            dt = np.sqrt(dx**2 + dy**2)

            cflag = CheckFlag(4, Flag[idx])
            if cflag is False and dt < dist:
                dist = dt
                sindex = idx

        Num = Num[sindex]
        RA = RA[sindex]
        Dec = Dec[sindex]
        XPos = XPos[sindex]
        YPos = YPos[sindex]
        Mag = Mag[sindex]
        Kron = Kron[sindex]
        FluxRad = FluxRad[sindex]
        IsoArea = IsoArea[sindex]
        AIm = AIm[sindex]
        E = E[sindex]
        Theta = Theta[sindex]
        Background = Background[sindex]
        Class = Class[sindex]
        Flag = Flag[sindex]

        Angle = np.abs(Theta)
        # AR = 1 - E
        # RKron = KronScale * AIm * Kron
        # Sky = Background

        XPos = round(XPos)
        YPos = round(YPos)

        XPos = int(XPos)
        YPos = int(YPos)

    else:
        print("Sextractor cat file not found ")

    print("Sextractor done")
    return E, Angle, XPos, YPos, Background


#  GALFIT functions


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
    """Print the GALFIT header parameters to a file specified
    by the filehandle.


    Parameters
    ----------

    hdl : str
        filehandler

    The rest of Parameters A, B, C, D, ..., S
    follows the logic order of the GALFIT file

    Returns
    -------
    bool

    """
    # k Check
    # the header for GALFIT

    lineZa = "============================================="
    lineZb = "=====================================================\n"

    lineZ = lineZa + lineZb

    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}    # Input Data image (FITS file)  \n".format(A)
    lineB = "B) {}    # Output data image block     \n".format(B)
    lineC = "C) {}    # Sigma image name   \n".format(C)
    lineD = "D) {}    # Input PSF image and     \n".format(D)
    lineE = "E) {}    # PSF fine sampling factor relative to data \n".format(E)
    lineF = "F) {}    # Bad pixel mask \n".format(F)
    lineG = "G) {}    # File with parameter constraints\n".format(G)
    lineH = "H) {} {} {} {}   # Image region to fit (xmin xmax ymin ymax) \n".format(
        xlo, xhi, ylo, yhi
    )
    lineI = "I) {} {}  # Size of the convolution box (x y)   \n".format(convx, convy)
    lineJ = "J) {}     # Magnitude photometric zeropoint   \n".format(J)
    lineK = "K) {} {}  # Plate scale (dx dy). [arcsec per pixel] \n".format(
        platedx, platedy
    )
    lineO = "O) {}  # Display type (regular, curses, both)   \n".format(Op)
    lineP = "P) {}  # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps \n".format(P)
    lineS = "S) {}  # Modify/create objects interactively?   \n".format(S)
    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS       \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:  \n"
    line3 = "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n"
    line4 = "#       ferrer, powsersic, sky, and isophote.       \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified:\n"
    line7 = "#      C0 (diskyness/boxyness),              \n"
    line8 = "#      Fn (n=integer, Azimuthal Fourier Modes),\n"
    line9 = "#      R0-R10 (PA rotation, for creating spiral structures).\n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number \n"
    line12 = "# column 2:         \n"
    line13 = "#   -- Parameter 0: the allowed functions are: sersic, nuker, expdisk\n"
    line14 = "#       edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky    \n"
    line15 = "#   -- Parameter 1-10: value of the initial parameters\n"
    line16 = "#   -- Parameter C0:   For diskiness/boxiness\n"
    line17 = "#      <0 = disky        \n"
    line18 = "#      >0 = boxy  \n"
    line19 = "#   -- Parameter Z:    Outputting image options, the options are: \n"
    line20 = "#      0 = normal, i.e. subtract final model from the data to create \n"
    line21 = "#      the residual image  \n"
    line22 = "#      1 = Leave in the model -- do not subtract from the data \n"
    line23 = "#      \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0) \n"
    line25 = "# column 4: comment    \n"
    line26 = " \n"

    line27a = "============================================="
    line27b = "=====================================================\n"

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
    """Print the GALFIT sky function parameters to a
    file specified by the filehandle.


    Parameters
    ----------

    hdl : str
        filehandler

    ncomp : int
        number of component
    sky : float
         value of sky background
    Z : int
        skip object. Z GALFIT component parameter
    fit : int
        leave parameter free/fixed during fit


    Returns
    -------
    bool

    """

    # k Check

    line00 = "# Object number: {}  \n".format(ncomp)
    line01 = " 0)   sky   #    Object type  \n"
    line02 = " 1) {}   {}   # sky background   [ADU counts] \n".format(sky, fit)
    line03 = " 2) 0.000      0     # dsky/dx (sky gradient in x) \n"
    line04 = " 3) 0.000      0     # dsky/dy (sky gradient in y) \n"
    line05 = " Z) {}    # Skip this model in output image?\n".format(Z)
    line06 = "\n"
    line07a = "=========================================="
    line07b = "======================================\n"
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


def PrintSersic(
    hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit, serfit
):
    """Print the GALFIT Sersic function parameters to a
    file specified by the filehandle.


    Parameters
    ----------

    hdl : str
        filehandler

    ncomp : int
        number of component

    xpos, ypos : int, int
                pixel position of the component's center

    magser : float
            magnitud of Sersic component
    reser : float
            effective radius in pixels
    nser : float
            Sersic index
    axratser : axis ratio of Sersic component

    angleser : angular position of component measured from Y-axis
               in degrees
    Z : int
        skip object. Z GALFIT component parameter
    fit : int
        leave parameter free/fixed during fit


    Returns
    -------
    bool

    """

    line00 = "# Object number: {}     \n".format(ncomp)
    line01 = " 0)     sersic      #  Object type \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  #  position x, y     [pixel]  \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}  {}    #  total magnitude    \n".format(magser, fit)
    line04 = " 4) {:.2f}  {}    #  R_e    [Pixels]   \n".format(reser, fit)
    line05 = " 5) {}     {}   #  Sersic exponent (deVauc=4, expdisk=1)  \n".format(
        nser, serfit
    )
    line06 = " 6)  0.0000       0   #  ---------------- \n"
    line07 = " 7)  0.0000       0   #  ---------------- \n"
    line08 = " 8)  0.0000       0   #  ---------------- \n"
    line09 = " 9) {:.2f}  {}  #  axis ratio (b/a)     \n".format(axratser, fit)
    line10 = "10) {:.2f}  {}  #  position angle (PA)     \n".format(angleser, fit)
    lineZ = " Z) {}    #  Skip this model in output image? \n".format(Z)
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
    """Print the GALFIT Gauss function parameters to a
    file specified by the filehandle.


    Parameters
    ----------

    hdl : str
        filehandler

    ncomp : int
        number of component

    xpos, ypos : int, int
                pixel position of the component's center
    magass : float
            magnitud of Gaussian component
    fwhm : float
            Full Width Half Maximum in pixels
    axratgass : axis ratio of Gauss component

    anglegass : angular position of component measured from Y-axis
               in degrees
    Z : int
        skip object. Z GALFIT component parameter
    fit : int
        leave parameter free/fixed during fit


    Returns
    -------
    bool

    """

    line00 = "# Object number: {}       \n".format(ncomp)
    line01 = " 0)     gaussian   #  Object type \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  #  position x, y     [pixel] \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}    {}   #  total magnitude   \n".format(magass, fit)
    line04 = " 4) {:.2f}    {}   #  FWHM  [Pixels] \n".format(fwhm, fit)
    line05 = " 9) {:.2f}    {}   #  axis ratio (b/a) \n".format(axratgass, fit)
    line06 = "10) {:.2f}  {}  #  position angle (PA) [Degrees: Up=0, Left=90]\n".format(
        anglegass, fit
    )
    line07 = " Z) {}    #  Skip this model in output image?  (yes=1, no=0) \n".format(Z)
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
    """Print the GALFIT exponential function parameters to a
    file specified by the filehandle.


    print GALFIT exponential function to filehandle
    Parameters
    ----------

    hdl : str
        filehandler

    ncomp : int
        number of component

    xpos, ypos : int, int
                pixel position of the component's center
    magexp : float
            magnitud of exponential component
    rsexp : float
            scale radius of exponential component
    axratexp : axis ratio of exponential component
    angleexp : angular position of component measured from Y-axis
               in degrees
    Z : int
        skip object. Z GALFIT component parameter
    fit : int
        leave parameter free/fixed during fit

    Returns
    -------
    bool

    """

    line00 = "# Object number: $ncomp      \n"
    line01 = " 0)     expdisk        # Object type \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  # position x, y [pixel]\n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}    {}   # total magnitude \n".format(magexp, fit)
    line04 = " 4) {:.2f}    {}   #      Rs  [Pixels]\n".format(rsexp, fit)
    line05 = " 9) {:.2f}        {}   # axis ratio (b/a) \n".format(axratexp, fit)
    line06 = "10) {:.2f}  {} # position angle (PA)   \n".format(angleexp, fit)
    line07 = " Z) {}      # Skip this model in output image?  (yes=1, no=0) \n".format(
        Z
    )
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
    """Get a cutout of the image
    # repeated

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
    """Get exposition time from the FITS image header
    # repeated

    """

    try:
        hdu = fits.open(Image)
        exptime = hdu[0].header["EXPTIME"]
        hdu.close()
    except Exception:
        exptime = 1
    return float(exptime)


def makeConstraints(consfile: str, numcomp: int) -> True:
    """Creates a contraints file for GALFIT with the number of components

    Parameters
    ----------
    consfile : str
               the new constraints file
    numcomp : int
             the number of components to constraint

    Returns
    -------
    bool

    """

    fout = open(consfile, "w")

    line00 = "# Component/    parameter   constraint	Comment           \n"
    line01 = "# operation	(see below)   range         \n"

    linempty = "                  \n"

    comp = np.arange(2, numcomp + 1)

    cad = "1"
    for idx, item in enumerate(comp):
        cad = cad + "_" + str(item)

    line02 = "    " + cad + "    " + "x" + "    " + "offset \n"
    line03 = "    " + cad + "    " + "y" + "    " + "offset \n"

    fout.write(line00)
    fout.write(line01)
    fout.write(linempty)
    fout.write(line02)
    fout.write(line03)

    fout.close()

    return True


def CheckFlag(val, check):
    """Check if flag is contained in val,
    returns True if found.

    This is useful for Sextractor Flags


    # repeated

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


def GetInfoEllip(regfile):
    """Extracts parameters' information
    from DS9 ellipse region.

    # repeated
    """

    if not os.path.exists(regfile):
        print("%s: reg filename does not exist!" % (regfile))
        sys.exit()

    f1 = open(regfile, "r")

    lines = f1.readlines()

    f1.close()

    flag = False
    found = False

    # reading reg file
    for line in lines:

        line = line.split("#")
        line = line[0]

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True
            found = True

        if flag is True:
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
        obj = v0
        xpos = v1
        ypos = v2
        rx = v3
        ry = v4
        angle = v5

        # avoids ds9 regions with Area = 0
        if rx < 1:
            rx = 1
        if ry < 1:
            ry = 1

        return obj, xpos, ypos, rx, ry, angle
        # return eps, theta, xpos, ypos
    else:
        print("ellipse region was not found in file. Exiting.. ")
        sys.exit()

    return 0, 0, 0, 0


def Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle):
    """converts DS9 ellipse parameters to geometrical
    parameters

    # repeated

    """

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle - 90
        xx = xpos
        yy = ypos
    else:
        q = rx / ry
        e = 1 - q
        Rkron = ry
        theta = angle  # + 90
        xx = xpos
        yy = ypos

    return xx, yy, Rkron, theta, e


def GetSize(x, y, R, theta, ell, ncol, nrow):
    """Gets the maximum and minimum pixels
    positions for an ellipse

    # repeated

    """
    # k Check
    q = 1 - ell
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    xmin = x - np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )

    xmax = x + np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )

    ymin = y - np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
    )

    ymax = y + np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
    )

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

    return (xmin, xmax, ymin, ymax)


def GetPmax(image, mask, xmin, xmax, ymin, ymax):
    """Obtains the pixel position with the highest value within
    the indicated region.

    # repeated
    """

    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    chuckimg = image[ymin - 1 : ymax, xmin - 1 : xmax]
    if mask.any():
        chuckmsk = mask[ymin - 1 : ymax, xmin - 1 : xmax]

        invmask = np.logical_not(chuckmsk)

        invmask = invmask * 1

        chuckimg = chuckimg * invmask

    maxy, maxx = np.where(chuckimg == np.max(chuckimg))

    xpos = maxx[0] + xmin - 1
    ypos = maxy[0] + ymin - 1

    return (xpos, ypos)


#############################################################################
#  End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
