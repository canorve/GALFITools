#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path

import mgefit
from mgefit.mge_fit_1d import mge_fit_1d
from scipy.special import gammaincinv

#check modify
def main():



    rad,me=np.genfromtxt("Blackout.txt",unpack=True)
    #rad,me=np.genfromtxt("Blueout.txt",unpack=True)
    #rad,me=np.genfromtxt("Redout.txt",unpack=True)
    #rad,me=np.genfromtxt("Greenout.txt",unpack=True)
    #rad,me=np.genfromtxt("Magentaout.txt",unpack=True)




    lrad= np.log10(rad)
    lme= np.log10(me)
     


    print("\nFitting 1-dim profile-----------------------------------\n")
    counts,sigma = fit_1d(rad,me)

    print('mge done formatting to output file')


    parfile="mse1dGALFIT.txt"

    # printing files

    mgeoutfile="mgegas.txt"


    fout1 = open(parfile, "w")
    fout2 = open(mgeoutfile, "w")


    outline2 = "# Mag Sig(pixels) FWHM(pixels) Re(pixels) q angle \n"
    fout2.write(outline2)


    totGauss = len(counts)

    #param values

    index = 0
    magzpt = 25
    exptime = 1
    anglegass = 0

    fit = 1
    serfit = 0
    skyfit = 0
    sky = 0

    Z=0

    image = "galaxy.fits"
    rmsname = "galaxy-rms.fits"
    psfname = "psf.fits"
    maskfile="mask.fits"
    outname = "galaxy"


    consfile="constraints.txt"

    T1 = "{}".format(image)
    T2 = outname + "-mge.fits"
    T3 = "{}".format(rmsname)

    xlo = 1
    ylo = 1

    xhi = 2000
    yhi = 2000 

    xpeak = 1000
    ypeak = 1000


    convbox = 100
    scale = 1

    K = gammaincinv(1,0.5)

    PrintHeader(fout1, T1, T2, T3, psfname, 1, maskfile, consfile, xlo, xhi, ylo,
                yhi, convbox, convbox, magzpt, scale, scale, "regular", 0, 0)




    print("total number of gaussians by mge_fit_sectors: ", totGauss)

    while index < totGauss: 


        TotCounts = counts[index]
        SigPix =  sigma[index]
        qobs = 1

        TotCounts = float(TotCounts)
        SigPix = float(SigPix)


        #C0 = TotCounts/(2*np.pi*qobs*SigPix**2)
        C0 = TotCounts/(np.sqrt(2*np.pi)*SigPix)

        Ftot = 2*np.pi*qobs*SigPix**2*C0

        mgemag = magzpt  + 2.5*np.log10(exptime)  - 2.5*np.log10(Ftot)
        #mgemag = magzpt  + 2.5*np.log10(exptime)  - 2.5*np.log10(TotCounts)

        FWHM = 2.35482*SigPix

        h  = np.sqrt(2) * SigPix
             
        Re = (K**(0.5))*h

        outline = "Mag: {:.2f}  Sig: {:.2f}  FWHM: {:.2f} Re: {:.2f}  q: {:.2f} angle: {:.2f} \n".format(mgemag, SigPix, FWHM, Re, qobs, anglegass)
        print(outline)

        outline2 = "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} \n".format(mgemag, SigPix, FWHM, Re,  qobs, anglegass)

        fout2.write(outline2)

        PrintSersic(fout1, index+1, xpeak + 1, ypeak + 1, mgemag, Re, 0.5, qobs, anglegass, Z, fit, serfit)

        index+=1

    PrintSky(fout1, index+1, sky, Z, skyfit)
    fout1.close()
    fout2.close()

    print("Done. Gaussians are stored in {}, and {} for galfit format ".format(mgeoutfile,parfile))


def fit_1d(x,y):
    """
    Usage example for mge_fit_1d().
    This example reproduces Figure 3 in Cappellari (2002)
    It takes <1s on a 2.5 GHz computer

    """
    #n = 300  # number of sampled points
    #x = np.geomspace(0.01, 300, n)  # logarithmically spaced radii
    #y = (1 + x)**-4  # The profile must be logarithmically sampled!
    plt.clf()
    m = mge_fit_1d(x, y, ngauss=16, plot=True)
    plt.pause(1)  # allow the plot to appear in certain situations

    (counts, sigma) = m.sol
    
    return counts, sigma

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
