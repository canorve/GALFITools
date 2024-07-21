#!/usr/bin/env python3


import numpy as np
from astropy.io import fits
import sys
import subprocess as sp
import os.path

import os

from galfitools.galin.MaskDs9 import GetAxis 
from galfitools.galin.MaskDs9 import checkCompHDU

from galfitools.galin.galfit import Galfit
from galfitools.galin.galfit import galfitLastFit 

from galfitools.galin.getStar import getStar

from galfitools.mge.mge2galfit import PrintHeader, PrintSky, PrintSersic
from galfitools.mge.mge2galfit import mge2gal


def makePSF(galfitFile: str, image: str, regfile: str, center: bool, psfout: str, sigma: str, twist: bool, numgauss: int)-> None: 


    #######################
    # reading options from galfit header

    galfit = Galfit(galfitFile)
    head = galfit.ReadHead()
    galsky = galfit.ReadSky()

    inputimage = head.inputimage

    imageout = head.outimage
    magzpt = head.mgzpt
    maskfile = head.maskimage

    sky = galsky.sky

    scale = head.scale
    psfile = head.psfimage

    sigfile = head.sigimage

    convbox =  head.convx
    convboxy = head.convy

    xlo = head.xmin
    ylo = head.ymin
    xhi = head.xmax
    yhi = head.ymax
 
    imsize = head.xmax - head.xmin +  1

    #################

    #inputimage of galfit header will be the output of getStar
    getStar(image, regfile, imsize, center, sky,  inputimage, sigma, sigfile)



    #checking the center position
    even = False

    if imsize % 2 == 0:
        even = True

    
    if even:
        lx = imsize/2

        xpeak = int(lx + 1) 
        ypeak = int(lx + 1) 
        
    else:
        lx = imsize/2

        xpeak = int(lx + 0.5) 
        ypeak = int(lx + 0.5) 


    ####

    psf = 2
    gauss = None
    freeser = False 
    freesky = False 

    xypos = [xpeak, ypeak]  
    ellip = 0.1 
    posang = 1 
    regmgefile = None



    mge2gal(galfitFile, regmgefile, center, psf, twist, gauss, freeser, freesky, numgauss, xypos=xypos, ellip = ellip, posang = posang) 



    print("calling to GALFIT to fit MGE")

    rungal = "galfit  {}".format("mseGALFIT.txt")
    errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)


    try: 
        lastfit_file = galfitLastFit(".")
    except:
        print("probably GALFIT has been unable to find a solution")
        print("Try to reduce the number of gaussians with -ng option")
        print("Exiting now")
        sys.exit(1)

    galfit = Galfit(lastfit_file)
    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()


    #printing output file to create psf model
    fout1 = open("psfmodel.txt", "w")

    head.outimage = psfout
    head.P = 1 

    PrintHeader(fout1, head.inputimage, head.outimage, head.sigimage, head.psfimage, head.psfsamp, head.maskimage, head.constraints, head.xmin, head.xmax, head.ymin,
                head.ymax, head.convx, head.convy, head.mgzpt, head.scale, head.scaley, 
                head.display, head.P, 0)


    index = 0


    for index, item in enumerate(galcomps.N):

            PrintSersic(fout1, index+1, galcomps.PosX[index], galcomps.PosY[index], 
                    galcomps.Mag[index], galcomps.Rad[index], galcomps.Exp[index], 
                    galcomps.AxRat[index], galcomps.PosAng[index], galcomps.skip[index], 
                    galcomps.MagFree[index], galcomps.ExpFree[index])




    galsky.skip = 1
 
    PrintSky(fout1, index+1, galsky.sky, galsky.skip, galsky.skyfree)


    fout1.close()


    print("calling to GALFIT to create final model")

    rungal = "galfit  {}".format("psfmodel.txt")
    errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
                universal_newlines=True)

    print("Done. PSF model file: ", psfout)





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

