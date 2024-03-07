#! /usr/bin/env python

import numpy as np
import argparse
import os
from astropy.io import fits
import subprocess as sp
import scipy
import sys
import copy

from scipy.special import gamma, gammainc, gammaincinv


from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps, GetRadAng


from galfitools.galin.MaskDs9 import GetAxis

from galfitools.galin.galfit import GalComps, GalHead
from galfitools.galin.galfit import numParFree, numParSkyFree



class HeadInfo():
    '''store the flags of the header of the galfit file'''

    inputimageflag =  False   
    outimageflag = False 
    sigimageflag = False 
    psfimageflag = False 
    maskimageflag = False 
    constraintsflag = False 
    xsizeflag = False 
    ysizeflag = False
    convxflag = False 
    convyflag = False 


    inputimage = "none.fits"
    outimage = "none-out.fits"
    sigimage = "none"
    psfimage = "none"
    psfsamp = 1
    maskimage = "none"
    constraints = "none"
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    convx = 1
    convy = 1
    mgzpt = 25
    scale = 1
    scaley = 1
    display = "regular"
    P = 0







def checkFile(galfitFile: str, dis: int) -> float:
    '''prints information of the galfit File '''




    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    galsky = galfit.ReadSky()

    headinfo  = HeadInfo()

    #copying info from header

    headinfo = copy2info(head, headinfo)


    galcomps = SelectGal(galcomps, dis, 1) #modify functions such that 0 will clean galcomps


    #================================================================================
    # IMAGE and GALFIT CONTROL PARAMETERS
    #A) A407.fits           # Input data image (FITS file)
    #B) A407-G1-mge.fits      # Output data image block
    #C) None                # Sigma image name (made from data if blank or "none") 
    #D) psf.fits            # Input PSF image and (optional) diffusion kernel
    #E) 1                   # PSF fine sampling factor relative to data 
    #F) mask.fits           # Bad pixel mask (FITS image or ASCII coord list)
    #G) cons.G1G2           # File with parameter constraints (ASCII file) 
    #H) 613  1290 882  1275 # Image region to fit (xmin xmax ymin ymax)
    #I) 71     71           # Size of the convolution box (x y)
    #J) 21.4710             # Magnitude photometric zeropoint 
    #K) 0.6800  0.6800      # Plate scale (dx dy)   [arcsec per pixel]
    #O) regular             # Display type (regular, curses, both)
    #P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps


    if os.path.exists(headinfo.inputimage):
        headinfo.inputimageflag = True

    if os.path.exists(headinfo.outimage):
        headinfo.outimageflag = True

    if (os.path.exists(headinfo.sigimage)) or (headinfo.sigimage == 'none') :
            headinfo.sigimageflag = True

    if (os.path.exists(headinfo.psfimage)) or (headinfo.psfimage == 'none') :
            headinfo.psfimageflag = True
    
    if (os.path.exists(headinfo.maskimage)) or (headinfo.maskimage == 'none'):
            headinfo.maskimageflag = True

    if (os.path.exists(headinfo.constraints)) or (headinfo.constraints == 'none'):
            headinfo.constraintsflag = True

    if headinfo.xmax > headinfo.xmin:
        headinfo.xsizeflag = True

    if headinfo.ymax > headinfo.ymin:
        headinfo.ysizeflag = True


    if headinfo.psfimageflag:

        (nx, ny) = GetAxis(headinfo.psfimage)
        

        if headinfo.convx > nx:
            headinfo.convxflag = True

        if headinfo.convy > ny:
            headinfo.convyflag = True

    

    galax = np.zeros(len(galcomps.N))
    mag = np.zeros(len(galcomps.N))

    mag.fill(99)


    mag =np.copy(galcomps.Mag) 
    cont=1

    dmax = dis

    #counting components per galaxy
    for idx, item in enumerate(galcomps.N):
         dx = (galcomps.PosX[idx] - galcomps.PosX)
         dy = (galcomps.PosY[idx] - galcomps.PosY)
         d = np.sqrt(dx**2 + dy**2)
         maskdis = (d <= dmax) & (galax == 0)
         if maskdis.any():
             galax[maskdis] = cont
             cont=cont+1
     
 

    galcomps.Active = True
    freepar = numParFree(galcomps) #computing the number of free parameters
    freeparsky = numParSkyFree(galsky) #computing the number of free parameters

    freepar = freepar + freeparsky
   

    return headinfo, galax, mag, freepar


def copy2info(head, head2):

    head2.inputimage = head.inputimage
    head2.outimage = head.outimage 
    head2.sigimage = head.sigimage 
    head2.psfimage = head.psfimage 
    head2.psfsamp = head.psfsamp
    head2.maskimage =  head.maskimage  
    head2.constraints = head.constraints 
    head2.xmin = head.xmin
    head2.xmax =  head.xmax
    head2.ymin =  head.ymin
    head2.ymax =  head.ymax
    head2.convx =  head.convx
    head2.convy = head.convy
    head2.mgzpt = head.mgzpt
    head2.scale = head.scale 
    head2.scaley = head.scaley
    head2.display = head.display
    head2.P =  head.P





    return head2


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



