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


#from scipy.optimize import bisect

#import matplotlib.pyplot as plt


from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps, GetRadAng

#from galfitools.galout.getRads import GetReff, GetMe

from galfitools.galin.MaskDs9 import GetAxis

from galfitools.galin.galfit import GalComps, GalHead
from galfitools.galin.galfit import numParFree


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

    if os.path.exists(headinfo.sigimage):
        headinfo.sigimageflag = True

    if os.path.exists(headinfo.psfimage):
        headinfo.psfimageflag = True


    if os.path.exists(headinfo.maskimage):
        headinfo.maskimageflag = True

    if os.path.exists(headinfo.constraints):
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
    for idx,item in enumerate(galcomps.N):
         dx = (galcomps.PosX[idx] - galcomps.PosX)
         dy = (galcomps.PosY[idx] - galcomps.PosY)
         d = np.sqrt(dx**2 + dy**2)
         maskdis = (d <= dmax) & (galax == 0)
         if maskdis.any():
             galax[maskdis] = cont
             cont=cont+1
     
 

    freepar = numParFree(galcomps) #computing the number of free parameters

   

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


class GetN:
    '''Class to compute the Sersic index from photometric parameters'''



    def GalNs(self, EffRad, EffRadfrac, F):



        sers = np.array([])

        for idx, f in  enumerate(F):


            n = self.ReRfrac(EffRad, EffRadfrac[idx], f)

            sers = np.append(sers, n)



        return sers 





    def MeMeanMe(self, me: float, meanme: float) -> float:

        a = 0.1
        b = 12


        Sersic = self.solveSer(a, b, me, meanme)


        return Sersic



    def solveSer(self, a: float, b: float, me: float, meanme: float) -> float:
        "return the sersic index. It uses Bisection"


        N = bisect(self.funMeMeanMe, a, b, args=(me, meanme))

        return N 


    def funMeMeanMe(self, n: float, me: float, meanme: float) -> float:
        

        k = gammaincinv(2*n, 0.5)

        fn = self.Fn(n)

        result = me - meanme  - 2.5*np.log10(fn)


        return result 
     
    def Fn(self, n: float) -> float:

        k = gammaincinv(2*n, 0.5)


        fn = ((n*np.exp(k))/(k**(2*n)))*gamma(2*n)


        return fn


    def ReRfrac(self, Re: float, Rfrac: float, frac: float) -> float:

        a = 0.1
        b = 12

        Sersic = self.solveSerRe(a, b, Re, Rfrac, frac)


        return Sersic



    def solveSerRe(self, a: float, b: float, Re: float, Rfrac: float, frac: float) -> float:
        "return the sersic index. It uses Bisection"


        N = bisect(self.funReRfrac, a, b, args=(Re, Rfrac, frac))

        return N 


    def funReRfrac(self, n: float, Re: float, Rfrac: float, frac: float) -> float:
        

        k = gammaincinv(2*n, 0.5)


        x = k*(Rfrac/Re)**(1/n)
        

        result = frac*gamma(2*n) - gamma(2*n)*gammainc(2*n, x) 


        return result 
     


    def MeM0(self, me: float, m0: float) -> float:
        '''Uses me and m0 to compute n. It is not very realiable 
        and takes longer than the other two methods'''
        a = 0
        b = 45


        K = self.solveKm0(a, b, me, m0)


        a = 0.2
        b = 40



        Sersic = self.solveSerK(a, b, K)


        return Sersic

    def solveKm0(self, a: float, b: float, me: float, m0: float) -> float:
        "return the sersic index. It uses Bisection"

        K = bisect(self.funMeM0, a, b, args=(me, m0))

        return K 


    def funMeM0(self, K: float, me: float, m0: float) -> float:

        result = me - m0 - 2.5*K/np.log(10)

        return result 

    def solveSerK(self, a: float, b: float, k: float) -> float: 


        sersic = bisect(self.funK, a, b, args=(k))

        return sersic

    def funK(self, n: float, k: float) -> float:
       

        result = gammaincinv(2*n, 0.5) - k



        return result 




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



