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

from scipy.optimize import bisect

import matplotlib.pyplot as plt


from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps,GetRadAng

from galfitools.galout.getRads import GetReff, GetMe



from galfitools.galin.galfit import GalComps, GalHead





def getCOW(galfitFile: str, dis: int, angle: float, 
            num_comp: int, plotname: str) -> float, int, float:
    '''plots the curve-of-growth from the galfit.XX file. Only for Sersic functions'''




    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()



    galcomps = SelectGal(galcomps, dis, num_comp)


    #taking the last component position angle for the whole galaxy
    maskgal = (galcomps.Active == True) 
    if angle:
        theta = angle
    else:
        theta = galcomps.PosAng[maskgal][-1]  



    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 


    N = numComps(comps,'all')
    #print('number of model components: ',N)

    if N == 0:
        print('not enough number of components to plot COW')
        print('exiting..')
        sys.exit(1)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta)
    #print(line)


    #lastmod


    frac = 0.95

    EffRadfrac, totmag = GetReff().GetReSer(head, comps, frac, theta) #obtains the radius at 95% of total flux



    #computing the Sersic indexes for different radius

    R = np.arange(0.1, EffRadfrac*10,.1) 


    cows, totmag = GetCOW().GetCOWrad(head, comps, theta, R)


    #ns = GetN().GalNs(EffRad, R, F) 


    if plot:

        plt.plot(R, cows)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Radius in pixels")
        plt.ylabel("Curve of Growth")
        plt.savefig(plotname)


    

    return totmag, N, theta 



### GetCOW
class GetCOW:
    '''class to obtain the Curve-of-Growth at a given radius '''


    def GetCOWrad(self, head, comps, theta, R):



        cowfluxs = np.array([])

        for r in R: #check if functions works without the for 

            cowflux = GetCOW().GetCOWFluxR(head, comps, theta, r)

            cowfluxs = np.append(cowfluxs, cowflux)



        maskgal = (comps.Active == True) 

        comps.Flux = 10**((galhead.mgzpt - comps.Mag)/2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + galhead.mgzpt


        cows = -2.5*np.log10(cowsflux) #+ galhead.mgzpt


        return cows, totmag 




    def GetCOWFluxR(self, galhead: GalHead, comps: GalComps, theta: float, rad: float) -> float:


        maskgal = (comps.Active == True) 

        comps.Flux = 10**((galhead.mgzpt - comps.Mag)/2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + galhead.mgzpt


        #rad refers to the radius where flux will be computed
        #comps.Rad refers to the effective radius of every component

        cowR = self.funCOWSer(rad, comps.Flux[maskgal], comps.Rad[maskgal], 
                comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal],  theta)


        return cowR, totmag


    def funCOWSer(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:
        

        fun = self.COWtotser(R, flux, rad, n, q, pa, theta) #- totFlux*eff

        return fun
    

    def COWtotser(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:

        ftotR = self.COWser(R, flux, rad, n, q, pa, theta) 

        return ftotR.sum()



    def COWser(self, R: float, Flux: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''sersic flux to a determined R'''
        
        k = gammaincinv(2*n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta) 

        X = k*(Rcor/Re)**(1/n) 

        Fr = Flux*gammainc(2*n, X) 
        
        return Fr


############################
############################
#





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



