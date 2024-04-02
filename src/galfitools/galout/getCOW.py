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





def getCOW(galfitFile: str, dis: int, angle: float, frac: float,
        num_comp: int, plotname: str, dpival: int, galfitF2: str, maxdiff: bool) -> float:
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



    EffRadfrac, totmag = GetReff().GetReSer(head, comps, frac, theta) #obtains the radius at 95% of total flux


    #computing the Sersic indexes for different radius

    R = np.arange(0.1, EffRadfrac,.1) 

    cows, totmag = GetCOW().GetCOWrad(head, comps, theta, R)

    if galfitF2:

        #it uses the same num_comp as the first file
        head2, comps2, theta2 = readGalfitF2(galfitF2, dis, angle, num_comp) 


        cows2, totmag2 = GetCOW().GetCOWrad(head2, comps2, theta, R) #same theta as galfitF1


        diff =  np.abs(cows - cows2)

        idx = np.where(diff == np.max(diff))


        xline = R[idx][0]


        ymin = cows[idx] 
        ymax = cows2[idx]




    plt.plot(R, cows, color='blue', label='model 1', dpi=dpival)
    plt.xlabel("Rad")
    plt.ylabel("curve of growth")
    plt.grid(True)
    plt.minorticks_on()

    xmin = 0.1 
    xmax = np.max(R) 
   


    plt.gca().invert_yaxis()
    
    plt.xlabel("Radius (pixels)")
    plt.title("Curve of Growth ")
    plt.ylabel("magnitude (< R) ")


    
    if galfitF2:

        plt.plot(R, cows2, color='green',label='model 2', dpi=dpival)

        plt.xlabel("Rad")
        plt.ylabel("curve of growth")

        if maxdiff:
            plt.vlines(xline, ymin, ymax, color='red',label='max diff')


    plt.hlines(totmag, xmin, xmax, color='black',label='total magnitude')
    
    plt.legend(loc='lower right')
    plt.savefig(plotname)

    return totmag, N, theta 



def readGalfitF2(galfitF2, dis, angle, num_comp):

    galfit = Galfit(galfitF2)

    head2 = galfit.ReadHead()
    galcomps2 = galfit.ReadComps()


    galcomps2 = SelectGal(galcomps2, dis, num_comp)


    #taking the last component position angle for the whole galaxy
    maskgal = (galcomps2.Active == True) 
    if angle:
        theta2 = angle
    else:
        theta2 = galcomps2.PosAng[maskgal][-1]  

    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps2 = conver2Sersic(galcomps2) 


    N = numComps(comps2,'all')
    #print('number of model components: ',N)

    if N == 0:
        print('not enough number of components to plot COW')
        print('exiting..')
        sys.exit(1)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta2)
    #print(line)


    return head2, comps2, theta2 


### GetCOW
class GetCOW:
    '''class to obtain the Curve-of-Growth at a given radius '''


    def GetCOWrad(self, head, comps, theta, R):


        cowfluxs = np.array([])

        for r in R: #check if functions works without the for 

            cowflux = self.GetCOWFluxR(head, comps, theta, r)

            cowfluxs = np.append(cowfluxs, cowflux)



        maskgal = (comps.Active == True) 

        comps.Flux = 10**((head.mgzpt - comps.Mag)/2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + head.mgzpt


        cows = -2.5*np.log10(cowfluxs) + head.mgzpt


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



        return cowR


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

        try:

            N = bisect(self.funMeMeanMe, a, b, args=(me, meanme))

        except:
            print("solution not found for the given range")
            N = 0

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

        try:
            N = bisect(self.funReRfrac, a, b, args=(Re, Rfrac, frac))

        except:
            print("solution not found for the given range")
            N = 0



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

        try:
            K = bisect(self.funMeM0, a, b, args=(me, m0))

        except:
            print("solution not found for the given range")
            K = 0



        return K 


    def funMeM0(self, K: float, me: float, m0: float) -> float:

        result = me - m0 - 2.5*K/np.log(10)

        return result 

    def solveSerK(self, a: float, b: float, k: float) -> float: 

        try:
            sersic = bisect(self.funK, a, b, args=(k))

        except:
            print("solution not found for the given range")
            sersic = 0



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



