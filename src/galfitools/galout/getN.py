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





def getN(galfitFile: str, dis: int, frac: float, angle: float, num_comp: int, plot: bool) -> float:
    '''gets the effective radius from a set of Sersics'''



    eff = 0.5

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

    #theta = 18.2534 


    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 


    N = numComps(comps,'all')
    #print('number of model components: ',N)

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta)
    #print(line)



    EffRad, totmag = GetReff().GetReSer(head, comps, eff, theta)

    line = 'Total Magnitude of the galaxy: {:.2f} \n'.format(totmag)
    #print(line)


    line = 'The radius at {:.0f}% of light is {:.2f} pixels \n'.format(eff*100,EffRad)
    #print(line)



    EffRadfrac, totmag = GetReff().GetReSer(head, comps, frac, theta)


    line = 'The radius at {:.0f}% of light is {:.2f} pixels \n'.format(frac*100, EffRadfrac)
    #print(line)



    comps2 =  copy.deepcopy(comps)
    meanme = GetMe().MeanMe(totmag, EffRad*head.scale)
    me = GetMe().Me(head, comps2, EffRad*head.scale, theta)

    line = 'Mean Surface Brightness at effective radius: {:.2f} mag/\" \n'.format(meanme)
    #print(line)


    line = 'Surface brightness at effective radius {:.2f} mag/\" \n'.format(me)
    #print(line)

    sersic = GetN().MeMeanMe(me, meanme)

    line = 'Sersic index with the method of Mean Surface Brightness at effective radius: {:.2f}  \n'.format(sersic)
    #print(line)


    sersic2 = GetN().ReRfrac(EffRad, EffRadfrac, frac)

    line = 'Sersic index with the method of fraction of light at {:.2f}: {:.2f}  \n'.format(frac,sersic2)
    #print(line)



    #computing the Sersic indexes for different radius

    Fa = np.arange(0.1, .5, .05)
    Fb = np.arange(.55,1,.05)
    F = np.concatenate((Fa,Fb))

    R = GetReff().GetRfracSer(head, comps, F, theta)


    ns = GetN().GalNs(EffRad, R, F) 


    if plot:

        plt.plot(F, ns)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Fraction of light")
        plt.ylabel("Sersic index")
        plt.savefig("Serind.png")


    #print("Sersic index computed at different fraction of light radius: ")


    #for idx, item in enumerate(F):
    #    line = 'Fraction of light: {:.2f} ; Sersic index: {:.2f} '.format(F[idx],ns[idx])
    #    print(line)

    line = '\nSersic index mean: {:.2f}  Standard deviation: {:.2f}  '.format(np.mean(ns),np.std(ns))
    #print(line)


    
    #separate in two functions:

    #return Effrad, totmag, meanme, me, N, theta 
    return sersic, np.mean(ns), np.std(ns), totmag, N, theta 


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



