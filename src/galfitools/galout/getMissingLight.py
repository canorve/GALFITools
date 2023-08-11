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


from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps, GetRadAng
from galfitools.galin.galfit import GalComps, GalHead




def getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad):
    '''gets the missing light from a two models'''

    #reading arguments parsing

    galfit1 = Galfit(galfitFile1)
    head1 = galfit1.ReadHead()
    galcomps1 = galfit1.ReadComps()


    galfit2 = Galfit(galfitFile2)
    head2 = galfit2.ReadHead()
    galcomps2 = galfit2.ReadComps()



    galcomps1 = SelectGal(galcomps1, dis, num_comp)
    galcomps2 = SelectGal(galcomps2, dis, num_comp)



    #taking the last component position angle for the whole galaxy

    maskgal1 = (galcomps1.Active == True) 
    maskgal2 = (galcomps2.Active == True) 


    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps1 = conver2Sersic(galcomps1) 
    comps2 = conver2Sersic(galcomps2) 


    N1 = numComps(comps1,'all')
    N2 = numComps(comps2,'all')

    if (N1 == 0) or (N2 == 0):
        print('not enough number of components of one of the two models to proceed')
        print('exiting..')
        sys.exit(1)


    Flux1, mag1 = GetMag().GetFlux(head1, comps1, rad)
    Flux2, mag2 = GetMag().GetFlux(head2, comps2, rad)


    #if Flux1 > Flux2:

    #    Fluxmiss = Flux1 - Flux2
    #else:
    #    print("Warning: Flux1 is less than Flux2 for the selected radius. Inverting")

    #    Fluxmiss = Flux2 - Flux1

    #magmiss = -2.5*np.log10(Fluxmiss) + head1.mgzpt


    magmiss = getMiss(head1, mag1, mag2)




    return magmiss, N1, N2




class GetMag:


    def GetFlux(self, head, comps, gamRad):


        #comps.Rad = comps.Rad*head.scale
        comps.Flux = 10**((head.mgzpt - comps.Mag)/2.5)



        k = gammaincinv(2*comps.Exp, 0.5)

        denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
        denom2 = (comps.Exp)*(k**(-2*comps.Exp))
        denom3 = (gamma(2*comps.Exp))*(comps.AxRat) #consider this
        #denom3 = (gamma(2*comps.Exp))

        denom = denom1*denom2*denom3 
        
        comps.Ie = comps.Flux/denom

        maskgal = (comps.Active == True) 


        Ftotr = self.Ftotser(gamRad, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal]) 

        mag = -2.5*np.log10(Ftotr) + head.mgzpt


        return Ftotr, mag 

     
    def Ftotser(self, R: float, Ie: list, rad: list, n: list, q: list, pa: list) -> float:

        FtotR = self.Fser(R, Ie, rad, n, q, pa ) 

        return FtotR.sum()



    def Fser(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list) -> float:
        '''sersic flux to a determined R'''


        k = gammaincinv(2*n, 0.5)

        #Rcor = GetRadAng(R, q, pa) 
        x = k * (R/Re)**(1/n)

        Fr = Ie*(Re**2)*(2*np.pi*q*n)*np.exp(k)*(k**(-2*n))*gammainc(2*n, x)*gamma(2*n)




        
        return Fr





def getMiss(head, mag1, mag2):


    Flux1 = 10**((head.mgzpt - mag1)/2.5)
    Flux2 = 10**((head.mgzpt - mag2)/2.5)

    if Flux1 > Flux2:

        Fluxmiss = Flux1 - Flux2
    else:
        print("Warning: Flux1 is less than Flux2 for the selected radius. Inverting")

        Fluxmiss = Flux2 - Flux1



    magMiss = -2.5*np.log10(Fluxmiss) + head.mgzpt

    return magMiss






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



