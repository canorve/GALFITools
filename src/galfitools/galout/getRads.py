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

from scipy.optimize import bisect, fmin

import matplotlib.pyplot as plt


from galfitools.galin.galfit import Galfit, conver2sersic, SelectGal, numComps




def getBreak(galfitFile: str, dis: int, eff: float, inicomp: int, quick: bool, random: int, num_comp: int, angle: float)-> float:
    '''gets the break radius from a set of Sersics'''


    assert (eff > 0) and (eff <= 1), 'effrad must be a value between 0 and 1'
   

    #head = ReadHead(galfitFile)
    #galcomps = ReadComps(galfitFile)


    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()



    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 



    comps.Flux = 10**((-comps.Mag)/2.5)

    k = gammaincinv(2*comps.Exp, 0.5)

    denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
    denom2 = (comps.Exp)*(k**(-2*comps.Exp))
    denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 

    denom = denom1*denom2*denom3 
    
    comps.Ie = comps.Flux/denom



    comps = SelectGal(comps, dis, num_comp)

    #taking the last component position angle for the whole galaxy

    maskgal = (comps.Active == True) 

    if angle:
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]  




    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)



    #########################
    ### computing the slope
    #########################


    R = np.arange(0.1,100,.1)

    gam = GetBreak().GalBreak(R, comps, theta) 


    plt.plot(R, gam)
    plt.grid(True)
    plt.minorticks_on()
    plt.savefig("Break.png")


    if quick:

        rbreak = GetBreak().FindBreak(comps, theta, inicomp) 

        line = 'The break radius  is {:.2f} pixels \n'.format(rbreak)
        #print(line)

    else:

        if random:

            radius = np.random.random(random)*comps.Rad[maskgal][inicomp] 
            #print('The initial search radius are: ',radius)

        else:
    
            radius = comps.Rad[maskgal]
            #print('The initial search radius are the effective radius of the components')


        rbreak = MulFindBreak(comps, theta, radius) 



    return rbreak, N, theta


def MulFindBreak(comps, theta, radius):

    maskgal = (comps.Active == True) 

    radsbreak = GetBreak().MulFindBreak(comps, theta, radius) 

    betas = np.array([])

    #print('finding global minium')

    for r in radsbreak:

        beta = GetBreak().BreakSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

        
        betas = np.append(betas, beta)


    radmask = (radsbreak < radius[-1]) #hope it works

    idx = np.where(betas == min(betas))  

    
    return radsbreak[idx][0]


class GetBreak:
    '''class to obtain the break radius for the whole galaxy'''



    def FullSlopeSer(self, R: float, Re: list, n: list, q: list, pa: list, theta: float) -> float:


        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta) 

        return SlptotR.sum()


    def GalBreak(self, R: list, comps: GalComps, theta: float) -> float:
        "plots the curvature of the function"

        maskgal = (comps.Active == True) #using active components only 

        kappa = np.array([])

        for r in R:

            beta = self.BreakSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

            kappa = np.append(kappa, beta)


        return kappa 


    def funGalSlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float, slope: float) -> float:


        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta)  - slope


        return fun
     

    def FindBreak(self, comps: GalComps, theta: float, initial_comp: int) -> float:
        "return the break radii of a set of Sersic functions"

        maskgal = (comps.Active == True) #using active components only 

        init = comps.Rad[maskgal][initial_comp] #component as initial value, hope it works 

        breakrad = scipy.optimize.fmin(self.funGalBreakSer, init, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta))


        return breakrad[0] 


    def MulFindBreak(self, comps: GalComps, theta: float, radius: list) -> float:
        "return the break radius evaluated at different effective radius"


        brads = np.array([])

        maskgal = (comps.Active == True) #using active components only 


        for idx, item in enumerate(radius):

            init = item 

            breakrad = scipy.optimize.fmin(self.funGalBreakSer, init, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta))


            line = 'Optimized radius: {:.2f} \n '.format(breakrad[0])
            print(line)


            brads = np.append(brads, breakrad[0])


        return brads 
 



    def funGalBreakSer(self, R, Ie, Re, n, q, pa, theta):


        beta = self.BreakSer(R, Ie, Re, n, q, pa, theta)

           

        return beta 


    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"

        maskgal = (comps.Active == True) #using active components only 

        a = 0.1 
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        Radslp = bisect(self.funGalSlopeSer, a, b, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta, slope))

        return Radslp 


    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varx = k*((R/Re)**(1/n) - 1)


        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr = (k/n)*((R/Re)**(1/n - 1))*(1/Re) 

        return varxpr

    def var_Xprim2(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr2 = (k/n)*((R/Re)**(1/n - 2))*(1/Re**2)*(1/n -1) 

        return varxpr2


    def var_S(self, R: float, Ie: list,  Re: list, n: list, X: list):

        S = Ie*np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list,  Re: list, n: list, X: list, Xprim: list):

        Sprim = -Ie*np.exp(-X)*Xprim

        return Sprim.sum()
    
    def var_Sprim2(self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list, Xprim2: list):

        Sprim2 = Ie*np.exp(-X)*Xprim**2 - Ie*np.exp(-X)*Xprim2 

        return Sprim2.sum()


    def BreakSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''Break from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)
        Xprim2 = self.var_Xprim2(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)
        Sprim2 = self.var_Sprim2(Rcor, Ie, Re, n, X, Xprim, Xprim2)

        Break = (Sprim/S)*(R/np.log10(np.e)) + (
                    (Sprim2*S - Sprim**2)/S**2)*(R**2/np.log10(np.e)) 

        return Break 

    def SlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''slope from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie,  Re, n, X, Xprim)

        Slp = (Sprim/S)*R 


        return Slp


############################
############################
############################



def getFWHM(galfitFile: str, dis: int, angle: float, num_comp: int):
    '''gets the FWHM from a set of Sersics'''


    #head = ReadHead(galfitFile)
    #galcomps = ReadComps(galfitFile)


    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
 


    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 



    comps.Flux = 10**((-comps.Mag)/2.5)

    k = gammaincinv(2*comps.Exp, 0.5)

    denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
    denom2 = (comps.Exp)*(k**(-2*comps.Exp))
    denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 

    denom = denom1*denom2*denom3 
    
    comps.Ie = comps.Flux/denom



    comps = SelectGal(comps, dis, num_comp)

    #taking the last component position angle for the whole galaxy

    maskgal = (comps.Active == True) 

    if angle:
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]  




    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute FWHM')
        print('exiting..')
        sys.exit(1)



    #########################
    ### computing the slope
    #########################



    fwhm = GetFWHM().FindFWHM(comps, theta) 


    return fwhm, N, theta 



class GetFWHM:
    '''class to obtain the FWHM for the whole galaxy'''



    def FullFWHMSer(self, R: float, Re: list, n: list, q: list, pa: list, theta: float) -> float:


        SlptotR = self.FWHMSer(R, Re, n, q, pa, theta) 

        return SlptotR.sum()


    def GalFWHM(self, R: list, comps: GalComps, theta: float) -> float:

        maskgal = (comps.Active == True) #using active components only 

        gam = np.array([])

        for r in R:

            slp = self.FWHM(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)
            gam = np.append(gam, slp)



        return gam


    def funGalFWHMSer(self, R: float, S0: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:


        fun = self.FWHMSer(R, Ie, Re, n, q, pa, theta)  -  S0/2  


        return fun
     



    def FindFWHM(self, comps: GalComps, theta: float) -> float:
        "return the fwhm of a set of Sersic functions. It uses Bisection"

        maskgal = (comps.Active == True) #using active components only 

        #find the surface brightness at S(R=0)
        X = self.var_X(0, comps.Rad[maskgal], comps.Exp[maskgal])
        S0 = self.var_S(0, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], X)
 



        a = 0.001 
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        Radslp = bisect(self.funGalFWHMSer, a, b, args=(S0,comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta ))

        return 2*Radslp 


    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varx = k*((R/Re)**(1/n) - 1)

        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr = (k/n)*((R/Re)**(1/n - 1))*(1/Re) 

        return varxpr

    def var_S(self, R: float, Ie: list,  Re: list, n: list, X: list):

        S = Ie*np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list,  Re: list, n: list, X: list, Xprim: list):

        Sprim = Ie*np.exp(-X)*Xprim

        return Sprim.sum()

    def FWHMSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''I(R) from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        #Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        #Sprim = self.var_Sprim(Rcor, Ie,  Re, n, X, Xprim)

        Ifwhm = S 


        return Ifwhm 

############################
############################
############################



def getKappa(galfitFile: str, dis: int, eff: float, inicomp: int, quick: bool, random: int, angle: float, num_comp: int) -> float:
    '''gets the Kappa radius from a set of Sersics'''


    assert (eff > 0) and (eff <= 1), 'effrad must be a value between 0 and 1'
   

    #head = ReadHead(galfitFile)
    #galcomps = ReadComps(galfitFile)


    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
 




    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 



    comps.Flux = 10**((-comps.Mag)/2.5)

    k = gammaincinv(2*comps.Exp, 0.5)

    denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
    denom2 = (comps.Exp)*(k**(-2*comps.Exp))
    denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 

    denom = denom1*denom2*denom3 
    
    comps.Ie = comps.Flux/denom



    comps = SelectGal(comps, dis, num_comp)

    #taking the last component position angle for the whole galaxy

    maskgal = (comps.Active == True) 
    if angle:
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]  




    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)

    line = 'Using a theta value of : {:.2f} degrees\n'.format(theta)
    print(line)



    #########################
    ### computing the slope
    #########################


    R = np.arange(0.1,100,.1)

    gam = GetKappa().GalKappa(R, comps, theta) 


    plt.plot(R, gam)
    plt.grid(True)
    plt.minorticks_on()
    plt.savefig("Kappa.png")


    if quick:

        rkappa = GetKappa().FindKappa(comps, theta, inicomp) 

        line = 'The Kappa radius  is {:.2f} pixels \n'.format(rbreak)
        #print(line)

    else:

        if random:

            radius = np.random.random(random)*comps.Rad[maskgal][inicomp] 
            #print('The initial search radius are: ',radius)

        else:
    
            radius = comps.Rad[maskgal]
            #print('The initial search radius are the effective radius of the components')


        rkappa = MulFindKappa(comps, theta, radius) 



    return rkappa, N, theta 


def MulFindKappa(comps, theta, radius):


    maskgal = (comps.Active == True) 

    radskappa = GetKappa().MulFindKappa(comps, theta, radius) 

    kappas = np.array([])

    #print('finding global minium:')

    for r in radskappa:

        kappa = GetKappa().KappaSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

        
        kappas = np.append(kappas, kappa)


    radmask = (radskappa < radius[-1]) #hope it works

    idx = np.where(kappas == min(kappas))  

    
    return radskappa[idx][0]



class GetKappa:
    '''class to obtain the Kappa radius for the whole galaxy'''



    def FullSlopeSer(self, R: float, Re: list, n: list, q: list, pa: list, theta: float) -> float:


        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta) 

        return SlptotR.sum()


    def GalKappa(self, R: list, comps: GalComps, theta: float) -> float:
        "plots the curvature of the function"

        maskgal = (comps.Active == True) #using active components only 

        kappa = np.array([])

        for r in R:

            beta = self.KappaSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

            gam = self.SlopeSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

            kap = np.abs(beta)/(1 + gam**2)**(3/2)

            kappa = np.append(kappa, kap)


        return kappa 


    def funGalSlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float, slope: float) -> float:


        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta)  - slope


        return fun
     

    def FindKappa(self, comps: GalComps, theta: float, initial_comp: int) -> float:
        "return the break radius of a set of Sersic functions"

        maskgal = (comps.Active == True) #using active components only 

        init = comps.Rad[maskgal][initial_comp] #component as initial value, hope it works 

        breakrad = scipy.optimize.fmin(self.funGalKappaSer, init, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta))


        return breakrad[0] 
  

    def funGalKappaSer(self, R, Ie, Re, n, q, pa, theta):



        beta = self.KappaSer(R, Ie, Re, n, q, pa, theta)

        gam = self.SlopeSer(R, Ie, Re, n, q, pa, theta)
           
        kappa = np.abs(beta)/(1 + gam**2)**(3/2)


        return -kappa 


    def MulFindKappa(self, comps: GalComps, theta: float, radius: list) -> float:
        "return the kappa radius evaluated at different effective radius"


        krads = np.array([])

        maskgal = (comps.Active == True) #using only active components 


        for idx, item in enumerate(radius):

            init = item 

            kapparad = scipy.optimize.fmin(self.funGalKappaSer, init, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta))

            
            line = 'Optimized radius: {:.2f} \n '.format(kapparad[0])
            print(line)

            krads = np.append(krads, kapparad[0])


        return krads 
 




    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"

        maskgal = (comps.Active == True) #using active components only 

        a = 0.1 
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        Radslp = bisect(self.funGalSlopeSer, a, b, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta, slope))

        return Radslp 


    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varx = k*((R/Re)**(1/n) - 1)


        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr = (k/n)*((R/Re)**(1/n - 1))*(1/Re) 

        return varxpr

    def var_Xprim2(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr2 = (k/n)*((R/Re)**(1/n - 2))*(1/Re**2)*(1/n -1) 

        return varxpr2


    def var_S(self, R: float, Ie: list,  Re: list, n: list, X: list):

        S = Ie*np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list,  Re: list, n: list, X: list, Xprim: list):

        Sprim = -Ie*np.exp(-X)*Xprim

        return Sprim.sum()
    
    def var_Sprim2(self, R: float, Ie: list, Re: list, n: list, X: list, Xprim: list, Xprim2: list):

        Sprim2 = Ie*np.exp(-X)*Xprim**2 - Ie*np.exp(-X)*Xprim2 

        return Sprim2.sum()


    def KappaSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''Kappa from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)
        Xprim2 = self.var_Xprim2(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)
        Sprim2 = self.var_Sprim2(Rcor, Ie, Re, n, X, Xprim, Xprim2)

        Kappa = (Sprim/S)*(R/np.log10(np.e)) + (
                    (Sprim2*S - Sprim**2)/S**2)*(R**2/np.log10(np.e)) 

        return Kappa 

    def SlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''slope from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie,  Re, n, X, Xprim)

        Slp = (Sprim/S)*R 


        return Slp



############################
############################
############################



def getReComp(galfitFile: str, dis: int, eff: float, angle: float, num_comp: int)-> float:
    '''gets the effective radius from a set of Sersics'''


    assert (eff > 0) and (eff <= 1), 'effrad must be a value between 0 and 1'


    #head = ReadHead(galfitFile)
    #galcomps = ReadComps(galfitFile)

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

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)



    EffRad, totmag = GetReff().GetReSer(head, comps, eff, theta)

    meanme = GetMe().MeanMe(totmag, EffRad*head.scale)
    me = GetMe().Me(head, comps, EffRad*head.scale, theta)

    return Effrad, totmag, meanme, me, N, theta 




class GetMe:
    '''Class to obtain the surface brightness at effective radius'''

    def MeanMe(self, magtot: float, effrad: float) -> float:

        meanme = magtot + 2.5*np.log10(2*np.pi*effrad**2)

        return meanme


    def Me(self, head, comps, EffRad, theta):

        comps.Rad = comps.Rad*head.scale
        comps.Flux = 10**((-comps.Mag)/2.5)

        k = gammaincinv(2*comps.Exp, 0.5)

        denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
        denom2 = (comps.Exp)*(k**(-2*comps.Exp))
        #denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 
        denom3 = (gamma(2*comps.Exp))

        denom = denom1*denom2*denom3 
        
        comps.Ie = comps.Flux/denom


        maskgal = (comps.Active == True) 
        Itotr = self.Itotser(EffRad, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta) 

        me = -2.5*np.log10(Itotr)

        return me

     
    def Itotser(self, R: float, Ie: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:

        ItotR = self.Iser(R, Ie, rad, n, q, pa, theta) 

        return ItotR.sum()



    def Iser(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''sersic flux to a determined R'''


        k = gammaincinv(2*n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta) 

        Ir = Ie*np.exp(-k*((Rcor/Re)**(1/n) - 1))

        
        return Ir



### Sersic components 
class GetReff:
    '''class to obtain the effective radius for the whole galaxy'''

    def GetReSer(self, galhead: GalHead, comps: GalComps, eff: float, theta: float) -> float:


        maskgal = (comps.Active == True) 

        comps.Flux = 10**((galhead.mgzpt - comps.Mag)/2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + galhead.mgzpt

        a = 0.1
        b = comps.Rad[maskgal][-1] * 1000  # hope it doesn't crash

        Reff = self.solveSerRe(a, b, comps.Flux[maskgal], comps.Rad[maskgal], 
                comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], totFlux, eff, theta)


        return Reff, totmag



    def solveSerRe(self, a: float, b: float, flux: list, rad: list, n: list, q: list, pa: list, totFlux: float, eff: float, theta: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"


        Re = bisect(self.funReSer, a, b, args=(flux, rad, n, q, pa, totFlux, eff, theta))

        return Re


    def funReSer(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, totFlux: float, eff: float, theta: float) -> float:
        

        fun = self.Ftotser(R, flux, rad, n, q, pa, theta) - totFlux*eff

        return fun
     
    def Ftotser(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:

        ftotR = self.Fser(R, flux, rad, n, q, pa, theta) 

        return ftotR.sum()



    def Fser(self, R: float, Flux: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''sersic flux to a determined R'''
        
        k = gammaincinv(2*n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta) 

        X = k*(Rcor/Re)**(1/n) 

        Fr = Flux*gammainc(2*n, X) 
        
        return Fr


############################
############################
############################



def getSlope(galfitFile: str, dis: int, eff: float, slope: float, angle: float, num_comp: int)-> float:
    '''gets the slope radius from a set of Sersics'''


    assert (eff > 0) and (eff <= 1), 'effrad must be a value between 0 and 1'

    #head = ReadHead(galfitFile)
    #galcomps = ReadComps(galfitFile)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()
 
    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 



    comps.Flux = 10**((-comps.Mag)/2.5)

    k = gammaincinv(2*comps.Exp, 0.5)

    denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
    denom2 = (comps.Exp)*(k**(-2*comps.Exp))
    denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 

    denom = denom1*denom2*denom3 
    
    comps.Ie = comps.Flux/denom



    comps = SelectGal(comps, dis, num_comp)

    #taking the last component position angle for the whole galaxy

    maskgal = (comps.Active == True) 

    if angle:
        theta = angle
    else:
        theta = comps.PosAng[maskgal][-1]  




    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)



    #########################
    ### computing the slope
    #########################


    R = np.arange(0.1,100,.1)

    gam = GetSlope().GalSlope(R, comps, theta) 


    plt.plot(R, gam)
    plt.grid(True)
    plt.minorticks_on()
    plt.savefig("slope.png")


    rgam = GetSlope().FindSlope(comps, theta, slope) 


    return rgam, N, theta 




class GetSlope:
    '''class to obtain the effective radius for the whole galaxy'''



    def FullSlopeSer(self, R: float, Re: list, n: list, q: list, pa: list, theta: float) -> float:


        SlptotR = self.SlopeSer(R, Re, n, q, pa, theta) 

        return SlptotR.sum()


    def GalSlope(self, R: list, comps: GalComps, theta: float) -> float:

        maskgal = (comps.Active == True) #using active components only 

        gam = np.array([])

        for r in R:

            slp = self.SlopeSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)
            gam = np.append(gam, slp)



        return gam


    def funGalSlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float, slope: float) -> float:


        fun = self.SlopeSer(R, Ie, Re, n, q, pa, theta)  - slope


        return fun
     



    def FindSlope(self, comps: GalComps, theta: float, slope: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"

        maskgal = (comps.Active == True) #using active components only 

        a = 0.1 
        b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

        Radslp = bisect(self.funGalSlopeSer, a, b, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta, slope))

        return Radslp 


    def var_X(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varx = k*((R/Re)**(1/n) - 1)


        return varx

    def var_Xprim(self, R: float, Re: list, n: list):

        k = gammaincinv(2*n, 0.5)

        varxpr = (k/n)*((R/Re)**(1/n - 1))*(1/Re) 


        return varxpr

    def var_S(self, R: float, Ie: list,  Re: list, n: list, X: list):

        S = Ie*np.exp(-X)

        return S.sum()

    def var_Sprim(self, R: float, Ie: list,  Re: list, n: list, X: list, Xprim: list):

        Sprim = Ie*np.exp(-X)*Xprim

        return Sprim.sum()

    def SlopeSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''slope from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie,  Re, n, X, Xprim)

        Slp = (Sprim/S)*R 


        return Slp








