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

from scipy.optimize import bisect, fmin, newton

import matplotlib.pyplot as plt

from scipy.interpolate import UnivariateSpline

from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps, GetRadAng


from galfitools.galin.galfit import GalComps, GalHead



def getBreak(galfitFile: str, dis: int, inicomp: int, quick: bool, random: int, num_comp: int, angle: float, plot: bool, ranx: list)-> float:
    '''gets the break radius from a set of Sersics'''


   

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


    if plot:

        if ranx:
            (xmin,xmax)=ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin,xmax,.1)

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

    #print('finding global minimum')

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




def getBreak2(galfitFile: str, dis: int, angle: float, num_comp: int, plot: bool, ranx: list) -> float:
    '''gets the break radius from a set of Sersics using another method'''


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


    maskgal = (comps.Active == True) 
    if angle:
        theta = angle
    else:
        #taking the last component position angle for the whole galaxy
        theta = comps.PosAng[maskgal][-1]  


    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute Re')
        print('exiting..')
        sys.exit(1)



    #########################
    ### computing the slope
    #########################
    if ranx:
        (xmin,xmax)=ranx[0], ranx[1]
    else:
        xmin = 0.1
        xmax = 100


    R = np.arange(xmin,xmax,.1)

    lrad= np.log10(R)
    slp = GetSlope().GalSlope(R, comps, theta) 


    if plot:
        plt.close()
        plt.plot(R, slp)
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("Nukslope.png")


    ###
    yspl = UnivariateSpline(lrad,slp,s=0,k=4)

    yspl2d = yspl.derivative(n=2)
    yspl1d = yspl.derivative(n=1)

    if plot:
        plt.close()
        plt.plot(lrad,yspl1d(lrad))
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("Nuk2d.png")




    idx = np.where(yspl1d(lrad) == max(yspl1d(lrad)))[0][0]

    brad = lrad[idx]




    rbreak = 10**brad



    return rbreak, N, theta  
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



def getKappa(galfitFile: str, dis: int, inicomp: int, quick: bool, random: int, angle: float, num_comp: int, plot: bool, ranx: list) -> float:
    '''gets the Kappa radius from a set of Sersics'''


   

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


    if plot:

        if ranx:
            (xmin,xmax)=ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin,xmax,.1)


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


        rkappa = MultiFindKappa(comps, theta, radius) 



    return rkappa, N, theta 



def MultiFindKappa(comps, theta, radius):


    maskgal = (comps.Active == True) 

    radskappa = GetKappa().MulFindKappa(comps, theta, radius) 

    kappas = np.array([])

    print('finding global minimum:')

    for r in radskappa:

        kappa = GetKappa().funGalKappaSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

        
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

            beta = self.BetaSer(r, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta)

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

        kapparad = scipy.optimize.fmin(self.funGalKappaSer, init, args=(comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta))


        return kapparad[0] 
  

    def funGalKappaSer(self, R, Ie, Re, n, q, pa, theta):



        beta = self.BetaSer(R, Ie, Re, n, q, pa, theta)

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


    def BetaSer(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''Kappa from sersic function to a determined R'''
            
        Rcor = GetRadAng(R, q, pa, theta) 

        X = self.var_X(Rcor, Re, n)
        Xprim = self.var_Xprim(Rcor, Re, n)
        Xprim2 = self.var_Xprim2(Rcor, Re, n)

        S = self.var_S(Rcor, Ie,  Re, n, X)
        Sprim = self.var_Sprim(Rcor, Ie, Re, n, X, Xprim)
        Sprim2 = self.var_Sprim2(Rcor, Ie, Re, n, X, Xprim, Xprim2)

        Beta = (Sprim/S)*(R/np.log10(np.e)) + (
                    (Sprim2*S - Sprim**2)/S**2)*(R**2/np.log10(np.e)) 

        return Beta 

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

    return EffRad, totmag, meanme, me, N, theta 




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


    def GetRfracSer(self, head, comps, F, theta):



        rads = np.array([])

        for f in F:


            r, totmag = GetReff().GetReSer(head, comps, f, theta)

            rads = np.append(rads, r)



        return rads 





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



def getSlope(galfitFile: str, dis: int, slope: float, angle: float, num_comp: int, plot: bool, ranx: list)-> float:
    '''gets the slope radius from a set of Sersics'''


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

    if plot:

        if ranx:
            (xmin,xmax)=ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin,xmax,.1)


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



############################
############################
############################



def getBulgeRad(galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx):
    '''gets the bulge radius or the radius where two models of surface brightness models are equal'''



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

    if angle:
        theta = angle
    else:
        theta = galcomps2.PosAng[maskgal2][-1]  

    #theta = 18.2534 


    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps1 = conver2Sersic(galcomps1) 
    comps2 = conver2Sersic(galcomps2) 


    N1 = numComps(comps1,'all')
    N2 = numComps(comps2,'all')


    if (N1 == 0) or (N2 == 0):
        print('not enough number of components of one of the two models to proceed')
        print('exiting..')
        sys.exit(1)




    if plot:

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100



        R = np.arange(xmin,xmax,.1)


        Idiff = getDiff(head1, comps1, comps2, R, theta)

        plt.plot(R, Idiff)
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("BulgeRad.png")





    #computing bulge radius

    rbulge = newton(getDiffx, 0, args=(head1, comps1, comps2, theta)) 


    return rbulge, N1, N2, theta 



def getDiff(head1, comps1, comps2, R, theta):


    Idiff = np.array([])

    for r in R:

        Ir1 = GetIr().Ir(head1, comps1, r, theta)
        Ir2 = GetIr().Ir(head1, comps2, r, theta)

        Ird = Ir1 - Ir2
        Idiff = np.append(Idiff, Ird)



    return Idiff  


def getDiffx(r, head1, comps1, comps2, theta):


    Ir1 = GetIr().Ir(head1, comps1, r, theta)
    Ir2 = GetIr().Ir(head1, comps2, r, theta)

    Irdx = Ir1 - Ir2


    return Irdx





class GetIr:


    def Ir(self, head, comps, R, theta):


        #comps.Rad = comps.Rad*head.scale
        comps.Flux = 10**((head.mgzpt - comps.Mag)/2.5)

        k = gammaincinv(2*comps.Exp, 0.5)

        denom1 = (2*np.pi*comps.Rad**2)*(np.exp(k))
        denom2 = (comps.Exp)*(k**(-2*comps.Exp))
        #denom3 = (gamma(2*comps.Exp))*(comps.AxRat) 
        denom3 = (gamma(2*comps.Exp))

        denom = denom1*denom2*denom3 
        
        comps.Ie = comps.Flux/denom

        maskgal = (comps.Active == True) 


        Itotr = self.Itotser(R, comps.Ie[maskgal], comps.Rad[maskgal], comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], theta) 

        #me = -2.5*np.log10(Itotr)

        return Itotr 

     
    def Itotser(self, R: float, Ie: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:

        ItotR = self.Iser(R, Ie, rad, n, q, pa, theta) 

        return ItotR.sum()



    def Iser(self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''sersic flux to a determined R'''


        k = gammaincinv(2*n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta) 

        Ir = Ie*np.exp(-k*((Rcor/Re)**(1/n) - 1))

        
        return Ir











