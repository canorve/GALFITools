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



#console scripts
def main() -> None: 
    '''gets the effective radius from a set of Sersics'''

    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getSlopeDat: gets the slope and store it in a file ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-er","--effrad", type=float, 
                        help="percentage of light to compute for radius. default=.5 for effective radius ", default=.5)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-s","--slope", type=float, 
                        help="value of slope to find. default=.5 ", default=.5)




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    slope = args.slope

    assert (eff > 0) and (eff <= 1), 'effrad must be a value between 0 and 1'
   
    num_comp =  args.numcomp



    head = ReadHead(galfitFile)


    galcomps = ReadComps(galfitFile)

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
    if args.angle:
        theta = args.angle
    else:
        theta = comps.PosAng[maskgal][-1]  




    N = numComps(comps,'all')
    print('number of model components: ',N)

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

    gam = GetSlope().GalSlope(R, comps, theta) 


    plt.plot(R, gam)
    plt.savefig("slope.png")


    #saving data

    stuff = [['nameservers','panel'], ['nameservers','panel']]
    with open("out.txt", "w") as o:
        for idx,item in enumerate(R):
            print("{} {}".format(R[idx], gam[idx]), file=o)



    rgam = GetSlope().FindSlope(comps, theta, slope) 

    line = 'The radius with slope {:.2f} is {:.2f} pixels \n'.format(slope,rgam)
    print(line)


    return None


class GalHead():
    '''store the header of galfit file'''

    inputimage = "galaxy.fits"
    outimage = "galaxy-out.fits"
    sigimage = "sig.fits"
    psfimage = "psf.fits"
    psfsamp = 1
    maskimage = "galaxy-mask.fits"
    constraints = "constraints"
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    convx = 0
    convy = 0
    mgzpt = 25
    scale = 1
    scaley = 1
    display = "regular"
    P = 0
    S = 0
   


def ReadHead(File: str) -> GalHead:
    '''reads header of galfit file'''
    inputf = File 

    galhead = GalHead() # class for header

    maskimage = ""
    #    skylevel=0

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    while index < len(lines):

        #================================================================================
        # IMAGE and GALFIT CONTROL PARAMETERS
        #A) tempfits/A2399-3-2.fits      # Input data image (FITS file)
        #B) A2399-215615.96-m073822.7-337-out.fits      # Output data image block
        #C) tempfits/none-3-2.fits      # Sigma image name (made from data if blank or "none")
        #D) psfs/PSF-1309-721.fits          # Input PSF image and (optional) diffusion kernel
        #E) 1                   # PSF fine sampling factor relative to data
        #F) mask-337            # Bad pixel mask (FITS image or ASCII coord list)
        #G) constraints         # File with parameter constraints (ASCII file)
        #H) 129  809  265  809  # Image region to fit (xmin xmax ymin ymax)
        #I) 60     60           # Size of the convolution box (x y)
        #J) 21.672              # Magnitude photometric zeropoint
        #K) 0.680  0.680        # Plate scale (dx dy)   [arcsec per pixel]
        #O) regular             # Display type (regular, curses, both)
        #P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        line = lines[index]
        (tmp) = line.split()
        if tmp[0] == "A)":     # input image
            try:
                galhead.inputimage = tmp[1]
            except IndexError:
                galhead.inputimage = "None"


        if tmp[0] == "B)":     # out image
            try:
                galhead.outimage = tmp[1]
            except IndexError:
                galhead.outimage = "None"


        if tmp[0] == "C)":   # sigma 
            try:
                galhead.sigimage = tmp[1]
            except IndexError:
                galhead.sigimage = "None"


        if tmp[0] == "D)":  # psf file 
            try:
                galhead.psfimage = tmp[1]
            except IndexError:
                galhead.psfimage = "None"


        if tmp[0] == "E)":  #psf sampling 
            galhead.psfsamp = int(tmp[1])

        if tmp[0] == "F)":     # mask image
            try:
                galhead.maskimage = tmp[1]
            except IndexError:
                galhead.maskimage = "None"

        if tmp[0] == "G)":  # constraints file 
            try:
                galhead.constraints = tmp[1]
            except IndexError:
                galhead.constraints = "None"



        if tmp[0] == "H)":     # region fit box
            galhead.xmin = int(tmp[1])
            galhead.xmax = int(tmp[2])
            galhead.ymin = int(tmp[3])
            galhead.ymax = int(tmp[4])

        if tmp[0] == "I)":     # convolution size 
            galhead.convx = int(tmp[1])
            galhead.convy = int(tmp[2])

        if tmp[0] == "J)":     # mgzpt
            galhead.mgzpt = float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            galhead.scale = float(tmp[1])
            galhead.scaley = float(tmp[2])

        if tmp[0] == "O)":     # display 
            galhead.display = tmp[1]

        if tmp[0] == "P)":     # optional output 
            galhead.P = int(tmp[1])

        if tmp[0] == "S)":     # optional output 
            galhead.S = int(tmp[1])



        index += 1

    return galhead

### class for Galfit components
class GalComps:
    '''stores the components of galfit file'''
    N = np.array([])
    NameComp = np.array([])        #0)
    PosX = np.array([])            #1)   
    PosY = np.array([])            #2)   
    Mag = np.array([])             #3)
    Rad = np.array([])             #4)
    Exp = np.array([])             #5)
    Exp2 = np.array([])            #6)  for moffat
    Exp3 = np.array([])            #7)  for moffat
                                   #8)  There is No 8 in any galfit model
    AxRat = np.array([])           #9)  AxisRatio
    PosAng = np.array([])          #10) position angle
    skip = np.array([])            #z)  skip model

    Active = np.array([])            #activate component  for galaxy

    Ie = np.array([])            # surface brightness at effective radius
    Flux = np.array([])         # Flux galax  

    # store the flags related to parameters
    FreePosX = np.array([])            #1)   
    FreePosY = np.array([])            #2)   
    FreeMag = np.array([])             #3)
    FreeRad = np.array([])             #4)
    FreeExp = np.array([])             #5)
    FreeExp2 = np.array([])            #6)  for moffat
    FreeExp3 = np.array([])            #7)  for moffat
                                   #8)  There is No 8 in any galfit model
    FreeAxRat = np.array([])           #9)  AxisRatio
    FreePosAng = np.array([])          #10) position angle

def numComps(galcomps: GalComps, name: str) -> int:
    '''obtains the number of components'''


    if name == 'all':
        nummask = (galcomps.Active == True) & ( (galcomps.NameComp == 'sersic') 
                    | (galcomps.NameComp == 'expdisk') | (galcomps.NameComp == 'gaussian')
                    | (galcomps.NameComp == 'devauc'))

    else:
        nummask = (galcomps.Active == True) & (galcomps.NameComp == name)




    N = galcomps.Active[nummask].size

    return N



def ReadComps(File: str) -> GalComps:
    '''reads all the components in the galfit file'''
    galcomps = GalComps()
     
    GalfitFile = open(File,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    N = 0


    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()

        #init values
        NameComp = "none"
        PosX = 0
        PosY = 0
        Mag = 99
        Rad = 0
        Exp = 0
        Exp2 = 0
        Exp3 = 0
        AxRat = 1
        PosAng = 0
        skip = 0

        flagcomp = False 

        freeX = 1
        freeY = 1
        Magfree = 1
        Radfree = 1
        Expfree = 1
        Exp2free = 1
        Exp3free = 1
        AxRatfree = 1
        PosAngfree = 1



        if (tmp[0] == "0)") and (tmp[1] != 'sky'):

            namec = tmp[1] 
            N = N + 1 
            NameComp = namec
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    PosX = float(tmp[1])
                    PosY = float(tmp[2])
                    freeX = int(tmp[3])
                    freeY = int(tmp[4])
                if tmp[0] == "3)" :    # mag 
                    Mag = float(tmp[1])
                    Magfree = int(tmp[2])
                if tmp[0] == "4)" :
                    Rad = float(tmp[1])
                    Radfree = int(tmp[2])
                if tmp[0] == "5)" : 
                    Exp = float(tmp[1])
                    Expfree = int(tmp[2])
                if tmp[0] == "6)": 
                    Exp2 = float(tmp[1])
                    Exp2free = int(tmp[2])
                if tmp[0] == "7)":  
                    Exp3 = float(tmp[1])
                    Exp3free = int(tmp[2])
                if tmp[0] == "9)":
                    AxRat = float(tmp[1])
                    AxRatfree = int(tmp[2])
                if tmp[0] == "10)":
                    PosAng = float(tmp[1])
                    PosAngfree = int(tmp[2])
                if tmp[0] == "Z)":
                    skip=int(tmp[1])

            galcomps.PosX = np.append(galcomps.PosX, PosX)
            galcomps.PosY = np.append(galcomps.PosY, PosY)
            galcomps.NameComp = np.append(galcomps.NameComp, NameComp)
            galcomps.N = np.append(galcomps.N, N)
                
            galcomps.Mag = np.append(galcomps.Mag, Mag)
            galcomps.Rad = np.append(galcomps.Rad, Rad)
            galcomps.Exp = np.append(galcomps.Exp, Exp)
            galcomps.Exp2 = np.append(galcomps.Exp2, Exp2)
            galcomps.Exp3 = np.append(galcomps.Exp3, Exp3)
            galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
            galcomps.PosAng = np.append(galcomps.PosAng, PosAng)
            galcomps.skip = np.append(galcomps.skip, skip)
            galcomps.Active = np.append(galcomps.Active, flagcomp)


            galcomps.FreePosX = np.append(galcomps.FreePosX, freeX)
            galcomps.FreePosY = np.append(galcomps.FreePosY, freeY)
            galcomps.FreeMag = np.append(galcomps.FreeMag, Magfree)
            galcomps.FreeRad = np.append(galcomps.FreeRad, Radfree)
            galcomps.FreeExp = np.append(galcomps.FreeExp, Expfree)
            galcomps.FreeExp2 = np.append(galcomps.FreeExp2, Exp2free)
            galcomps.FreeExp3 = np.append(galcomps.FreeExp3, Exp3free)
            galcomps.FreeAxRat = np.append(galcomps.FreeAxRat, AxRatfree)
            galcomps.FreePosAng = np.append(galcomps.FreePosAng, PosAngfree)



        index += 1

    GalfitFile.close()

    return galcomps 


class GalSky:
    '''stores the value of the GALFIT file'''

    sky = 0 
    dskyx = 0
    dskyy = 0
    skip = 0 

    skyfree = 1 
    dskyxfree = 0 
    dskyyfree = 0 
 



def ReadSky(File: str) -> GalSky:
    '''reads the sky value of the galfit file'''
    
    galsky = GalSky()

    GalfitFile = open(File,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    N = 0


    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()

        #init values
        NameComp = "sky"
        sky = 0
        dskyx = 0
        dskyy = 0
        skip = 0

        freesky = 1
        freedskyx = 0
        freedskyy = 0


        if (tmp[0] == "0)") and (tmp[1] == NameComp):

            namec = tmp[1] 

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    sky = float(tmp[1])
                    freesky = float(tmp[2])
                if tmp[0] == "2)" :    # mag 
                    dskyx = float(tmp[1])
                    freedskyx = int(tmp[2])
                if tmp[0] == "3)" :
                    dskyy = float(tmp[1])
                    freedskyy = int(tmp[2])
                if tmp[0] == "Z)":
                    skip=int(tmp[1])

            galsky.NameComp = np.append(galcomps.NameComp, NameComp)
                
            galsky.sky = np.append(galsky.sky, sky)
            galsky.dskyx = np.append(galsky.dskyx, dskyx)
            galsky.dskyy = np.append(galsky.dskyy, dskyy)
            galskyskip = np.append(galsky.skip, skip)

            galsky.skyfree = np.append(galsky.skyfree, skyfree)
            galsky.dskyxfree = np.append(galsky.dskyxfree, dskyxfree)
            galsky.dskyyfree = np.append(galsky.dskyyfree, dskyyfree)
 

        index += 1

    GalfitFile.close()

 
    return galsky



def SelectGal(galcomps: GalComps, distmax: float, n_comp: int) -> GalComps:
    '''changes Flag to true for those components who belongs
        to the same galaxy of n_comp'''

    galcomps.Active.fill(False)

    idx = np.where(galcomps.N ==  n_comp)

    assert idx[0].size !=0, 'component not found'

    n_idx = idx[0].item(0)


    galcomps.Active[n_idx] = True #main component

    posx = galcomps.PosX[n_idx] 
    posy = galcomps.PosY[n_idx]      

    dx = galcomps.PosX - posx
    dy = galcomps.PosY - posy

    dist = np.sqrt((dx)**2 + (dy)**2)

    maskdist = (dist <= distmax ) 

    galcomps.Active[maskdist] = True

    return galcomps
 



def conver2Sersic(galcomps: GalComps) -> GalComps:
    ''' function to convert exponential, gaussian params to Sersic params'''

    comps =  copy.deepcopy(galcomps)

    maskdev = comps.NameComp == "devauc"
    maskexp = comps.NameComp == "expdisk"
    maskgas = comps.NameComp == "gaussian"


    K_GAUSS = 0.6931471805599455 #constant k for gaussian
    K_EXP = 1.6783469900166612 # constant k for expdisk
    SQ2 = np.sqrt(2) 

    #for gaussian functions
    if maskgas.any():
        comps.Exp[maskgas] = 0.5 
        comps.Rad[maskgas] = comps.Rad[maskgas]/2.354 #converting to sigma 
        comps.Rad[maskgas] = SQ2*(K_GAUSS**0.5)*comps.Rad[maskgas] #converting to Re 


    #for de vaucouleurs
    if maskdev.any():
        comps.Exp[maskdev] = 4

    #for exponential disks
    if maskexp.any():
        comps.Exp[maskexp] = 1
        comps.Rad[maskexp] = K_EXP*comps.Rad[maskexp] #converting to Re


    return comps



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




def GetRadAng(R: float, q: list, pa: list, theta: float) -> float:
    '''Given an ellipse and an angle it returns the radius in angle direction. 
    Theta are the values for the galaxy and the others for every component'''


    #changing measured angle from y-axis to x-axis
    # and changing to rads:
    newpa = (pa + 90)*np.pi/180 #angle of every component
    theta = (theta + 90)*np.pi/180 #angle of direction of R 

    #bim = q * R

    ecc = np.sqrt(1 - q**2)

    alpha = theta - newpa #this is the direction 


    bell =  R*np.sqrt(1 - (ecc*np.cos(alpha))**2)


    aell = bell/q  #rad to evalue for every component



    return aell 







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



