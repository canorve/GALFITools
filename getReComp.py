#! /usr/bin/env python

import numpy as np
import argparse
import os
from astropy.io import fits
import subprocess as sp
import scipy


from scipy.optimize import bisect



#console scripts
def main() -> None: 
    '''gets the effective radius from a set of Sersics'''

    parser = argparse.ArgumentParser(description = "getReComp: gets the effective radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    args = parser.parse_args()
    galfitFile = args.GalfitFile
    dis = args.dis

    head = ReadHead(galfitFile)


    comps = ReadComps(galfitFile)



    num_comp = 1 #select galaxy by the number of component 
    comps = SelectGal(comps,dis,num_comp)
    

    EffRad = GetReGauss(head,comps)

    line = 'The effective radius is {:.2f} pixels \n'.format(EffRad)
    
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
            galhead.inputimage=tmp[1]

        if tmp[0] == "B)":     # out image
            galhead.outimage = tmp[1]

        if tmp[0] == "C)":   # sigma 
            galhead.sigimage = tmp[1]

        if tmp[0] == "D)":  # psf file 
            galhead.psfimage = tmp[1]

        if tmp[0] == "E)":  #psf sampling 
            galhead.psfsamp = int(tmp[1])

        if tmp[0] == "F)":     # mask image
            try:
                galhead.maskimage = tmp[1]
            except IndexError:
                galhead.maskimage = "None"

        if tmp[0] == "G)":  # psf file 
            galhead.constraints = tmp[1]


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

    Activate = np.array([])            #activate component  for galaxy

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
            galcomps.Activate = np.append(galcomps.Activate, flagcomp)


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

    idx = np.where(galcomps.N ==  n_comp)

    assert idx[0].size !=0, 'component not found'

    n_idx = idx[0].item(0)


    galcomps.Activate[n_idx] = True #main component

    posx = galcomps.PosX[n_idx] 
    posy = galcomps.PosY[n_idx]      

    dx = galcomps.PosX - posx
    dy = galcomps.PosY - posy

    dist = np.sqrt((dx)**2 + (dy)**2)

    maskdist = (dist <= distmax ) 

    galcomps.Activate[maskdist] = True

    return galcomps
 
def GetReGauss(galhead: GalHead, galcomps: GalComps) -> float:


    maskgal = (galcomps.Activate == True) & (galcomps.NameComp == "gaussian")
    #computing extra parameters of 
    galcomps.sig = galcomps.Rad / 2.354      
    galcomps.Flux = 10**((galhead.mgzpt - galcomps.Mag)/2.5)
    galcomps.Io = galcomps.Flux/(2*np.pi*galcomps.AxRat*galcomps.sig**2) 

    
    totFlux = galcomps.Flux[maskgal].sum()

    a = 0.1
    b = galcomps.Rad[-1] * 10 #hope it doesn't crash



    Reff = solveGaussRe(a, b, galcomps.Io[maskgal], galcomps.sig[maskgal], 
                        galcomps.AxRat[maskgal], totFlux)

    return Reff


def Fgauss(R: float, Io:float, sig: float, q: float) -> float:
    '''gauss flux to a determined R'''
    k = 0.6931471805599455 # k sersic constant for gaussian n = 0.5
    h = np.sqrt(2)*sig
    X = (R/h)**2 
    #it is not necessary to multiply by gamma since gamma(2*0.5)=1
    Fr = 2*np.pi*q*(sig**2)*Io*scipy.special.gammainc(1,X) 

    return Fr

def Ftotgauss(R: float, Io:float, sig: float, q: float) -> float:

   #it is not necessary to multiply by gamma since gamma(2*0.5)=1
    ftotR = Fgauss(R, Io, sig, q) 

    return ftotR.sum()


#def Ftotgauss(Io:float, sig: float, q: float) -> float:
#    Ftot = 2*np.pi*q*(sig**2)*Io
#    return Ftot

def funReGauss(R: float, Io: list, sig: list, q: list, totFlux: float) -> float:
    

    fun = Ftotgauss(R, Io, sig, q) - totFlux/2

    return fun
   

def solveGaussRe(a: float, b: float, Io: list, sig: list, q: list, totFlux: float) -> float:
    "return the Re of a set of Sersic functions. It uses Bisection"


    Re = bisect(funReGauss,a,b,args=(Io,sig,q,totFlux))

    return Re




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


