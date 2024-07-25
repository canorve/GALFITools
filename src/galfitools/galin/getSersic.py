#!/usr/bin/env python3


import numpy as np
from astropy.io import fits
import sys
import subprocess as sp
import os.path


from galfitools.galout.getPeak import getPeak
from galfitools.galout.PhotDs9 import photDs9 

from galfitools.mge.mge2galfit import GetInfoEllip
from galfitools.galin.galfit import GalComps, GalHead

def getSersic(image: float, regfile: str, center: bool, maskfile: str, 
                zeropoint: float, sky: float, noprint: float, 
                bulgetot: float, bards9: str) -> GalComps:


    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)

    mag, exptime = photDs9(image, regfile, maskfile, zeropoint, sky)


    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(regfile)


    if bards9:

        Xbar, Ybar, AxRatbar, PAbar = getPeak(image, bards9, center, maskfile)
        magbar, exptimebar = photDs9(image, bards9, maskfile ,zeropoint, sky)
        objbar, xposbar, yposbar, rxbar, rybar, anglebar = GetInfoEllip(bards9)





    if bulgetot:
        Fluxtot = 10**(-mag/2.5)
        FluxBulge = Fluxtot*bulgetot
        FluxDisk = Fluxtot - FluxBulge
        mag = -2.5*np.log10(FluxBulge)
        mag2 = -2.5*np.log10(FluxDisk)

        if bards9:
            Fluxbar = 10**(-magbar/2.5) #assuming bar flux containts bulge flux
            FluxDisk = Fluxtot - Fluxbar
            FluxBulge = Fluxbar*.3   #wild guess
            Fluxbar = Fluxbar*.7   #wild guess
            mag = -2.5*np.log10(FluxBulge)
            mag2 = -2.5*np.log10(FluxDisk)
            magbar = -2.5*np.log10(Fluxbar)

            if rxbar >= rybar:
                Rebar = rxbar/2 #wild guess
            else:
                Rebar = rybar/2 #same



    if rx >= ry:
        Re = rx/2 #wild guess
    else:
        Re = ry/2 #same

    #wild guesses for n and Re
    n = 2
    skip = 0


    N = 1
    #store in GalComps data class
    galcomps = GalComps()


    if (bulgetot): 
        if not(noprint):
            print("The initial parameters for the Sersic component based on the DS9 ellipse region are: ")
            print("")
            print("WARNING: these are initial parameters. True values will be computed by GALFIT")
            print("")


            print("# Component number: 1")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
            print("3) {:.2f}   1       # Integrated magnitude ".format(mag))
            print("4) {:.2f}   1       # R_e (effective radius) ".format(Re*bulgetot)) #wild guess
            print("5) {:.2f}   1       # Sersic index n  ".format(n))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}   1       # Axis Ratio (b/a)  ".format(1))
            print("10) {:.2f}  1       # Position angle (PA)  ".format(0))
            print("Z) {}  # Skip this model in output image?  ".format(skip))

            print("")

            if bards9:
                print("# Component number: 1b bar")
                print("0) sersic # Component type")
                print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
                print("3) {:.2f}    1      # Integrated magnitude ".format(magbar))
                print("4) {:.2f}    1      # R_e (effective radius) ".format(Rebar))
                print("5) {:.2f}    0      # Sersic index n  ".format(0.5))
                print("6) 0.0000   0      # ----  ")
                print("7) 0.0000   0      # ----  ")
                print("8) 0.0000   0      # ----  ")
                print("9) {:.2f}   1       # Axis Ratio (b/a)  ".format(AxRatbar))
                print("10) {:.2f}  1       # Position angle (PA)  ".format(PAbar))
                print("Z) {}  # Skip this model in output image?  ".format(skip))
                print("")


            print("# Component number: 2")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1   # Position x, y".format(X, Y))
            print("3) {:.2f}    1      # Integrated magnitude ".format(mag2))
            print("4) {:.2f}    1      # R_e (effective radius) ".format(Re))
            print("5) {:.2f}    0      # Sersic index n  ".format(1))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}    1      # Axis Ratio (b/a)  ".format(AxRat))
            print("10) {:.2f}   1      # Position angle (PA)  ".format(PA))
            print("Z) {}  # Skip this model in output image?  ".format(skip))



            print("")
            print("for parameter constraints file: ") 
            if bards9:

                print(" 1_2_3    x    offset ")
                print(" 1_2_3    y    offset ")
            else:

                print(" 1_2    x    offset ")
                print(" 1_2    y    offset ")
 
        #first component: bulge
        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)
                        
        galcomps.Mag = np.append(galcomps.Mag, mag)
        galcomps.Rad = np.append(galcomps.Rad, Re*bulgetot)
        galcomps.Exp = np.append(galcomps.Exp, n)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, 1)
        galcomps.PosAng = np.append(galcomps.PosAng, 0)
        galcomps.skip = np.append(galcomps.skip, skip)

        #free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 1)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
 


        if bards9:

            #alternative component: bar
            N = N + 1
            galcomps.PosX = np.append(galcomps.PosX, X)
            galcomps.PosY = np.append(galcomps.PosY, Y)
            galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
            galcomps.N = np.append(galcomps.N, N)
                            
            galcomps.Mag = np.append(galcomps.Mag, magbar)
            galcomps.Rad = np.append(galcomps.Rad, Rebar)
            galcomps.Exp = np.append(galcomps.Exp, 0.5)
            galcomps.Exp2 = np.append(galcomps.Exp2, 0)
            galcomps.Exp3 = np.append(galcomps.Exp3, 0)
            galcomps.AxRat = np.append(galcomps.AxRat, AxRatbar)
            galcomps.PosAng = np.append(galcomps.PosAng, PAbar)
            galcomps.skip = np.append(galcomps.skip, skip)

            #free parameters
            galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
            galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
            galcomps.MagFree = np.append(galcomps.MagFree, 1)
            galcomps.RadFree = np.append(galcomps.RadFree, 1)
            galcomps.ExpFree = np.append(galcomps.ExpFree, 0)
            galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
            galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
            galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
     

        #second component: disk
        N = N + 1
        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)
                        
        galcomps.Mag = np.append(galcomps.Mag, mag2)
        galcomps.Rad = np.append(galcomps.Rad, Re)
        galcomps.Exp = np.append(galcomps.Exp, 1)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
        galcomps.PosAng = np.append(galcomps.PosAng, PA)
        galcomps.skip = np.append(galcomps.skip, skip)


        #free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 0)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
 


    else:

        if not(noprint):
            print("The initial parameters for the Sersic component based on the DS9 ellipse region are: ")
            print("")
            print("WARNING: these are initial parameters. True values will be computed by GALFIT")
            print("")

            print("# Component number: 1")
            print("0) sersic # Component type")
            print("1) {:.2f} {:.2f} 1 1  # Position x, y".format(X, Y))
            print("3) {:.2f}    1      # Integrated magnitude ".format(mag))
            print("4) {:.2f}    1      # R_e (effective radius) ".format(Re))
            print("5) {:.2f}    1      # Sersic index n  ".format(n))
            print("6) 0.0000   0      # ----  ")
            print("7) 0.0000   0      # ----  ")
            print("8) 0.0000   0      # ----  ")
            print("9) {:.2f}    1      # Axis Ratio (b/a)  ".format(AxRat))
            print("10) {:.2f}   1      # Position angle (PA)  ".format(PA))
            print("Z) {}  # Skip this model in output image?  ".format(skip))



        galcomps.PosX = np.append(galcomps.PosX, X)
        galcomps.PosY = np.append(galcomps.PosY, Y)
        galcomps.NameComp = np.append(galcomps.NameComp, "sersic")
        galcomps.N = np.append(galcomps.N, N)
                        
        galcomps.Mag = np.append(galcomps.Mag, mag)
        galcomps.Rad = np.append(galcomps.Rad, Re)
        galcomps.Exp = np.append(galcomps.Exp, n)
        galcomps.Exp2 = np.append(galcomps.Exp2, 0)
        galcomps.Exp3 = np.append(galcomps.Exp3, 0)
        galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
        galcomps.PosAng = np.append(galcomps.PosAng, PA)
        galcomps.skip = np.append(galcomps.skip, skip)

        #free parameters
        galcomps.PosXFree = np.append(galcomps.PosXFree, 1)
        galcomps.PosYFree = np.append(galcomps.PosYFree, 1)
        galcomps.MagFree = np.append(galcomps.MagFree, 1)
        galcomps.RadFree = np.append(galcomps.RadFree, 1)
        galcomps.ExpFree = np.append(galcomps.ExpFree, 1)
        galcomps.Exp2Free = np.append(galcomps.Exp2Free, 0)
        galcomps.Exp3Free = np.append(galcomps.Exp3Free, 0)
        galcomps.AxRatFree = np.append(galcomps.AxRatFree, 1)
 


    return galcomps 



