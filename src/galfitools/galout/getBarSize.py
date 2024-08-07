#! /usr/bin/env python

import numpy as np
import argparse
import os
from astropy.io import fits
import subprocess as sp
import scipy
import sys
import copy



from galfitools.galin.galfit import Galfit, conver2Sersic, SelectGal, numComps, GetRadAng

from galfitools.galin.galfit import GalComps, GalHead

from galfitools.galout.getRads import getKappa2
from galfitools.galout.getRads import getBreak2


def getBarSize(galfitFile: str, dis: int, num_comp: int, plot: bool, ranx: list, out: str) -> float:
    '''gets the Kappa radius (maximum curvature) from a set of Sersics using another method'''


    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()


    #convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps) 



    comps = SelectGal(comps, dis, num_comp)

    maskgal = (comps.Active == True) 

    theta = comps.PosAng[maskgal][1] #it assumes bar is positioned as second galfit component  

    AxRat = comps.AxRat[maskgal][1]
    X = comps.PosX[maskgal][1]
    Y = comps.PosY[maskgal][1]

    N = numComps(comps,'all')

    if N == 0:
        print('not enough number of components to compute bar size')
        print('exiting..')
        sys.exit(1)



    #########################
    ### computing the slope
    #########################
    if ranx:
        (xmin,xmax)=ranx[0], ranx[1]
    else:
        Re = comps.Rad[maskgal][1] #it assumes bar is positioned as second galfit component  

        #it assumes bar size is in this range. Hopefully it founds the solution there:
        xmin = 1
        xmax = 2.5*Re 

        ranx=[xmin,xmax]


    rbreak, N, theta =  getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)

    rkappa, N2, theta2 =  getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)

    #bar size is just the average of these two radius

    rbar = (rbreak + rkappa)/2


    #now it creates the ellipse region file 
    fout = open(out, "w")

    line = "# Region file format: DS9 version 4.1 \n"
    fout.write(line)
    line = 'global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' 
    fout.write(line)
    line = 'physical\n'
    fout.write(line)


    rbarminor = rbar*AxRat

    elline = "ellipse({:.2f}, {:.2f}, {:.2f}, {:.2f} {:.2f}) \n".format(X, Y, 
                rbarminor, rbar, theta)
    fout.write(elline)

    fout.close()

    

    return rbar, N, theta  


