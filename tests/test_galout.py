

import pytest
import os


import subprocess as sp

from galfitools.galout.getRads import getBreak
from galfitools.galout.getRads import getBreak2
from galfitools.galout.getRads import getFWHM
from galfitools.galout.getRads import getKappa
from galfitools.galout.getRads import getReComp
from galfitools.galout.getRads import getSlope
from galfitools.galout.getRads import getBulgeRad
from galfitools.galout.getMissingLight import getMissLight
from galfitools.galout.getN import getN
from galfitools.galout.showcube import displayCube

from galfitools.galout.PhotDs9 import photDs9 




def test_getBreak():



    galfitFile = "galfit.galout" 

    path="tests/"
    galfitFile=path+galfitFile
 


    dis = 10 

    inicomp = 2 

    quick = False 
    random = False 

    num_comp =  1 

    angle = None 


    ranx = None 
    plot = False 


    result = 248.62

    tol = 1e-2

    rbreak, N, theta = getBreak(galfitFile, dis, inicomp, quick, random, num_comp, angle, plot, ranx)
    

    diffrbreak = abs(rbreak -result)


    assert diffrbreak < tol



    return None




def test_getBreak2():

    galfitFile = "galfit.galout" 
    path="tests/"
    galfitFile=path+galfitFile
 
    dis = 10 
    angle = None 
    num_comp = 1 
    plot = False 
    ranx = None 


    rbreak, N, theta =  getBreak2(galfitFile, dis, angle, num_comp, plot, ranx)

    result = 42.7 

    tol = 1e-2


    diffrbreak = abs(rbreak -result)


    assert diffrbreak < tol





    return None

def test_getFWHM():

    galfitFile = "galfit.galout" 
    path="tests/"
    galfitFile=path+galfitFile
 
    dis = 10 
    angle = None 

   
    num_comp = 1 

    fwhm, N, theta = getFWHM(galfitFile, dis, angle, num_comp)


    result = 8.81 

    tol = 1e-2


    difffwhm = abs(fwhm - result)


    assert difffwhm < tol



    return None


def test_getKappa():

    galfitFile = "galfit.galout" 
    dis = 10 
    path="tests/"
    galfitFile=path+galfitFile
 

    inicomp = 2 
    quick = False 
    random = False 
    num_comp =  1 
    angle = None 
    ranx = None 
    plot = False 





    rkappa, N, theta = getKappa(galfitFile, dis, inicomp, quick, random, angle, num_comp, plot, ranx) 

    result = 2.62 
    tol = 1e-2


    diffrkappa = abs(rkappa - result)


    assert diffrkappa < tol




    return None


def test_getReComp():


    galfitFile = "galfit.galout" 
    path="tests/"
    galfitFile=path+galfitFile
 
    dis = 10 
    eff = .5  
    num_comp = 1 
    angle = None 

    EffRad, totmag, meanme, me, N, theta = getReComp(galfitFile, dis, eff, angle, num_comp)

    result = 97.82 
    tol = 1e-2


    diffEffRad = abs(EffRad - result)


    assert diffEffRad < tol


    return None



def test_getSlope():


    galfitFile = "galfit.galout" 
    dis = 10
    path="tests/"
    galfitFile=path+galfitFile
 
    slope = 0.5 

    num_comp = 1 

    angle = None 

    ranx = None 
    plot = False 




    rgam, N, theta = getSlope(galfitFile, dis, slope, angle, num_comp, plot, ranx)

    result = 2.62 
    tol = 1e-2


    diffrgam = abs(rgam - result)


    assert diffrgam < tol





    return None



def test_getBulgeRad():

    galfitFile1 = "galfit.1ser"
    galfitFile2 = "galfit.3ser"
    path="tests/"
    galfitFile1=path+galfitFile1
    galfitFile2=path+galfitFile2
 
    dis = 10 

    num_comp =  1 

    angle = None 

    ranx = None 
    plot = False 


    rbulge, N1, N2, theta = getBulgeRad(galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx)

    result = 3.86
    tol = 1e-2


    diffrbulge= abs(rbulge - result)


    assert diffrbulge < tol





    return None

def test_getMissLight():

    galfitFile1 = "galfit.1ser"
    galfitFile2 = "galfit.3ser"
    path="tests/"
    galfitFile1=path+galfitFile1
    galfitFile2=path+galfitFile2
 

    dis = 10 

    num_comp = 1 

    rad = 3.86  



    magmiss, N1, N2 = getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad)

    result = 18.66 
    tol = 1e-2


    diffmiss = abs(magmiss - result)


    assert diffmiss < tol




def test_getN():

    galfitFile = "galfit.galout" 
    dis = 10 
    path="tests/"
    galfitFile=path+galfitFile
 
   
    num_comp =  1

    frac = .2 #irrelevant variable, it is not used anymore but it was kept 

    angle = None

    plot = False 


    sersic, meanser, stdser, totmag, N, theta = getN(galfitFile, dis, frac, angle, num_comp, plot)




    result1 = 3.48 
    result2 = 2.66 

    tol = 1e-2


    diffser = abs(sersic - result1)

    assert diffser < tol



    diffser2 = abs(meanser - result2)

    assert diffser2 < tol



    return None


def test_displayCube():

    cubeimage =  "A1656-showcube.fits"
    namecube = "cube.png" 
    dpival = 100 
    brightness = 0
    contrast = 1 
    cmap = "viridis" 
    scale = 1 
    noplot = False 

    path="tests/"
    cubeimage=path+cubeimage
    namecube=path+namecube
 
    displayCube(cubeimage, namecube, dpival, brightness, contrast, cmap, scale, noplot)

    assert os.path.isfile(namecube)

    if os.path.isfile(namecube):
        os.remove(namecube)




    return None


def test_photDs9():

    ImageFile =   "A671.gtMakeMask.maskds9.masksky.fits"
    RegFile = "maskds9.reg" 
    path="tests/"
    ImageFile=path+ImageFile
    RegFile=path+RegFile
 


    zeropoint = 25 
    sky = 0 


    mag = photDs9(ImageFile, RegFile, zeropoint, sky)



    result =  11.71
    tol = 1e-2


    diffmag = abs(mag - result)


    assert diffmag < tol




    return None





