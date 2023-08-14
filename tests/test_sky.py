

import pytest
import os


import subprocess as sp



from galfitools.sky.GalfitSky import galfitSky

from galfitools.sky.Sky import sky

from galfitools.sky.SkyDs9 import SkyDs9
from galfitools.sky.SkyRing import SkyRing

def test_galfitSky():

    imgname =   "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits" 
    path="tests/"
    imgname = path+imgname
    maskfile = path+maskfile
    mgzpt =  21.817
    scale = 0.68 

    X = 169 
    Y = 544 

    initsky = 376 


    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)



    assert os.path.isfile("sky.txt")

    if os.path.isfile("sky.txt"):
        os.remove("sky.txt")





    return None


def test_sky():



    imgname =   "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits" 
    filereg = "ds9.sky.reg"
    path="tests/"
    imgname = path+imgname
    maskfile = path+maskfile
    filereg = path+filereg
 
    mean, sig = sky(imgname, maskfile, filereg)

    tol = 1e-2

    result1 =  369.717 
    result2 = 4.407 



    diffsky1 = abs(mean -result1)

    assert diffsky1 < tol

    diffsky2 = abs(sig -result2)

    assert diffsky2 < tol



    return None



def test_SkyDs9():

    ImageFile = "A671.gtMakeMask.maskds9.masksky.fits" 
    RegFile = "skyDs9.reg" 

    path="tests/"
    ImageFile = path+ImageFile
    RegFile = path+RegFile
 

    mean, sig = SkyDs9(ImageFile, RegFile) 


    tol = 1e-2

    result1 = 1146.627
    result2 = 6.919



    diffsky1 = abs(mean -result1)

    assert diffsky1 < tol

    diffsky2 = abs(sig -result2)

    assert diffsky2 < tol



    return None




def test_skyring():

    image = "A1656-1-3.sbProf.fits" 
    ds9regfile = "skyring.reg" 
    mask = "mask.sbProf.fits" 

    path="tests/"
    image = path+image 
    ds9regfile= path+ds9regfile
    mask=path+mask

    width = 20 
    center = False



    ##end input
    mean, std, median, rad = SkyRing(image, mask, ds9regfile, width, center)

    tol = 1e-2


    #update numbers
    result1 = 370.41 
    result2 = 4.43



    diffsky1 = abs(mean -result1)

    assert diffsky1 < tol

    diffsky2 = abs(std -result2)

    assert diffsky2 < tol


    #assert os.path.isfile("ring.fits")

    #if os.path.isfile("ring.fits"):
    #    os.remove("ring.fits")

    #assert os.path.isfile("ringmask.fits")

    #if os.path.isfile("ringmask.fits"):
    #    os.remove("ringmask.fits")




    return None








