

import pytest
import os


import subprocess as sp



from galfitools.sky.GalfitSky import galfitSky

from galfitools.sky.Sky import sky



def test_galfitSky():

    imgname =   "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits" 

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

    mean, sig = sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:") 

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))

    tol = 1e-3

    result1 =  369.717 
    result2 = 4.407 



    diffsky1 = abs(mean -result1)

    assert diffsky1 < tol

    diffsky2 = abs(sig -result2)

    assert diffsky2 < tol



    return None



