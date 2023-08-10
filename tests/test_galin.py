

import pytest
import os


import subprocess as sp



from galfitools.galin.getStar import getStar
from galfitools.galin.initgal import InitGal
from galfitools.galin.MaskDs9 import maskDs9

from galfitools.galin.MakeMask import makeMask
from galfitools.galin.MaskSky import skyRem
from galfitools.galin.xy2fits import xy2fits



def test_getStar():


    image = "A671.gtMakeMask.maskds9.masksky.fits"
    regfile = "ds9.getStar.reg"
    imsize = 70 
    center = False 

    sky = 1153 
    imout = "testgetStar.fits" 

    sigma = None 
    sigout = None 

    getStar(image, regfile, imsize, center, sky, imout, sigma, sigout)


    assert os.path.isfile(imout)

    if os.path.isfile(imout):
        os.remove(imout)



    return None


def test_InitGal():

    GalfitFile = "galfit.initgal"

    number = 1 
    param3 = [1,50] 
    param4 = None 
    param5 = None 
    param6 = None 
    param7 = None 
    param8 = None 
    param9 = None 
    param10 = None 


    numcomp = 1


    InitGal(GalfitFile, number, param3, param4, param5, param6, param7, param8, param9, param10, numcomp)


    fileout="galfit-1.gal"


    assert os.path.isfile(fileout)

    if os.path.isfile(fileout):
        os.remove(fileout)


    return None



def test_maskDs9():


    MaskFile = "tempmask.fits" 
    RegFile = "maskds9.reg" 
    fill = 100 
    image = "A671.gtMakeMask.maskds9.masksky.fits"

    bor_flag = False 
    borValue = 100 


    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue) 

    assert os.path.isfile(tempmask)

    if os.path.isfile(tempmask):
        os.remove(tempmask)


    return None


def test_makeMask():

    sexfile = "cold.gtMakeMask"
    image =  "A671.gtMakeMask.maskds9.masksky.fits"
    maskfile = "tempmakemask.fits" 
    scale = 1 
    satfileout = None 


    makeMask(sexfile, image, maskfile, scale, satfileout)

    assert os.path.isfile(maskfile)

    if os.path.isfile(maskfile):
        os.remove(maskfile)


    return None



def test_skyRem():


    image =  "A671.gtMakeMask.maskds9.masksky.fits"
    mask = "tempmasksky.fits" 
    sky_mean = 1150
    sky_sig = 14
    nsig = 1 

    bor_flag = False 
    borValue = 100 


    skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)

    assert os.path.isfile(mask)

    if os.path.isfile(mask):
        os.remove(mask)



    return None

def test_xy2fits():


    ImageFile=  "A671.gtMakeMask.maskds9.masksky.fits"

    AsciiFile= "maskscii.txt" 
    Value = 100 

    xy2fits().MakeFits(ImageFile, AsciiFile, Value)

    maskfile = "maskscii.fits"

    assert os.path.isfile(maskfile)

    if os.path.isfile(maskfile):
        os.remove(maskfile)



    return None









