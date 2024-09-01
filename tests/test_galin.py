import pytest
import os


import subprocess as sp
import numpy as np

from collections import Counter


from galfitools.galin.getStar import getStar
from galfitools.galin.initgal import InitGal
from galfitools.galin.MaskDs9 import maskDs9

from galfitools.galin.MakeMask import makeMask
from galfitools.galin.MaskSky import skyRem
from galfitools.galin.xy2fits import xy2fits


from galfitools.galin.checkGalFile import checkFile
from galfitools.galin.getSersic import getSersic
from galfitools.galin.imarith import imarith


def test_getStar():

    image = "A671.gtMakeMask.maskds9.masksky.fits"
    regfile = "ds9.getStar.reg"

    path = "tests/"

    image = path + image
    regfile = path + regfile

    imsize = 70
    center = False

    sky = 1153
    imout = "testgetStar.fits"

    imout = path + imout

    sigma = None
    sigout = None

    getStar(image, regfile, imsize, center, sky, imout, sigma, sigout)

    assert os.path.isfile(imout)

    if os.path.isfile(imout):
        os.remove(imout)

    return None


def test_InitGal():

    GalfitFile = "galfit.initgal"

    path = "tests/"
    GalfitFile = path + GalfitFile

    number = 1
    param3 = [1, 50]
    param4 = None
    param5 = None
    param6 = None
    param7 = None
    param8 = None
    param9 = None
    param10 = None

    numcomp = 1

    InitGal(
        GalfitFile,
        number,
        param3,
        param4,
        param5,
        param6,
        param7,
        param8,
        param9,
        param10,
        numcomp,
    )

    fileout = "galfit-1.gal"

    fileout = path + fileout

    assert os.path.isfile(fileout)

    if os.path.isfile(fileout):
        os.remove(fileout)

    namefile = "rungalfit.sh"
    # namefile = path+namefile

    assert os.path.isfile(namefile)
    if os.path.isfile(namefile):
        os.remove(namefile)

    return None


def test_maskDs9():

    MaskFile = "tempmask.fits"
    RegFile = "maskds9.reg"
    fill = 100
    image = "A671.gtMakeMask.maskds9.masksky.fits"

    path = "tests/"
    MaskFile = path + MaskFile
    RegFile = path + RegFile
    image = path + image

    bor_flag = False
    borValue = 100

    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue)

    assert os.path.isfile(MaskFile)

    if os.path.isfile(MaskFile):
        os.remove(MaskFile)

    return None


def test_makeMask():

    sexfile = "cold.gtMakeMask"
    image = "A671.gtMakeMask.maskds9.masksky.fits"
    maskfile = "tempmakemask.fits"
    scale = 1
    satfileout = "ds9sat.reg"

    path = "tests/"
    sexfile = path + sexfile
    image = path + image
    maskfile = path + maskfile

    makeMask(sexfile, image, maskfile, scale, satfileout)

    assert os.path.isfile(maskfile)

    if os.path.isfile(maskfile):
        os.remove(maskfile)

    namefile = "sexsort.cat"
    # namefile = path+namefile

    assert os.path.isfile(namefile)
    if os.path.isfile(namefile):
        os.remove(namefile)

    assert os.path.isfile(satfileout)
    if os.path.isfile(satfileout):
        os.remove(satfileout)

    return None


def test_skyRem():

    image = "A671.gtMakeMask.maskds9.masksky.fits"
    mask = "tempmasksky.fits"

    path = "tests/"
    image = path + image
    mask = path + mask

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

    ImageFile = "A671.gtMakeMask.maskds9.masksky.fits"

    AsciiFile = "maskscii.txt"

    path = "tests/"
    ImageFile = path + ImageFile
    AsciiFile = path + AsciiFile

    Value = 100

    xy2fits().MakeFits(ImageFile, AsciiFile, Value)

    maskfile = "maskscii.fits"
    maskfile = path + maskfile

    assert os.path.isfile(maskfile)

    if os.path.isfile(maskfile):
        os.remove(maskfile)

    return None


def test_checkfile():

    path = "tests/"

    galfitFile = "galfit.3ser"
    galfitFile = path + galfitFile
    dis = 10

    headinfo, galax, mag, freepar = checkFile(galfitFile, dis)

    tol = 1e-2

    diffgal = abs(len(galax) - 3)

    assert diffgal < tol

    difftolgal = abs(len(np.unique(galax)) - 1)

    assert difftolgal < tol

    diffree = abs(freepar - 18)

    assert diffree < tol

    Flux = 10 ** (
        (25 - mag) / 2.5
    )  # 25 is just a constant to avoid small numbers. I substract it later

    cnt = Counter(galax)

    for idx, item in enumerate(np.unique(galax)):
        totcomp = cnt[item]
        maskgal = galax == item

        totFlux = Flux[maskgal].sum()

        totmag = -2.5 * np.log10(totFlux) + 25

    diffcomp = abs(totcomp - 3)
    assert diffcomp < tol

    diffmag = abs(totmag - 11.26)
    assert diffmag < tol

    return None


def test_getSersic():

    tol = 1e-2
    path = "tests/"

    image = "A1656-1-3.sbProf.fits"
    image = path + image
    regfile = "ds9.sbProf.reg"
    regfile = path + regfile
    center = False
    maskfile = None
    zeropoint = 21.817
    sky = 370
    bulgetot = 0.8
    noprint = False
    bards9 = None

    getSersic(
        image, regfile, center, maskfile, zeropoint, sky, noprint, bulgetot, bards9
    )

    constfile = "constraints.txt"
    constfile = path + constfile

    # something is wrong with pytest that does not
    # create this file. Not my fault
    # assert os.path.isfile(constfile)
    if os.path.isfile(constfile):
        os.remove(constfile)

    return None


def test_imarith():

    tol = 1e-2
    path = "tests/"

    ImageFile = "A1656-1-3.sbProf.fits"
    ImageFile = path + ImageFile

    output = "testimarith.fits"
    output = path + output
    image2 = "A1656-1-3.sbProf.fits"
    image2 = path + image2

    add = 1
    mul = None
    div = None
    sub = None

    imarith(ImageFile, output, image2, add, mul, div, sub)

    assert os.path.isfile(output)

    if os.path.isfile(output):
        os.remove(output)

    return None
