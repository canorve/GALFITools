import pytest
import os


import subprocess as sp


from galfitools.sky.GalfitSky import galfitSky

# from galfitools.sky.Sky import sky # deprecated

from galfitools.sky.SkyDs9 import SkyDs9
from galfitools.sky.SkyRing import SkyRing


def test_galfitSky():

    imgname = "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits"
    path = "tests/"
    imgname = path + imgname
    maskfile = path + maskfile
    mgzpt = 21.817
    scale = 0.68

    X = 169
    Y = 544

    initsky = 376

    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)

    assert os.path.isfile("sky.txt")

    if os.path.isfile("sky.txt"):
        os.remove("sky.txt")

    return None


def test_SkyDs9():

    ImageFile = "A671.gtMakeMask.maskds9.masksky.fits"
    RegFile = "skyDs9.reg"

    path = "tests/"
    ImageFile = path + ImageFile
    RegFile = path + RegFile
    maskfile = "none"

    mean, sig, ms, mstd = SkyDs9(
        ImageFile, RegFile, maskfile, outliers=True, mgzpt=21.471, scale=0.68
    )

    tol = 1e-2

    result1 = 1146.627
    result2 = 6.919
    result3 = 21.102
    result4 = 0.014

    diffsky1 = abs(mean - result1)

    assert diffsky1 < tol

    diffsky2 = abs(sig - result2)

    assert diffsky2 < tol

    diffsky3 = abs(ms - result3)

    assert diffsky3 < tol

    diffsky4 = abs(mstd - result4)

    assert diffsky4 < tol

    return None


def test_skyring():

    image = "A1656-1-3.sbProf.fits"
    ds9regfile = "skyring.reg"
    mask = "mask.sbProf.fits"

    path = "tests/"
    image = path + image
    ds9regfile = path + ds9regfile
    mask = path + mask

    width = 20
    center = False

    ##end input
    mean, std, median, ms, mstd, rad = SkyRing(
        image, mask, ds9regfile, width, center, outliers=True, mgzpt=21.817, scale=0.68
    )

    tol = 1e-2

    result1 = 370.415
    result2 = 4.429
    result3 = 21.48
    result4 = 0.01

    diffsky1 = abs(mean - result1)

    assert diffsky1 < tol

    diffsky2 = abs(std - result2)

    assert diffsky2 < tol

    diffsky3 = abs(ms - result3)

    assert diffsky3 < tol

    diffsky4 = abs(mstd - result4)

    assert diffsky4 < tol

    # assert os.path.isfile("ring.fits")

    if os.path.isfile("tests/skyring.fits"):
        os.remove("tests/skyring.fits")

    if os.path.isfile("skyring.fits"):
        os.remove("skyring.fits")

    # assert os.path.isfile("ringmask.fits")

    # if os.path.isfile("ringmask.fits"):
    #    os.remove("ringmask.fits")

    return None
