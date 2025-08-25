import pytest
import os


import subprocess as sp
import numpy as np

from galfitools.galout.getRads import getBreak
from galfitools.galout.getRads import getBreak2
from galfitools.galout.getRads import getFWHM
from galfitools.galout.getRads import getKappa
from galfitools.galout.getRads import getKappa2
from galfitools.galout.getRads import getReComp
from galfitools.galout.getRads import getSlope
from galfitools.galout.getRads import getBulgeRad
from galfitools.galout.getMissingLight import getMissLight
from galfitools.galout.getN import getN

# import galfitools.galout.getN as mod
import galfitools.galout.getN as getn_mod
from galfitools.galout.getMeRad import getMeRad
from galfitools.galout.getBT import getBT
from galfitools.galout.showcube import displayCube

from galfitools.galout.PhotDs9 import photDs9


from galfitools.galout.fitlog2csv import log2csv

from galfitools.galout.getCOW import getCOW
from galfitools.galout.getCOW import readGalfitF2 as mod
import galfitools.galout.getCOW as cow

from galfitools.galout.getPeak import getPeak

from galfitools.galout.getBarSize import getBarSize


def test_getBreak():

    galfitFile = "galfit.galout"

    path = "tests/"
    galfitFile = path + galfitFile

    dis = 10

    inicomp = 2

    quick = False
    random = False

    num_comp = 1

    angle = None

    ranx = None
    plot = False

    result = 248.62

    tol = 1e-2

    rbreak, N, theta = getBreak(
        galfitFile, dis, inicomp, quick, random, num_comp, angle, plot, ranx
    )

    diffrbreak = abs(rbreak - result)

    assert diffrbreak < tol

    return None


def test_getBreak2():

    galfitFile = "galfit.galout"
    path = "tests/"
    galfitFile = path + galfitFile

    dis = 10
    angle = None
    num_comp = 1
    plot = False
    ranx = None

    rbreak, N, theta = getBreak2(galfitFile, dis, angle, num_comp, plot, ranx)

    result = 42.7

    tol = 1e-2

    diffrbreak = abs(rbreak - result)

    assert diffrbreak < tol

    return None


def test_getFWHM():

    galfitFile = "galfit.galout"
    path = "tests/"
    galfitFile = path + galfitFile

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
    path = "tests/"
    galfitFile = path + galfitFile

    inicomp = 2
    quick = False
    random = False
    num_comp = 1
    angle = None
    ranx = None
    plot = False

    rkappa, N, theta = getKappa(
        galfitFile, dis, inicomp, quick, random, angle, num_comp, plot, ranx
    )

    result = 2.62
    tol = 1e-2

    diffrkappa = abs(rkappa - result)

    assert diffrkappa < tol

    return None


def test_getKappa2():

    galfitFile = "galfit.galout"
    dis = 10
    path = "tests/"
    galfitFile = path + galfitFile

    # inicomp = 2
    # quick = False
    # random = False
    num_comp = 1
    angle = None
    ranx = None
    plot = False

    # rkappa, N, theta = getKappa2(galfitFile, dis, inicomp, quick, random, angle, num_comp, plot, ranx)

    rkappa, N, theta = getKappa2(galfitFile, dis, angle, num_comp, plot, ranx)

    result = 2.62
    tol = 1e-1

    diffrkappa = abs(rkappa - result)

    assert diffrkappa < tol

    return None


def test_getReComp():

    galfitFile = "galfit.galout"
    path = "tests/"
    galfitFile = path + galfitFile

    dis = 10
    eff = 0.5
    num_comp = 1
    angle = None

    EffRad, totmag, meanme, me, N, theta = getReComp(
        galfitFile, dis, eff, angle, num_comp
    )

    result = 97.82
    tol = 1e-2

    diffEffRad = abs(EffRad - result)

    assert diffEffRad < tol

    return None


def test_getMeRad():

    galfitFile = "galfit.galout"
    path = "tests/"
    galfitFile = path + galfitFile

    dis = 5
    rad = 10
    num_comp = 1
    angle = None

    totmag, meanmerad, merad, N, theta = getMeRad(galfitFile, dis, rad, angle, num_comp)

    result = 19.16
    result2 = 19.66
    tol = 1e-2

    diffmeanme = abs(meanmerad - result)
    diffme = abs(merad - result2)

    assert diffmeanme < tol
    assert diffme < tol

    return None


def test_getBT():

    galfitFile = "galfit.barsize"
    path = "tests/"
    galfitFile = path + galfitFile

    dis = 5
    num_comp = 1

    bulge_total, totmag, N = getBT(galfitFile, dis, num_comp)

    result = 0.44
    tol = 1e-2

    diffbt = abs(bulge_total - result)

    assert diffbt < tol

    return None


def test_getSlope():

    galfitFile = "galfit.galout"
    dis = 10
    path = "tests/"
    galfitFile = path + galfitFile

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
    path = "tests/"
    galfitFile1 = path + galfitFile1
    galfitFile2 = path + galfitFile2

    dis = 10

    num_comp = 1

    angle = None

    ranx = None
    plot = False

    rbulge, N1, N2, theta = getBulgeRad(
        galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx
    )

    result = 3.86
    tol = 1e-2

    diffrbulge = abs(rbulge - result)

    assert diffrbulge < tol

    return None


def test_getMissLight():

    galfitFile1 = "galfit.1ser"
    galfitFile2 = "galfit.3ser"
    path = "tests/"
    galfitFile1 = path + galfitFile1
    galfitFile2 = path + galfitFile2

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
    path = "tests/"
    galfitFile = path + galfitFile

    num_comp = 1

    frac = 0.2  # irrelevant variable, it is not used anymore but it was kept

    angle = None

    plot = False

    sersic, meanser, stdser, totmag, N, theta = getN(
        galfitFile, dis, frac, angle, num_comp, plot
    )

    result1 = 3.48
    result2 = 2.66

    tol = 1e-2

    diffser = abs(sersic - result1)

    assert diffser < tol

    diffser2 = abs(meanser - result2)

    assert diffser2 < tol

    return None


def _fake_comps(active, posang):
    return type(
        "Comps",
        (),
        {
            "Active": np.array(active, dtype=int),
            "PosAng": np.array(posang, dtype=float),
        },
    )()


def test_getN_nominal(monkeypatch):
    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return type("H", (), {"scale": 0.5})()

        def ReadComps(self):
            return _fake_comps([1, 0, 1], [10.0, 20.0, 30.0])

    class FakeGetReff:
        def GetReSer(self, head, comps, frac, theta):
            return (10.0 if np.isclose(frac, 0.5) else 8.0), 12.34

        def GetRfracSer(self, head, comps, F, theta):
            return 20.0 * F

    class FakeGetMe:
        def MeanMe(self, totmag, r):
            return 20.0

        def Me(self, head, comps, r, theta):
            return 18.0

    class FakeGetN:
        def MeMeanMe(self, me, meanme):
            return 3.5

        def GalNs(self, Re, R, F):
            return np.full_like(F, 4.0, dtype=float)

    monkeypatch.setattr(getn_mod, "Galfit", FakeGalfit)
    monkeypatch.setattr(getn_mod, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(getn_mod, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(getn_mod, "numComps", lambda comps, mode: 3)
    monkeypatch.setattr(getn_mod, "GetReff", lambda: FakeGetReff())
    monkeypatch.setattr(getn_mod, "GetMe", lambda: FakeGetMe())
    monkeypatch.setattr(getn_mod, "GetN", lambda: FakeGetN())

    sersic, ns_mean, ns_std, totmag, N, theta = getn_mod.getN(
        "galfit.01", dis=10, frac=0.2, angle=None, num_comp=1, plot=False, const=0
    )
    assert sersic == 3.5
    assert ns_mean == 4.0
    assert ns_std == 0.0
    assert totmag == 12.34
    assert N == 3
    assert theta == 30.0


def test_getN_with_plot_calls_savefig(monkeypatch):
    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return type("H", (), {"scale": 1.0})()

        def ReadComps(self):
            return _fake_comps([1], [42.0])

    class FakeGetReff:
        def GetReSer(self, head, comps, frac, theta):
            return (5.0, 10.0)

        def GetRfracSer(self, head, comps, F, theta):
            return 10.0 * F

    class FakeGetMe:
        def MeanMe(self, totmag, r):
            return 1.0

        def Me(self, head, comps, r, theta):
            return 0.5

    class FakeGetN:
        def MeMeanMe(self, me, meanme):
            return 2.0

        def GalNs(self, Re, R, F):
            return np.ones_like(F) * 2.5

    # define FakePlt BEFORE patching
    calls = {"savefig": None}

    class FakePlt:
        def plot(self, *a, **k):
            pass

        def grid(self, *a, **k):
            pass

        def minorticks_on(self, *a, **k):
            pass

        def xlabel(self, *a, **k):
            pass

        def ylabel(self, *a, **k):
            pass

        def savefig(self, path):
            calls["savefig"] = path

    monkeypatch.setattr(getn_mod, "Galfit", FakeGalfit)
    monkeypatch.setattr(getn_mod, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(getn_mod, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(getn_mod, "numComps", lambda comps, mode: 1)
    monkeypatch.setattr(getn_mod, "GetReff", lambda: FakeGetReff())
    monkeypatch.setattr(getn_mod, "GetMe", lambda: FakeGetMe())
    monkeypatch.setattr(getn_mod, "GetN", lambda: FakeGetN())
    monkeypatch.setattr(getn_mod, "plt", FakePlt())  # only in this test

    sersic, ns_mean, ns_std, totmag, N, theta = getn_mod.getN(
        "galfit.01", dis=5, frac=0.3, angle=10.0, num_comp=1, plot=True, const=0.1
    )
    assert calls["savefig"] == "Serind.png"
    assert (sersic, ns_mean, ns_std, totmag, N, theta) == (2.0, 2.5, 0.0, 10.0, 1, 10.0)


def test_getN_no_components_exits(monkeypatch):
    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return type("H", (), {"scale": 1.0})()

        def ReadComps(self):
            return _fake_comps([0], [0.0])

    monkeypatch.setattr(getn_mod, "Galfit", FakeGalfit)
    monkeypatch.setattr(getn_mod, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(getn_mod, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(getn_mod, "numComps", lambda comps, mode: 0)

    with pytest.raises(SystemExit) as exc:
        getn_mod.getN(
            "galfit.01", dis=10, frac=0.2, angle=10.0, num_comp=1, plot=False, const=0.0
        )
    assert exc.value.code == 1


def test_displayCube():

    cubeimage = "A1656-showcube.fits"
    namecube = "cube.png"
    dpival = 100
    brightness = 0
    contrast = 1
    cmap = "viridis"
    scale = 1
    noplot = False

    path = "tests/"
    cubeimage = path + cubeimage
    namecube = path + namecube

    displayCube(cubeimage, namecube, dpival, brightness, contrast, cmap, scale, noplot)

    assert os.path.isfile(namecube)

    if os.path.isfile(namecube):
        os.remove(namecube)

    return None


def test_photDs9():

    ImageFile = "A671.gtMakeMask.maskds9.masksky.fits"
    RegFile = "maskds9.reg"
    path = "tests/"
    ImageFile = path + ImageFile
    RegFile = path + RegFile

    maskfile = "none"

    zeropoint = 25
    sky = 0

    mag, exptime = photDs9(ImageFile, RegFile, maskfile, zeropoint, sky)

    result = 11.71
    tol = 1e-2

    diffmag = abs(mag - result)

    assert diffmag < tol

    return None


def test_fitlog2csv():

    path = "tests/"
    file = "fit.log"
    file = path + file

    fileout = "fitlog.csv"

    num = None

    fileout = path + fileout

    log2csv(num, fileout, path=path)

    assert os.path.isfile(fileout)

    if os.path.isfile(fileout):
        os.remove(fileout)

    return None


def test_getCOW():

    path = "tests/"
    galfitFile = "galfit.1ser"
    galfitFile = path + galfitFile
    dis = 10
    plotfile = "cow.png"
    plotfile = path + plotfile

    dpival = 100

    frac = 0.95

    maxdiff = False

    num_comp = 1

    angle = None

    galfitF2 = None

    totmag, N, theta = getCOW(
        galfitFile, dis, angle, frac, num_comp, plotfile, dpival, galfitF2, maxdiff
    )

    assert os.path.isfile(plotfile)

    if os.path.isfile(plotfile):
        os.remove(plotfile)

    tol = 1e-2
    diff = abs(totmag - 11.96)
    assert diff < tol

    return None


def test_getPeak():

    path = "tests/"

    image = "A671.gtMakeMask.maskds9.masksky.fits"
    image = path + image

    regfile = "ds9.getStar.reg"
    regfile = path + regfile
    center = False
    maskfile = None

    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)

    tol = 1e-2

    diffx = abs(X - 1444)
    assert diffx < tol

    diffy = abs(Y - 1061)
    assert diffy < tol

    diffax = abs(AxRat - 1)
    assert diffax < tol

    diffPA = abs(PA - (-90.00))
    assert diffPA < tol

    return None


def test_getBarSize():

    path = "tests/"

    galfitFile = "galfit.barsize"
    galfitFile = path + galfitFile

    dis = 10
    out = "testbar.reg"
    out = path + out
    num_comp = 1
    plot = False
    ranx = None

    rbar, N, theta = getBarSize(galfitFile, dis, num_comp, plot, ranx, out)

    tol = 1e-2

    diffbar = abs(rbar - 131)
    assert diffbar < tol

    difftheta = abs(theta - 57.15)
    assert difftheta < tol

    assert os.path.isfile(out)

    if os.path.isfile(out):
        os.remove(out)

    return None


def test_readGalfitF2_with_angle(monkeypatch):
    """Angle provided → theta equals provided angle; arrays behave correctly."""

    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return {"scale": 0.5}

        def ReadComps(self):
            return type(
                "Comps",
                (),
                {
                    "Active": np.array([1, 1]),
                    "PosAng": np.array([10.0, 20.0]),
                },
            )()

    monkeypatch.setattr(cow, "Galfit", FakeGalfit)
    monkeypatch.setattr(cow, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(cow, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(cow, "numComps", lambda comps, mode: 2)

    head, comps, theta = cow.readGalfitF2("dummy.fits", dis=10, angle=45.0, num_comp=2)
    assert head["scale"] == 0.5
    assert theta == 45.0
    assert isinstance(comps.PosAng, np.ndarray)


def test_readGalfitF2_without_angle(monkeypatch):
    """Angle None → use last active component PosAng (boolean mask indexing works)."""

    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return {"scale": 1.0}

        def ReadComps(self):
            return type(
                "Comps",
                (),
                {
                    "Active": np.array([1, 0, 1]),
                    "PosAng": np.array([15.0, 25.0, 35.0]),
                },
            )()

    monkeypatch.setattr(cow, "Galfit", FakeGalfit)
    monkeypatch.setattr(cow, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(cow, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(cow, "numComps", lambda comps, mode: 3)

    head, comps, theta = cow.readGalfitF2("dummy2.fits", dis=5, angle=None, num_comp=3)
    assert head["scale"] == 1.0
    assert theta == 35.0  # last active PosAng


def test_readGalfitF2_no_components_exits(monkeypatch):
    """numComps == 0 → sys.exit(1). Provide angle to skip PosAng indexing."""

    class FakeGalfit:
        def __init__(self, f):
            self.f = f

        def ReadHead(self):
            return {"scale": 1.0}

        def ReadComps(self):
            return type(
                "Comps",
                (),
                {
                    "Active": np.array([0]),  # no active comps
                    "PosAng": np.array([0.0]),
                },
            )()

    monkeypatch.setattr(cow, "Galfit", FakeGalfit)
    monkeypatch.setattr(cow, "SelectGal", lambda comps, dis, num: comps)
    monkeypatch.setattr(cow, "conver2Sersic", lambda comps: comps)
    monkeypatch.setattr(cow, "numComps", lambda comps, mode: 0)

    with pytest.raises(SystemExit) as exc:
        # angle provided → avoid PosAng[mask][-1] on empty selection
        cow.readGalfitF2("dummy.fits", dis=10, angle=30.0, num_comp=1)
    assert exc.value.code == 1
