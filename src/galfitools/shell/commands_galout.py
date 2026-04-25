# commands_galout.py
import argparse
from galfitools.galin.galfit import Galfit
from galfitools.galout.fitlog2csv import log2csv
from galfitools.galout.getBarSize import getBarSize
from galfitools.galout.getCOW import getCOW
from galfitools.galout.getMissingLight import getMissLight
from galfitools.galout.getN import getN
from galfitools.galout.getBT import getBT
from galfitools.galout.magcorr import magCorr
from galfitools.galout.getPeak import getPeak
from galfitools.galout.getRads import (
    getBreak,
    getBreak2,
    getBulgeRad,
    getFWHM,
    getKappa,
    getKappa2,
    getReComp,
    getSlope,
)
from galfitools.galout.getMeRad import getMeRad
from galfitools.galout.PhotDs9 import photDs9
from galfitools.galout.showcube import displayCube
from galfitools.shell.prt import printWelcome

from galfitools.galout.getChiNu import getChiNu


def mainPhotDs9(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="computes photometry from a DS9 region file: Box, Ellipses and Polygons"
    )
    p.add_argument("ImageFile", help="fits image")
    p.add_argument("RegFile", help="DS9 region file")
    p.add_argument(
        "-zp", "--zeropoint", help="magnitude zero point", type=float, default=25
    )
    p.add_argument("-m", "--mask", help="mask image", type=str)
    p.add_argument("-sk", "--sky", help="sky background value", type=float, default=0)
    p.add_argument("-p", "--plate", help="plate scale", type=float, default=1)
    a = p.parse_args(argv)
    mag, sb, exptime = photDs9(
        a.ImageFile, a.RegFile, a.mask, a.zeropoint, a.plate, a.sky
    )
    print(f"the exposition time is: {exptime} \n")
    print(f"the magnitude and surface brightness are: \n")
    print(f"  mag   mag/arcsec^2 ")
    print(f"   {mag:.2f}     {sb:.2f}  \n")

    return 0


def mainGetBreak(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBreak: gets the break radius from a set of Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-ni", "--numinitial", type=int, default=2)
    p.add_argument("-q", "--quick", action="store_true")
    p.add_argument("-r", "--random", type=int)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rbreak, N, theta = getBreak(
        a.GalfitFile,
        a.dis,
        a.numinitial,
        a.quick,
        a.random,
        a.numcomp,
        a.angle,
        a.plot,
        a.ranx,
    )
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    plate = Galfit(a.GalfitFile).ReadHead().scale

    print(f"The break radius is:  \n")
    print(f" pixels   arcsec ")
    print(f" {rbreak:.2f}   {rbreak*plate:.2f}  \n")

    return 0


def mainGetBreak2(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBreak2: gets the break radius from a set of Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rbreak, N, theta = getBreak2(
        a.GalfitFile, a.dis, a.angle, a.numcomp, a.plot, a.ranx
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The break radius is:  \n")
    print(f" pixels   arcsec ")
    print(f" {rbreak:.2f}   {rbreak*plate:.2f}  \n")

    return 0


def mainGetBarSize(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBarlength: gets the bar length from Sersic and Ferrer models. Note: 3 components only: Bulge, Bar, Disk"
    )
    p.add_argument("GalfitFile", help="GALFIT file")
    p.add_argument(
        "-d", "--dis", type=int, default=3, help="maximum distance among components"
    )
    p.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=1,
        help="number of component to be selected",
    )
    p.add_argument(
        "-o", "--out", type=str, default="bar.reg", help="output DS9 ellipse region"
    )
    p.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="plots a file of kappa and break radius",
    )
    p.add_argument(
        "-r",
        "--red",
        action="store_true",
        help="If activated, DS9 region ellipse is red",
    )
    p.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="range of radius to search for barlength",
    )

    a = p.parse_args(argv)
    rbar, N, theta = getBarSize(
        a.GalfitFile, a.dis, a.numcomp, a.plot, a.ranx, a.out, a.red
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)

    print(f"The bar length is:  \n")
    print(f"  pixels   arcsec ")
    print(f"   {rbar:.2f}     {rbar*plate:.2f}  \n")

    return 0


def mainFWHM(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getFWHM: gets the FWHM from a set of Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    a = p.parse_args(argv)
    fwhm, N, theta = getFWHM(a.GalfitFile, a.dis, a.angle, a.numcomp)
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The FWHM is:  \n")
    print(f"  pixels   arcsec ")
    print(f"  {fwhm:.2f}   {fwhm*plate:.2f}  \n")

    return 0


def mainKappa(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getKappa: gets the Kappa radius from a set of Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-ni", "--numinitial", type=int, default=2)
    p.add_argument("-q", "--quick", action="store_true")
    p.add_argument("-r", "--random", type=int)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rkappa, N, theta = getKappa(
        a.GalfitFile,
        a.dis,
        a.numinitial,
        a.quick,
        a.random,
        a.angle,
        a.numcomp,
        a.plot,
        a.ranx,
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The Kappa radius is:  \n")
    print(f"  pixels   arcsec ")
    print(f"   {rkappa:.2f}     {rkappa*plate:.2f}  \n")

    return 0


def mainKappa2(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getKappa2: maximum curvature radius from Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rkappa, N, theta = getKappa2(
        a.GalfitFile, a.dis, a.angle, a.numcomp, a.plot, a.ranx
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The Kappa radius is:  \n")
    print(f"  pixels   arcsec ")
    print(f"   {rkappa:.2f}     {rkappa*plate:.2f}  \n")

    return 0


def mainGetReComp(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getReComp: gets the effective radius from Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=3)
    p.add_argument("-fr", "--fracrad", type=float, default=0.5)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument(
        "-mc", "--mecorr", help="Surface brightness correction", type=float, default=10
    )
    a = p.parse_args(argv)
    EffRad, totmag, meanme, me, N, theta = getReComp(
        a.GalfitFile, a.dis, a.fracrad, a.angle, a.numcomp, a.mecorr
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f": Total Magnitude of the galaxy: \n")
    print(f"   mag ")
    print(f"  {totmag:.2f}\n")

    print(f": Surface brightness at radius of {a.fracrad*100:.0f}% of light (μr): \n")
    print(f"   mag/''^2 ")
    print(f"   {me:.2f}\n")

    if a.fracrad == 0.5:
        print(f":  Mean Surface Brightness at effective radius (<μ>e)\n")
        print(f"   mag/''^2")
        print(f"   {meanme:.2f}  \n")

    print(f"The radius at {a.fracrad*100:.0f}% of light is:\n")
    print(f"  pixels   arcsec ")
    print(f"   {EffRad:.2f}     {EffRad*plate:.2f}  \n")

    return 0


def mainGetMeRad(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getMeRad: surface brightness at a given radius from Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=5)
    p.add_argument(
        "-r",
        "--rad",
        help="radius where the surface brightness is evaluated",
        type=float,
        default=10,
    )
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument(
        "-mc", "--mecorr", help="Surface brightness correction", type=float, default=10
    )
    a = p.parse_args(argv)
    totmag, meanmerad, merad, N, theta = getMeRad(
        a.GalfitFile, a.dis, a.rad, a.angle, a.numcomp, a.mecorr
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"Total Magnitude:  \n")
    print(f"  mag ")
    print(f"  {totmag:.2f}  \n")

    print(f"The radius used is:  \n")
    print(f"  pixels   arcsec ")
    print(f"   {a.rad:.2f}     {a.rad*plate:.2f}  \n")

    print(f"Surface brightness at this radius is (μr): \n")
    print(f"  mag/''^2 ")
    print(f"   {merad:.2f}  \n")

    print(f"Mean surface brightness at this radius is (<μ>r): \n")
    print(f"  mag/''^2 ")
    print(f"   {meanmerad:.2f}  \n")

    return 0


def maingetSlope(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getSlope: gets the slope radius from a set of Sersics"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-a", "--angle", type=float)
    p.add_argument("-s", "--slope", type=float, default=0.5)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rgam, N, theta = getSlope(
        a.GalfitFile, a.dis, a.slope, a.angle, a.numcomp, a.plot, a.ranx
    )
    plate = Galfit(a.GalfitFile).ReadHead().scale
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The radius with slope {a.slope:.2f} is:  \n")
    print(f"  pixels   arcsec ")
    print(f"   {rgam:.2f}     {rgam*plate:.2f}  \n")

    return 0


def maingetN(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getN: Sérsic index from surface brightness at effective radius"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-rf", "--radfrac", type=float, default=0.2)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-c", "--const", type=float, default=0)
    a = p.parse_args(argv)
    sersic, meanser, stdser, totmag, N, theta = getN(
        a.GalfitFile, a.dis, a.radfrac, a.angle, a.numcomp, a.plot, a.const
    )
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"Total Magnitude of the galaxy: \n")
    print(f"    mag     \n")
    print(f"   {totmag:.2f}   \n")

    print(f"Sersic index with the method of mean surface brightness at Re: \n")
    print(f"   {sersic:.2f}     {stdser:.2f}  \n")

    print(
        "Sersic index with the method of fraction of light\n(evaluated at different radius)\n"
    )

    print(f"Sersic index mean and standard deviation:\n")
    print(f"  mean   std")
    print(f"  {meanser:.2f}     {stdser:.2f}  \n")

    return 0


def maingetBT(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(description="getBT: Bulge to Total luminosity ratio")
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    a = p.parse_args(argv)
    bulge_total, totmag, N = getBT(a.GalfitFile, a.dis, a.numcomp)
    print("number of model components: ", N)
    print(f"Total Magnitude of the galaxy: \n")
    print(f"    mag     \n")
    print(f"   {totmag:.2f}   \n")

    print(f"Bulge to total luminosity ratio: \n")
    print(f"   {bulge_total:.2f}  \n")

    return 0


def mainGetBulgeRad(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBulgeRad: radius where two SB models are equal"
    )
    p.add_argument("GalfitFile1")
    p.add_argument("GalfitFile2")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-p", "--plot", action="store_true")
    p.add_argument("-rx", "--ranx", nargs=2, type=float)
    a = p.parse_args(argv)
    rbulge, N1, N2, theta = getBulgeRad(
        a.GalfitFile1, a.GalfitFile2, a.dis, a.numcomp, a.angle, a.plot, a.ranx
    )
    plate = Galfit(a.GalfitFile1).ReadHead().scale
    print("number of model components for the bulge: ", N1)
    print("number of model components for the rest of the galaxy: ", N2)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"The bulge radius is:\n")
    print(f"  pixels  arcsec ")
    print(f"   {rbulge:.2f}    {rbulge*plate:.2f}  \n")

    return 0


def mainMissingLight(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getMissLight: missing light between two SB models"
    )
    p.add_argument("GalfitFile1")
    p.add_argument("GalfitFile2")
    p.add_argument("rad", type=float)
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    a = p.parse_args(argv)
    magmiss, N1, N2 = getMissLight(
        a.GalfitFile1, a.GalfitFile2, a.dis, a.numcomp, a.rad
    )
    print("number of model components coreless model: ", N1)
    print("number of model components core model: ", N2)

    print(f"The missing light is:  \n")
    print(f"  mag ")
    print(f"   {magmiss:.2f}  \n")

    return 0


def mainShowCube(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(description="show the cube fits of the galfit output")
    p.add_argument("cubeimage")
    p.add_argument("-o", "--outimage", type=str, default="cube.png")
    p.add_argument("-br", "--brightness", type=float, default=0)
    p.add_argument("-co", "--contrast", type=float, default=1)
    p.add_argument("-cm", "--cmap", type=str, default="viridis")
    p.add_argument("-dpi", "--dotsinch", type=int, default=100)
    p.add_argument("-s", "--scale", type=float, default=1)
    p.add_argument("-np", "--noplot", action="store_true")
    a = p.parse_args(argv)
    displayCube(
        a.cubeimage,
        a.outimage,
        a.dotsinch,
        a.brightness,
        a.contrast,
        a.cmap,
        a.scale,
        a.noplot,
    )
    return 0


def maingetCOW(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getCOW: plot curve-of-growth from galfit.XX (Sersic only)"
    )
    p.add_argument("GalfitFile")
    p.add_argument("-d", "--dis", type=int, default=10)
    p.add_argument("-pf", "--plotfile", type=str, default="cow.png")
    p.add_argument("-g", "--galfitF2", type=str)
    p.add_argument("-md", "--maxdiff", action="store_true")
    p.add_argument("-fr", "--fracrad", type=float, default=0.95)
    p.add_argument("-n", "--numcomp", type=int, default=1)
    p.add_argument("-pa", "--angle", type=float)
    p.add_argument("-dpi", "--dotsinch", type=int, default=100)
    a = p.parse_args(argv)
    totmag, N, theta = getCOW(
        a.GalfitFile,
        a.dis,
        a.angle,
        a.fracrad,
        a.numcomp,
        a.plotfile,
        a.dotsinch,
        a.galfitF2,
        a.maxdiff,
    )
    print("number of model components: ", N)
    if a.angle:
        print(f"position angle: {theta:.2f} degrees\n")
    else:
        print(f"position angle: {theta:.2f} degrees, taken from component {N}\n")

    print(f"Total Magnitude of the galaxy: {totmag:.2f} \n")
    print("plot file: ", a.plotfile)
    return 0


def mainFitlog2CSV(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(description="converts fit.log file into a CSV file")
    p.add_argument("-o", "--fileout", type=str, default="fitlog.csv")
    p.add_argument("-n", "--num", type=int, help="fit number; default is last fit")
    p.add_argument("-p", "--path", type=str, help="path where fit.log is located")
    a = p.parse_args(argv)
    log2csv(a.num, a.fileout, path=a.path)
    print("conversion to CSV file done.")
    return 0


def maingetPeak(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="Obtains center, axis ratio and PA from DS9 region"
    )
    p.add_argument("Image")
    p.add_argument("RegFile")
    p.add_argument("-c", "--center", action="store_true")
    p.add_argument("-m", "--mask", type=str)
    a = p.parse_args(argv)
    X, Y, AxRat, PA = getPeak(a.Image, a.RegFile, a.center, a.mask)
    print("peak is at (x, y) = ")
    print(f"  {X:.2f}  {Y:.2f}  \n")

    print("axis ratio q = {:.2f} ".format(AxRat))
    print("position angle = {:.2f}".format(PA))
    return 0


def mainmagCorr(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="corrects GALFIT file for magnitude correction and K corrections"
    )
    parser.add_argument("galfile", help="GALFIT input file")
    parser.add_argument(
        "-o", "--out", default="magcorr.gal", type=str, help="output GALFIT file"
    )

    parser.add_argument(
        "-a",
        "--Aext",
        default=0,
        type=float,
        help="Magnitude correction for Extinction",
    )

    parser.add_argument("-k", "--Kcor", type=float, default=0, help="K correction")

    args = parser.parse_args(argv)

    magCorr(args.galfile, args.out, A=args.Aext, K=args.Kcor)

    return 0


def maingetChiNu(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(description="Computes the Chinu within a radius")
    parser.add_argument("galfile", help="GALFIT input file")

    parser.add_argument(
        "-fr",
        "--fracrad",
        default=0.99,
        type=float,
        help="Fraction of light radius where the chinu will be computed. Computed from Sersic functions only. Default = 0.99",
    )

    parser.add_argument(
        "-r",
        "--RegFile",
        help=" DS9 ellipse region file to use instead of fracrad. Must be based on the GALFIT output image cube.",
        type=str,
    )

    parser.add_argument(
        "-nc",
        "--numcomp",
        type=int,
        default=1,
        help="Number of component inside galfit file. Default=1",
    )

    args = parser.parse_args(argv)

    chinu, aic, bic, numparfree = getChiNu(
        args.galfile, args.numcomp, args.fracrad, args.RegFile
    )

    print(f"Number of free parameters:")
    print(f"  {numparfree}")

    print(f"Reduced chi-square within the {args.fracrad*100:.2f} % light radius:")
    print(f"  {chinu:.2f}")

    print(f"Akaike information criterion:")
    print(f"  {aic:.2f}")

    print(f"Bayesian information criterion:")
    print(f"  {bic:.2f}")

    return 0
