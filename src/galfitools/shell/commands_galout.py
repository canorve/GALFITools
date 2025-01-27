import argparse

from galfitools.galin.galfit import Galfit
from galfitools.galout.fitlog2csv import log2csv
from galfitools.galout.getBarSize import getBarSize
from galfitools.galout.getCOW import getCOW
from galfitools.galout.getMissingLight import getMissLight
from galfitools.galout.getN import getN
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
from galfitools.galout.PhotDs9 import photDs9
from galfitools.galout.showcube import displayCube
from galfitools.shell.prt import printWelcome


def mainPhotDs9():
    """
    Calls the photDs9 function based on argument parsing.
    This function serves as an example of an API for photDs9.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="computes photometry from a Ds9 region "
        + "file: Box, Ellipses and Polygons "
    )

    parser.add_argument(
        "ImageFile", help="the image file where the photometry will be computed"
    )
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument(
        "-zp",
        "--zeropoint",
        type=float,
        help="The value of the zero point. Default = 25 ",
        default=25,
    )

    parser.add_argument("-m", "--mask", type=str, help="the mask file")

    parser.add_argument(
        "-sk",
        "--sky",
        default=0,
        type=float,
        help="Sky background value to be removed from image "
        + "before photometry. Default = 0",
    )

    args = parser.parse_args()

    ImageFile = args.ImageFile
    RegFile = args.RegFile
    maskfile = args.mask
    zeropoint = args.zeropoint
    sky = args.sky

    mag, exptime = photDs9(ImageFile, RegFile, maskfile, zeropoint, sky)

    line = "the exposition time is is: {} \n".format(exptime)
    print(line)

    line = "the magnitude from the ds9 region is: {:.2f} \n".format(mag)
    print(line)


# console scripts
def mainGetBreak() -> None:
    """
    Calls the getBreak function based on argument parsing.
    This function serves as an example of an API for getBreak.
    """

    printWelcome()
    # reading argument parsing

    parser = argparse.ArgumentParser(
        description="getBreak: gets the break radius from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of "
        + "all components, default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Position angle of the major axis of the galaxy. Default ="
        + " it will take the angle of the last components",
    )

    parser.add_argument(
        "-ni",
        "--numinitial",
        type=int,
        help="Number of component where it'll obtain the initial parameter "
        + "to search break radius or to generated random initial radius. ",
        default=2,
    )

    parser.add_argument(
        "-q",
        "--quick",
        action="store_true",
        help="evaluate in position only (given by -ni parameter)",
    )

    parser.add_argument(
        "-r",
        "--random",
        type=int,
        help="Number of random radius as initial parameters to search for the "
        + "minimum. It will generated random radius from 0 to effective radius"
        + " of the component indicated by parameter -ni ",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of double derivative vs. radius ",
    )
    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for the plot x-axis: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    inicomp = args.numinitial

    quick = args.quick
    random = args.random

    num_comp = args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot

    rbreak, N, theta = getBreak(
        galfitFile, dis, inicomp, quick, random, num_comp, angle, plot, ranx
    )

    print("number of model components: ", N)

    line = "Using a theta value of: {:.2f} degrees\n".format(theta)
    print(line)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rbreak_arc = rbreak * plate

    line = 'The break radius is {:.2f} pixels or {:.2f} "  \n'.format(
        rbreak, rbreak_arc
    )
    print(line)


def mainGetBreak2():
    """
    Calls the getBreak2 function based on argument parsing.
    This function serves as an example of an API for getBreak2.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="getBreak2: gets the break radius from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )
    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Position angle of the major axis of the galaxy. Default = it will"
        + " take the angle of the last components",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of double derivative vs. radius ",
    )

    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="x-axis range to search for the Break radius: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    angle = args.angle
    num_comp = args.numcomp
    plot = args.plot
    ranx = args.ranx

    rbreak, N, theta = getBreak2(galfitFile, dis, angle, num_comp, plot, ranx)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rbreak_arc = rbreak * plate

    print("number of model components: ", N)

    line = "Using a theta value of: {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The break radius is {:.2f} pixels or {:.2f} " \n'.format(rbreak, rbreak_arc)
    print(line)


def mainGetBarSize():
    """
    Calls the getBarSize function based on argument parsing.
    This function serves as an example of an API for getBarSize.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="getBarSize: gets the bar size from Sersic models: Bulge,"
        + " Bar and disk. It assumes that bar is galfit component 2 "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile",
        help="Galfit File containing the Sersics components bulge, bar, disk",
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )
    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-o",
        "--out",
        type=str,
        default="bar.reg",
        help="output DS9 ellipse region file",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of double derivatives and Kappa radius ",
    )

    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="x-axis range to search for the Break and Kappa radius: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    out = args.out
    num_comp = args.numcomp
    plot = args.plot
    ranx = args.ranx

    rbar, N, theta = getBarSize(galfitFile, dis, num_comp, plot, ranx, out)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rbar_arc = rbar * plate

    print("number of model components: ", N)

    line = "Using a theta value of: {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The bar size is {:.2f} pixels or {:.2f} " \n'.format(rbar, rbar_arc)
    print(line)


def mainFWHM() -> None:
    """
    Calls the getFWHM function based on argument parsing.
    This function serves as an example of an API for getFWHM.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="getFWHM: gets the FWHM from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Position angle of the major axis of the galaxy. Default ="
        + " it will take the angle of the last components",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    angle = args.angle

    num_comp = args.numcomp

    fwhm, N, theta = getFWHM(galfitFile, dis, angle, num_comp)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    fwhm_arc = fwhm * plate

    print("number of model components: ", N)

    line = "Using a theta value of : {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The FWHM is {:.2f} pixels or {:.2f} " \n'.format(fwhm, fwhm_arc)
    print(line)


# console scripts
def mainKappa() -> None:
    """
    Calls the getKappa function based on argument parsing.
    This function serves as an example of an API for getKappa.
    """

    printWelcome()

    # reading argument parsing

    parser = argparse.ArgumentParser(
        description="getKappa: gets the Kappa radius from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Position angle of the major axis of the galaxy. Default ="
        + " it will take the angle of the last components",
    )

    parser.add_argument(
        "-ni",
        "--numinitial",
        type=int,
        help="Number of component where it'll obtain the initial parameter"
        + " to search break radius or to generated random initial radius. ",
        default=2,
    )

    parser.add_argument(
        "-q",
        "--quick",
        action="store_true",
        help="evaluate in position only (given by -ni parameter",
    )

    parser.add_argument(
        "-r",
        "--random",
        type=int,
        help="Number of random radius as initial parameters to search "
        + "for the minimum. It will generated random radius from 0 to "
        + "effective radius of the component indicated by parameter -ni ",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of maximum curvature vs. radius ",
    )
    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for x-axis: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    inicomp = args.numinitial

    quick = args.quick
    random = args.random
    angle = args.angle

    ranx = args.ranx
    plot = args.plot

    num_comp = args.numcomp

    rkappa, N, theta = getKappa(
        galfitFile, dis, inicomp, quick, random, angle, num_comp, plot, ranx
    )

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rkappa_arc = rkappa * plate

    print("number of model components: ", N)

    line = 'The Kappa radius  is {:.2f} pixels or {:.2f} " \n'.format(
        rkappa, rkappa_arc
    )
    print(line)


def mainKappa2():
    """
    Calls the getKappa2 function based on argument parsing.
    This function serves as an example of an API for getKappa2.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="getKappa2: gets the Kappa radius (maximum curvature radius)"
        + " from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )
    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Position angle of the major axis of the galaxy. Default = "
        + "it will take the angle of the last components",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of double derivative vs. radius ",
    )

    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="x-axis range to search for the Break radius: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    angle = args.angle
    num_comp = args.numcomp
    plot = args.plot
    ranx = args.ranx

    rkappa, N, theta = getKappa2(galfitFile, dis, angle, num_comp, plot, ranx)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rkappa_arc = rkappa * plate

    print("number of model components: ", N)

    line = "Using a theta value of: {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The kappa radius is {:.2f} pixels or {:.2f} " \n'.format(rkappa, rkappa_arc)
    print(line)


def mainGetReComp() -> None:
    """
    Calls the getReComp function based on argument parsing.
    This function serves as an example of an API for getReComp.
    """

    printWelcome()
    # reading arguments parsing
    parser = argparse.ArgumentParser(
        description="getReComp: gets the effective radius from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )
    parser.add_argument(
        "-fr",
        "--fracrad",
        type=float,
        help="fraction of light radius. default=.5 for effective radius ",
        default=0.5,
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Angle of the major axis of the galaxy. Default = it will take the "
        + " angle of the last components. Angle measured from Y-Axis as"
        + "same as GALFIT. ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    eff = args.fracrad
    num_comp = args.numcomp
    angle = args.angle

    EffRad, totmag, meanme, me, N, theta = getReComp(
        galfitFile, dis, eff, angle, num_comp
    )

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    EffRad_arc = EffRad * plate

    print("number of model components: ", N)

    line = "Using a theta value of : {:.2f} degrees \n".format(theta)
    print(line)

    line = "Total Magnitude of the galaxy: {:.2f} \n".format(totmag)
    print(line)

    line1 = "Surface brightness at radius of "
    line2 = "{:.0f}% of light (\u03BCr): {:.2f} mag/'' \n".format(eff * 100, me)
    line = line1 + line2
    print(line)

    line1 = "Mean Surface Brightness at effective radius (<\u03BC>e):"
    line2 = " {:.2f} mag/'' \n".format(meanme)
    line = line1 + line2
    if eff == 0.5:
        print(line)

    line = 'The radius at {:.0f}% of light is {:.2f} pixels or {:.2f} " \n'.format(
        eff * 100, EffRad, EffRad_arc
    )
    print(line)


def maingetSlope() -> None:
    """
    Calls the getSlope function based on argument parsing.
    This function serves as an example of an API for getSlope.
    """

    printWelcome()
    # reading argument parsing

    parser = argparse.ArgumentParser(
        description="getSlope: gets the slope radius from a set of Sersics "
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of "
        + "all components, default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-a",
        "--angle",
        type=float,
        help="Angle of the major axis of the galaxy. Default = "
        + "it will take the angle of the last components",
    )

    parser.add_argument(
        "-s",
        "--slope",
        type=float,
        help="value of slope to find. default=.5 ",
        default=0.5,
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of first derivative vs. radius ",
    )
    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for x-axis: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    slope = args.slope

    num_comp = args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot

    rgam, N, theta = getSlope(galfitFile, dis, slope, angle, num_comp, plot, ranx)

    galfit = Galfit(galfitFile)

    head = galfit.ReadHead()

    plate = head.scale
    rgam_arc = rgam * plate

    print("number of model components: ", N)

    line = "Using a theta value of : {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The radius with slope {:.2f} is {:.2f} pixels, or {:.2f} " \n'.format(
        slope, rgam, rgam_arc
    )
    print(line)


def maingetN() -> None:
    """
    Calls the getN function based on argument parsing.
    This function serves as an example of an API for getN.
    """

    printWelcome()
    # reading arguments parsing

    parser = argparse.ArgumentParser(
        description="getN: computes the Sersic index from surface brightness"
        + " at effective radius"
    )

    # required arguments
    parser.add_argument(
        "GalfitFile", help="Galfit File containing the Sersics or gaussians components"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components, "
        + "default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Angle of the major axis of the galaxy. Default = it will take "
        + "the angle of the last components. Angle measured from Y-Axis as "
        + "same as GALFIT. ",
    )

    parser.add_argument(
        "-rf",
        "--radfrac",
        type=float,
        help="fraction of light radius. Default = .2 ",
        default=0.2,
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of Sersic index vs. fraction of light ",
    )

    parser.add_argument(
        "-c",
        "--const",
        type=float,
        help="Substract constant from plot. Default = 0",
        default=0,
    )

    # parsing variables

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    num_comp = args.numcomp

    frac = args.radfrac

    angle = args.angle

    plot = args.plot
    const = args.const

    sersic, meanser, stdser, totmag, N, theta = getN(
        galfitFile, dis, frac, angle, num_comp, plot, const
    )

    print("number of model components: ", N)

    line = "Using a theta value of : {:.2f} degrees \n".format(theta)
    print(line)

    line = "Total Magnitude of the galaxy: {:.2f} \n".format(totmag)
    print(line)

    line = "Sersic index with the method of Mean Surface "
    line2 = "Brightness at effective radius: {:.2f} \n".format(sersic)
    line = line + line2
    print(line)

    line = "Sersic index with the method of fraction of light\n"
    line2 = "(evaluated at different radius)\n"
    line = line + line2
    print(line)

    line = "Sersic index mean: {:.2f}  Standard deviation: {:.2f}  ".format(
        meanser, stdser
    )
    print(line)


def mainGetBulgeRad() -> None:
    """
    Calls the getBulgeRad function based on argument parsing.
    This function serves as an example of an API for getBulgeRad.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="getBulgeRad: gets the bulge radius or the radius where "
        + "two models of surface brightness models are equal"
    )

    # required arguments
    parser.add_argument(
        "GalfitFile1",
        help="Galfit File containing the coreless surface brightness model",
    )
    parser.add_argument(
        "GalfitFile2", help="Galfit File containing the core surface brightness model"
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all components,"
        + " default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Angle of the major axis of the galaxy. Default = "
        + "it will take the angle of the last components. Angle "
        + "measured from Y-Axis as same as GALFIT. ",
    )

    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="makes plot of GalfitFile1 - GalfitFile2 vs. radius ",
    )
    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for x-axis: xmin - xmax ",
    )

    args = parser.parse_args()

    galfitFile1 = args.GalfitFile1
    galfitFile2 = args.GalfitFile2

    dis = args.dis

    num_comp = args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot

    rbulge, N1, N2, theta = getBulgeRad(
        galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx
    )

    galfit = Galfit(galfitFile1)

    head = galfit.ReadHead()

    plate = head.scale
    rbulge_arc = rbulge * plate

    print("number of model components for the bulge: ", N1)
    print("number of model components for the rest of the galaxy: ", N2)

    line = "Using a theta value of: {:.2f} degrees\n".format(theta)
    print(line)

    line = 'The bulge radius is {:.2f} pixels or {:.2f} "  \n'.format(
        rbulge, rbulge_arc
    )
    print(line)


def mainMissingLight() -> None:
    """
    Calls the getMissLight function based on argument parsing.
    This function serves as an example of an API for getMissLight.
    """

    printWelcome()
    # reading arguments parsing

    parser = argparse.ArgumentParser(
        description="getMissLight: computes the missing light from two "
        + "surface brightness models"
    )

    # required arguments
    parser.add_argument(
        "GalfitFile1",
        help="Galfit File containing the coreless surface brightness model",
    )
    parser.add_argument(
        "GalfitFile2", help="Galfit File containing the core surface brightness model"
    )

    parser.add_argument(
        "rad",
        type=float,
        help="upper limit of radius to integrate the missing light in pixels",
    )

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all "
        + "components, default = 1 ",
        default=1,
    )

    # parsing variables

    args = parser.parse_args()

    galfitFile1 = args.GalfitFile1
    galfitFile2 = args.GalfitFile2

    dis = args.dis

    num_comp = args.numcomp

    rad = args.rad

    magmiss, N1, N2 = getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad)

    print("number of model components coreless model: ", N1)
    print("number of model components core model: ", N2)

    line = "the missing light is {:.2f} mag \n".format(magmiss)
    print(line)


def mainShowCube():
    """
    Calls the displayCube function based on argument parsing.
    This function serves as an example of an API for displayCube.
    """

    printWelcome()
    parser = argparse.ArgumentParser(
        description="show the cube fits of the galfit output"
    )

    parser.add_argument("cubeimage", help="the cube GALFIT image")
    parser.add_argument(
        "-o", "--outimage", type=str, help="the output png file", default="cube.png"
    )

    ####
    parser.add_argument(
        "-br",
        "--brightness",
        type=float,
        help="brightness of the image. Only for galaxy and model. "
        + "Default = 0. Preferible range goes from -1 to 1",
        default=0,
    )
    parser.add_argument(
        "-co",
        "--contrast",
        type=float,
        help="contrast of the image. Only for galaxy and model. Default = 1. "
        + "Preferible range goes from 0 to 1",
        default=1,
    )

    parser.add_argument(
        "-cm",
        "--cmap",
        type=str,
        help="cmap to be used for the cube image ",
        default="viridis",
    )
    parser.add_argument(
        "-dpi",
        "--dotsinch",
        type=int,
        help="dots per inch used for images files ",
        default=100,
    )
    parser.add_argument(
        "-s",
        "--scale",
        type=float,
        help="plate scale of the image. Default = 1",
        default=1,
    )
    parser.add_argument(
        "-np", "--noplot", action="store_true", help="it doesn't show plotting window"
    )

    args = parser.parse_args()

    cubeimage = args.cubeimage
    namecube = args.outimage
    dpival = args.dotsinch
    brightness = args.brightness
    contrast = args.contrast
    cmap = args.cmap
    scale = args.scale
    noplot = args.noplot

    displayCube(cubeimage, namecube, dpival, brightness, contrast, cmap, scale, noplot)


def maingetCOW() -> None:
    """
    Calls the getStar function based on argument parsing.
    This function serves as an example of an API for getStar.
    """

    printWelcome()
    # reading arguments parsing

    parser = argparse.ArgumentParser(
        description="getCOW: plots the curve-of-growth from the galfit.XX "
        + "file. Only for Sersic functions"
    )

    # required arguments
    parser.add_argument("GalfitFile", help="GALFIT File containing the Sersics")

    parser.add_argument(
        "-d", "--dis", type=int, help="Maximum distance among components", default=10
    )

    parser.add_argument(
        "-pf", "--plotfile", type=str, help="name of the plot file", default="cow.png"
    )

    parser.add_argument(
        "-g",
        "--galfitF2",
        type=str,
        help="Second GALFIT file to add to the plot (optional) ",
    )

    parser.add_argument(
        "-md",
        "--maxdiff",
        action="store_true",
        help="plot the maximum difference between model 1 and 2 (a vertical line) ",
    )

    parser.add_argument(
        "-fr",
        "--fracrad",
        type=float,
        help="fraction of light radius. This is the upper limit of "
        + "X-Axis. default=.95 ",
        default=0.95,
    )

    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        help="Number of component where it'll obtain center of all "
        + "components, default = 1 ",
        default=1,
    )

    parser.add_argument(
        "-pa",
        "--angle",
        type=float,
        help="Angle of the major axis of the galaxy. Default = "
        + "it will take the angle of the last components. Angle "
        + "measured from Y-Axis as same as GALFIT. ",
    )

    parser.add_argument(
        "-dpi",
        "--dotsinch",
        type=int,
        help="dots per inch used for images files ",
        default=100,
    )

    # parsing variables

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    plotfile = args.plotfile

    dpival = args.dotsinch

    frac = args.fracrad

    maxdiff = args.maxdiff

    num_comp = args.numcomp

    angle = args.angle

    galfitF2 = args.galfitF2

    totmag, N, theta = getCOW(
        galfitFile, dis, angle, frac, num_comp, plotfile, dpival, galfitF2, maxdiff
    )

    print("number of model components: ", N)

    line = "Using a theta value of : {:.2f} degrees \n".format(theta)
    print(line)

    line = "Total Magnitude of the galaxy: {:.2f} \n".format(totmag)
    print(line)

    print("plot file: ", plotfile)


def mainFitlog2CSV():
    """
    Calls the getStar function based on argument parsing.
    This function serves as an example of an API for getStar.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="converts fit.log file into a comma separated values file"
    )

    parser.add_argument(
        "-o", "--fileout", type=str, help="output file", default="fitlog.csv"
    )

    parser.add_argument(
        "-n",
        "--num",
        type=int,
        help="number of the fit in fit.log to convert to csv. Default: last fit",
    )

    parser.add_argument("-p", "--path", type=str, help="path where fit.log is located")

    args = parser.parse_args()

    fileout = args.fileout
    num = args.num
    path = args.path

    log2csv(num, fileout, path=path)

    print("conversion to CSV file done.")


def maingetPeak():
    """
    Calls the getStar function based on argument parsing.
    This function serves as an example of an API for getStar.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="Obtains the center, axis ratio and angular "
        + "position from DS9 region"
    )

    parser.add_argument("Image", help="image fits file")
    parser.add_argument("RegFile", help="DS9 ellipse region file")

    parser.add_argument(
        "-c", "--center", action="store_true", help="takes center of ds9 region file"
    )
    parser.add_argument("-m", "--mask", type=str, help="the mask file")

    args = parser.parse_args()

    image = args.Image
    regfile = args.RegFile
    center = args.center
    maskfile = args.mask

    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)

    print("peak is at (x, y) = ", X, Y)

    print("axis ratio q = {:.2f} ".format(AxRat))
    print("angular position = {:.2f}".format(PA))
