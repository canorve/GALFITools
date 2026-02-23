import argparse
from collections import Counter
import numpy as np

from galfitools.galin.checkGalFile import checkFile
from galfitools.galin.getBoxSizeDs9 import getBoxSizeDs9
from galfitools.galin.getSersic import getSersic
from galfitools.galin.getStar import getStar
from galfitools.galin.imarith import imarith
from galfitools.galin.initgal import InitGal
from galfitools.galin.MakeMask import makeMask
from galfitools.galin.MakePSF import makePSF
from galfitools.galin.MaskDs9 import maskDs9
from galfitools.galin.MaskSky import skyRem
from galfitools.galin.xy2fits import xy2fits
from galfitools.shell.prt import printWelcome
from galfitools.galin.sersic2ferrer import Sersic2Ferrer


def mainGetStar(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="gets an image slice centered on the object peak"
    )
    parser.add_argument("image", help="the image file to obtain the slice")
    parser.add_argument(
        "Ds9regFile", help="the DS9 ellipse region file containing the star"
    )
    parser.add_argument("size", type=int, help="the size of the new image in pixels")
    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="use the center in DS9 region; else find x,y peak within DS9 ellipse",
    )
    parser.add_argument(
        "-s", "--sky", type=float, help="the sky background to be removed"
    )
    parser.add_argument(
        "-o", "--out", type=str, help="the image output.", default="star.fits"
    )
    parser.add_argument("-sig", "--sigma", type=str, help="introduce the sigma image")
    parser.add_argument(
        "-so",
        "--sigout",
        type=str,
        help="the sigma image output.",
        default="sigma.fits",
    )
    args = parser.parse_args(argv)

    getStar(
        args.image,
        args.Ds9regFile,
        args.size,
        args.center,
        args.sky,
        args.out,
        args.sigma,
        args.sigout,
    )
    print(f"Done. Object fits file: {args.out} created ")
    if args.sigma:  # pragma: no cover
        print(f"Done. sigma fits file: {args.sigout} created ")
    return 0


def mainInitGal(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="Creates GALFIT input files with different initial parameters"
    )
    parser.add_argument("inFile", help="the galfit file galfit.XX")
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=1,
        help="number of files generated. Default=1",
    )
    parser.add_argument(
        "-p3", "--param3", nargs=2, type=float, help="values for parameter 3) [min max]"
    )
    parser.add_argument(
        "-p4", "--param4", nargs=2, type=float, help="values for parameter 4) [min max]"
    )
    parser.add_argument(
        "-p5", "--param5", nargs=2, type=float, help="values for parameter 5) [min max]"
    )
    parser.add_argument(
        "-p6", "--param6", nargs=2, type=float, help="values for parameter 6) [min max]"
    )
    parser.add_argument(
        "-p7", "--param7", nargs=2, type=float, help="values for parameter 7) [min max]"
    )
    parser.add_argument(
        "-p8", "--param8", nargs=2, type=float, help="values for parameter 8) [min max]"
    )
    parser.add_argument(
        "-p9", "--param9", nargs=2, type=float, help="values for parameter 9) [min max]"
    )
    parser.add_argument(
        "-p10",
        "--param10",
        nargs=2,
        type=float,
        help="values for parameter 10) [min max]",
    )
    parser.add_argument(
        "-nc",
        "--numcomp",
        type=int,
        help="component number whose parameters will be changed",
    )
    args = parser.parse_args(argv)

    InitGal(
        args.inFile,
        args.number,
        args.param3,
        args.param4,
        args.param5,
        args.param6,
        args.param7,
        args.param8,
        args.param9,
        args.param10,
        args.numcomp,
    )
    print("rungalfit.sh has been created")
    return 0


def mainMakeMask(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="creates mask file from a SExtractor catalog"
    )
    parser.add_argument("Sexfile", help="SExtractor catalog file")
    parser.add_argument("ImageFile", help="Image file")
    parser.add_argument(
        "-o",
        "--maskout",
        type=str,
        default="masksex.fits",
        help="output mask file name",
    )
    parser.add_argument(
        "-sf", "--satds9", type=str, default="ds9sat.reg", help="ds9 saturation file"
    )
    parser.add_argument(
        "-s", "--scale", type=float, default=1, help="scale factor for ellipses"
    )
    args = parser.parse_args(argv)

    makeMask(args.Sexfile, args.ImageFile, args.maskout, args.scale, args.satds9)
    print("Done. Mask image created ")
    return 0


def mainMaskDs9(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="creates/modifies a mask image for GALFIT using DS9 regions"
    )
    parser.add_argument("MaskFile", help="Mask image file to modify or create")
    parser.add_argument("RegFile", help="DS9 region file")
    parser.add_argument(
        "-f", "--fill", type=int, default=0, help="value to fill DS9 regions (0=remove)"
    )
    parser.add_argument("-i", "--image", type=str, help="image to obtain the size")
    parser.add_argument(
        "-b",
        "--border",
        action="store_true",
        help="mask borders when their value is zero",
    )
    parser.add_argument(
        "-bv",
        "--borValue",
        type=float,
        default=0,
        help="border value if different from zero",
    )
    parser.add_argument("-sm", "--skymean", type=float, help="sky mean for sky patch")
    parser.add_argument(
        "-sd", "--skystd", type=float, default=1, help="sky std for sky patch"
    )
    args = parser.parse_args(argv)

    maskDs9(
        args.MaskFile,
        args.RegFile,
        args.fill,
        args.image,
        args.border,
        args.borValue,
        skymean=args.skymean,
        skystd=args.skystd,
    )
    print(f"Done. Mask {args.MaskFile} created (or modified)")
    return 0


def mainMaskSky(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="creates a mask using image and sky mean/sigma"
    )
    parser.add_argument("ImageFile", help="original data image")
    parser.add_argument("MaskFile", help="name of the new Mask file")
    parser.add_argument(
        "-sm", "--skymean", type=float, default=0, help="sky background mean"
    )
    parser.add_argument(
        "-ss", "--skysigma", type=float, default=0, help="sky background sigma"
    )
    parser.add_argument(
        "-ns",
        "--numbersig",
        type=float,
        default=1,
        help="multiplier for sigma to remove sky background",
    )
    parser.add_argument(
        "-r", "--region", type=str, help="DS9 region to remove from mask"
    )
    parser.add_argument(
        "-b",
        "--border",
        action="store_true",
        help="mask borders when their value is zero",
    )
    parser.add_argument(
        "-bv",
        "--borValue",
        type=float,
        default=0,
        help="border value if different from zero",
    )
    args = parser.parse_args(argv)

    skyRem(
        args.ImageFile,
        args.MaskFile,
        args.skymean,
        args.skysigma,
        args.numbersig,
        args.borValue,
        args.border,
    )
    if args.region:  # pragma: no cover
        maskDs9(
            args.MaskFile, args.region, 0, None, False, 0
        )  # remove ds9 region if provided
    print(f"Done. Mask file {args.MaskFile} created")
    return 0


def mainxy2fits(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="convert ASCII x,y positions to FITS mask"
    )
    parser.add_argument("ImageFile", help="The Image file")
    parser.add_argument("AsciiMask", help="The ASCII file with the x,y positions")
    parser.add_argument(
        "-c", "--val", type=int, default=1, help="value for masked pixels (counts)"
    )
    args = parser.parse_args(argv)

    xy2fits().MakeFits(args.ImageFile, args.AsciiMask, args.val)
    print("Ascii -> Fits done ")
    return 0


def maincheckFile(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="check parameters and file names inside a GALFIT input file"
    )
    parser.add_argument("GalfitFile", help="GALFIT input file")
    parser.add_argument(
        "-d",
        "--dis",
        type=int,
        default=10,
        help="maximum distance in pixels among components (Default=10)",
    )
    args = parser.parse_args(argv)

    headinfo, galax, mag, freepar = checkFile(args.GalfitFile, args.dis)

    testflag = False
    if not (headinfo.inputimageflag):  # pragma: no cover
        print(f"File {headinfo.inputimage} not found ")
        testflag = True
    if not (headinfo.sigimageflag):  # pragma: no cover
        print(f"File {headinfo.sigimage} not found ")
        print("GALFIT will create a sigma file ")
    if not (headinfo.psfimageflag):  # pragma: no cover
        print(f"File {headinfo.psfimage} not found ")
        testflag = True
    if headinfo.psfimageflag:  # pragma: no cover
        if not (headinfo.convxflag):  # pragma: no cover
            print("Warning: x-size of convolution is smaller than psf image x-axis")
            testflag = True
        if not (headinfo.convyflag):  # pragma: no cover
            print("Warning: y-size of convolution is smaller than psf image y-axis")
            testflag = True
    if not (headinfo.maskimageflag):  # pragma: no cover
        print(f"File {headinfo.maskimage} not found ")
        testflag = True
    if not (headinfo.constraintsflag):  # pragma: no cover
        print(f"File {headinfo.constraints} not found ")
        testflag = True
    if not (headinfo.xsizeflag):  # pragma: no cover
        print("format of size in x-axis not valid. xmax < xmin")
        testflag = True
    if not (headinfo.ysizeflag):  # pragma: no cover
        print("Format of size in y-axis not valid. ymax < ymin")
        testflag = True

    if not testflag:  # pragma: no cover
        print("No issues with the header information")
    else:
        print(
            "There are some issues with the header information. Check messages above."
        )

    print("Total number of model components: ", len(galax))
    print("Total number of galaxies: ", len(np.unique(galax)))
    print("Components per galaxy: ")

    cnt = Counter(galax)
    Flux = 10 ** ((25 - mag) / 2.5)
    for item in np.unique(galax):
        maskgal = galax == item
        totFlux = Flux[maskgal].sum()
        totmag = -2.5 * np.log10(totFlux) + 25
        print(
            f"galaxy {int(item)} has {cnt[item]} components and a total mag of: {totmag:.2f} "
        )

    print(f"Total number of free parameters: {freepar} ")
    return 0


def mainGetBoxSizeDs9(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="Compute Box size from a DS9 region file for GALFIT header"
    )
    parser.add_argument("RegFile", help="DS9 region file containing the box region")
    args = parser.parse_args(argv)

    xmin, xmax, ymin, ymax = getBoxSizeDs9(args.RegFile)
    print(f"xmin, xmax, ymin, ymax: {xmin} {xmax} {ymin} {ymax} ")
    return 0


def maingetSersic(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="prints the Sérsic function from a DS9 ellipse region"
    )
    parser.add_argument("Image", help="image FITS file")
    parser.add_argument("RegFile", help="DS9 ellipse region file")
    parser.add_argument(
        "-zp", "--zeropoint", type=float, default=25, help="photometric zero point"
    )
    parser.add_argument(
        "-sk", "--sky", type=float, default=0, help="sky background to subtract"
    )
    parser.add_argument(
        "-bt",
        "--bulgetot",
        type=float,
        help="bulge-to-total ratio (prints two Sérsics: bulge and disk)",
    )
    parser.add_argument(
        "-c", "--center", action="store_true", help="take center from DS9 region file"
    )
    parser.add_argument(
        "-n",
        "--noprint",
        action="store_true",
        help="do not print Sérsic functions to stdout",
    )
    parser.add_argument("-m", "--mask", type=str, help="mask file")
    parser.add_argument(
        "-o", "--out", default="sersic.init", type=str, help="output GALFIT file"
    )
    parser.add_argument(
        "-b",
        "--bards9",
        type=str,
        help="DS9 ellipse region containing the bar (requires --bulgetot)",
    )
    args = parser.parse_args(argv)

    getSersic(
        args.Image,
        args.RegFile,
        args.center,
        args.mask,
        args.zeropoint,
        args.sky,
        args.noprint,
        args.bulgetot,
        args.bards9,
        args.out,
    )
    return 0


def mainSersic2Ferrer(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="converts a Bar Sérsic function (2 component) to a Ferrer function"
    )
    parser.add_argument("galfile", help="GALFIT input file")
    parser.add_argument(
        "-a", "--alpha", action="store_true", help="keep Ferrer alpha parameter as free"
    )
    parser.add_argument(
        "-b", "--beta", action="store_true", help="keep Ferrer beta parameter as free"
    )
    parser.add_argument(
        "-o", "--out", default="ferrer.init", type=str, help="output GALFIT file"
    )

    args = parser.parse_args(argv)

    Sersic2Ferrer(
        args.galfile,
        args.alpha,
        args.beta,
        args.out,
    )

    return 0


def main_imarith(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="imarith: arithmetic operations on images"
    )
    parser.add_argument("ImageFile", help="the input image")
    parser.add_argument(
        "-o", "--output", type=str, default="output.fits", help="the output image"
    )
    parser.add_argument(
        "-i",
        "--image2",
        type=str,
        help=(
            "second input image; if provided, operations use ImageFile and image2; "
            "otherwise they use the given constants"
        ),
    )
    parser.add_argument("-a", "--add", type=float, help="add constant to image pixels")
    parser.add_argument("-d", "--div", type=float, help="divide all pixels by constant")
    parser.add_argument(
        "-m", "--mul", type=float, help="multiply all pixels by constant"
    )
    parser.add_argument(
        "-s", "--sub", type=float, help="subtract constant from all pixels"
    )
    args = parser.parse_args(argv)

    imarith(
        args.ImageFile, args.output, args.image2, args.add, args.mul, args.div, args.sub
    )
    print("done")
    return 0


def mainMakePSF(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="Makes a PSF model of a star using MGE"
    )
    parser.add_argument("image", help="image containing the star to be modeled")
    parser.add_argument("GalfitFile", help="GALFIT file to obtain header options")
    parser.add_argument(
        "Ds9regFile", help="DS9 ellipse region file containing the star"
    )
    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="use DS9 center; else find (x,y) peak within DS9 ellipse",
    )
    parser.add_argument(
        "-o", "--out", type=str, default="psf.fits", help="PSF model image"
    )
    parser.add_argument("-sig", "--sigma", type=str, help="sigma image")
    parser.add_argument(
        "-t", "--twist", action="store_true", help="use twist option for MGE"
    )
    parser.add_argument(
        "-ng", "--numgauss", type=int, help="number of Gaussians used for GALFIT"
    )
    args = parser.parse_args(argv)

    makePSF(
        args.GalfitFile,
        args.image,
        args.Ds9regFile,
        args.center,
        args.out,
        args.sigma,
        args.twist,
        args.numgauss,
    )
    return 0
