import argparse

from galfitools.shell.prt import printWelcome
from galfitools.sky.GalfitSky import galfitSky

# from galfitools.sky.Sky import sky
from galfitools.sky.SkyDs9 import SkyDs9
from galfitools.sky.SkyRing import SkyRing


def mainGalfitSky():
    """
    Calls the galfitSky function based on argument parsing.
    This function serves as an example of an API.
    """

    printWelcome()
    parser = argparse.ArgumentParser(description="computes the sky using GALFIT")

    parser.add_argument("image", help="the image file")
    parser.add_argument("mask", help="the GALFIT mask file")

    parser.add_argument(
        "-s", "--scale", type=float, help="the plate scale. default = 1", default=1
    )

    parser.add_argument(
        "-zp",
        "--mgzpt",
        type=float,
        help="the magnitud zero point. default=25",
        default=25,
    )

    parser.add_argument(
        "-x", "--xpos", type=float, help="the x position. default=1", default=1
    )
    parser.add_argument(
        "-y", "--ypos", type=float, help="the y position. default=1", default=1
    )

    parser.add_argument(
        "-is",
        "--initsky",
        type=float,
        help="the initial sky value default=0",
        default=0,
    )

    args = parser.parse_args()

    imgname = args.image
    maskfile = args.mask

    mgzpt = args.mgzpt
    scale = args.scale

    X = args.xpos
    Y = args.ypos

    initsky = args.initsky

    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)


"""
# Deprecated
def mainSky():
    ''' 
    Calls the sky function based on argument parsing.
    This function serves as an example of an API.
    ''' 

    printWelcome()

    parser = argparse.ArgumentParser(
        description="computes sky from a ds9 region box file"
    )

    parser.add_argument("image", help="the image file")
    parser.add_argument("maskfile", help="the GALFIT Mask image file ")
    parser.add_argument("Ds9regFile", help="the DS9 box region file")

    args = parser.parse_args()

    imgname = args.image
    maskimage = args.maskfile
    filereg = args.Ds9regFile

    mean, sig = sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:")

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))
"""


def mainSkyDs9():
    """
    Calls the SkyDs9 unction based on argument parsing.
    This function serves as an example of an API.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="SkyDs9: computes sky background from a Ds9 region file:"
        + " Box, Ellipses and Polygons "
    )

    parser.add_argument(
        "ImageFile", help="the image file where the photometry will be computed"
    )
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument("-m", "--mask", type=str, help="the mask file")

    parser.add_argument(
        "-ol",
        "--outliers",
        action="store_true",
        help="Removes the top 80%% and botttom 20%% of the pixel values to "
        + "compute sky within ring",
    )

    args = parser.parse_args()

    ImageFile = args.ImageFile
    RegFile = args.RegFile
    maskfile = args.mask
    outliers = args.outliers

    mean, sig = SkyDs9(ImageFile, RegFile, maskfile, outliers)

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))


def mainSkyRing():
    """
    Calls the SkyRing function based on argument parsing.
    This function serves as an example of an API.
    """

    printWelcome()

    parser = argparse.ArgumentParser(
        description="SkyRing: computes the sky using gradient over rings"
    )

    # required arguments
    parser.add_argument("Image", help="Fits image of the objects")

    parser.add_argument("Ds9regFile", help="the DS9 ellipse region file")

    parser.add_argument("-m", "--mask", type=str, help="Fits Mask image")

    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="use the center of the ellipse. Otherwise it will use the (x,y)"
        + " position with the highest value of the ellipse",
    )

    parser.add_argument(
        "-ol",
        "--outliers",
        action="store_true",
        help="Removes the top 80%% and botttom 20%% of the pixel values to"
        + " compute sky within ring",
    )

    # arguments with inputs

    parser.add_argument(
        "-w",
        "--width",
        type=int,
        help="width of the ring for the gradient method. Default = 20. ",
        default=20,
    )

    args = parser.parse_args()

    image = args.Image
    mask = args.mask
    ds9regfile = args.Ds9regFile
    width = args.width
    center = args.center
    outliers = args.outliers

    print("Major axis of ellipse is used as initial radius.")

    mean, std, median, rad = SkyRing(image, mask, ds9regfile, width, center, outliers)

    line = (
        "Total sky:  mean = "
        + " {:.2f}; std={:.2f}; median = {:.2f} at radius {:.2f} ".format(
            mean, std, median, rad
        )
    )
    print(line)
