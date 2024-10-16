import argparse

from galfitools.mge.mge2galfit import mge2gal
from galfitools.mge.SbProf import sbProf
from galfitools.shell.prt import printWelcome


# check modify: remove initgauss
def mainMGE():
    """
    Calls the mge2gal function based on argument parsing.
    This function serves as an example of an API.
    """
    printWelcome()
    parser = argparse.ArgumentParser(
        description="fits mge of Cappellari and formats to GALFIT"
    )

    # parser.add_argument("image", help="the Mask image file to modify or create")
    parser.add_argument("GalfitFile", help="GALFIT file to obtain the header options")
    parser.add_argument(
        "Ds9regFile",
        help="the DS9 ellipse region file containing the galaxy."
        + " Ignored if -xy option is used",
    )
    parser.add_argument(
        "-t", "--twist", action="store_true", help="uses twist option for mge "
    )
    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="uses the center given in DS9 region file, "
        + "otherwise it will found the x,y peaks within DS9 ellipse",
    )
    parser.add_argument(
        "-p", "--psf", type=float, help="the value of PSF sigma ", default=0
    )
    parser.add_argument(
        "-gas",
        "--gauss",
        action="store_true",
        help="uses gauss function for galfit file",
    )

    parser.add_argument(
        "-fser",
        "--freeser",
        action="store_true",
        help="leaves the sersic index as a free parameter to fit",
    )
    parser.add_argument(
        "-fsk",
        "--freesky",
        action="store_true",
        help="leaves the sky as a free parameter to fit",
    )

    parser.add_argument(
        "-ng",
        "--numgauss",
        type=int,
        help="number of gaussians that will be used for galfit. ",
    )

    parser.add_argument(
        "-xy",
        "--xypos",
        nargs=2,
        type=int,
        help="provides the (x y) position center of the object to fit",
    )

    parser.add_argument("-e", "--ellip", type=float, help="ellipticity of object")
    parser.add_argument(
        "-pa",
        "--posang",
        type=float,
        help="position angle of object. Measured from Y-axis",
    )

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    regfile = args.Ds9regFile
    twist = args.twist
    # regu = args.regu
    center = args.center
    psf = args.psf
    gauss = args.gauss

    # psfile = args.psfile
    # sigfile = args.sigfile

    freeser = args.freeser
    freesky = args.freesky

    numgauss = args.numgauss

    xypos = args.xypos
    ellip = args.ellip
    posang = args.posang

    # mge2gal(args)

    mge2gal(
        galfitFile,
        regfile,
        center,
        psf,
        twist,
        gauss,
        freeser,
        freesky,
        numgauss,
        xypos=xypos,
        ellip=ellip,
        posang=posang,
    )


def mainSbProf():
    """
    Calls the sbProf function based on argument parsing.
    This function serves as an example of an API.
    """

    printWelcome()
    parser = argparse.ArgumentParser(
        description="SbProf: creates a surface brightness profile"
        + " from a ellipse ds9 region"
    )

    # required arguments
    parser.add_argument("Image", help="image fits file")
    parser.add_argument("Ds9Region", help="Ds9 ellipse region file")

    parser.add_argument("-q", "--axrat", type=float, help="axis ratio")
    parser.add_argument(
        "-pa", "--angle", type=float, help="angular position (same as GALFIT)"
    )

    parser.add_argument(
        "-mz", "--mgzpt", type=float, help="Magnitud zero point", default=25
    )
    parser.add_argument("-m", "--mask", type=str, help="mask fits file")

    parser.add_argument(
        "-s", "--sky", type=float, help="sky value. Default = 0", default=0
    )
    parser.add_argument("-p", "--plate", type=float, help="plate scale ", default=1)
    parser.add_argument(
        "-o", "--output", type=str, help="output file", default="sb.png"
    )

    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="uses the center given in DS9 region file,"
        + "otherwise it will found the x,y peaks within DS9 ellipse",
    )

    parser.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for x-axis: xmin - xmax ",
    )
    parser.add_argument(
        "-ry",
        "--rany",
        nargs=2,
        type=float,
        help="provide a range for y-axis: ymin - ymax  ",
    )

    parser.add_argument(
        "-lx", "--logx", action="store_true", help="turn the X-axis to logarithm "
    )

    parser.add_argument(
        "-px", "--pix", action="store_true", help="turn the top x-axis in pixels "
    )
    parser.add_argument(
        "-g", "--grid", action="store_true", help="display a grid in the plot "
    )

    parser.add_argument(
        "-r", "--rad", type=float, help="value for a vertical line to add into the plot"
    )
    parser.add_argument(
        "-r2",
        "--rad2",
        type=float,
        help="value for a second vertical line to add into the plot",
    )

    args = parser.parse_args()

    sbProf(args)

    print("Done")
