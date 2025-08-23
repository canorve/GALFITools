import argparse

from galfitools.mge.mge2galfit import mge2gal
from galfitools.mge.SbProf import sbProf
from galfitools.shell.prt import printWelcome


# ---- helpers (small, testable) ----
def _build_parser_mge() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="fits mge of Cappellari and formats to GALFIT"
    )
    p.add_argument("GalfitFile", help="GALFIT file to obtain the header options")
    p.add_argument(
        "Ds9regFile",
        help=(
            "the DS9 ellipse region file containing the galaxy. "
            "Ignored if -xy option is used"
        ),
    )
    p.add_argument(
        "-t", "--twist", action="store_true", help="uses twist option for mge "
    )
    p.add_argument(
        "-c",
        "--center",
        action="store_true",
        help=(
            "uses the center given in DS9 region file, "
            "otherwise it will found the x,y peaks within DS9 ellipse"
        ),
    )
    p.add_argument("-p", "--psf", type=float, help="the value of PSF sigma ", default=0)
    p.add_argument(
        "-gas",
        "--gauss",
        action="store_true",
        help="uses gauss function for galfit file",
    )
    p.add_argument(
        "-fser",
        "--freeser",
        action="store_true",
        help="leaves the sersic index as a free parameter to fit",
    )
    p.add_argument(
        "-fsk",
        "--freesky",
        action="store_true",
        help="leaves the sky as a free parameter to fit",
    )
    p.add_argument(
        "-ng",
        "--numgauss",
        type=int,
        help="number of gaussians that will be used for galfit.",
    )
    p.add_argument(
        "-xy",
        "--xypos",
        nargs=2,
        type=int,
        help="provides the (x y) position center of the object to fit",
    )
    p.add_argument("-e", "--ellip", type=float, help="ellipticity of object")
    p.add_argument(
        "-pa",
        "--posang",
        type=float,
        help="position angle of object. Measured from Y-axis",
    )
    return p


# check modify: remove initgauss
def mainMGE():
    """
    Calls the mge2gal function based on argument parsing.
    This function serves as an example of an API.
    """
    printWelcome()
    parser = _build_parser_mge()

    args = parser.parse_args(argv)

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
    return 0


def _build_parser_sbprof() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SbProf: creates a surface brightness profile from a ellipse ds9 region"
    )
    p.add_argument("Image", help="image fits file")
    p.add_argument("Ds9Region", help="Ds9 ellipse region file")
    p.add_argument("-q", "--axrat", type=float, help="axis ratio")
    p.add_argument(
        "-pa", "--angle", type=float, help="angular position (same as GALFIT)"
    )
    p.add_argument(
        "-mz", "--mgzpt", type=float, help="Magnitude zero point", default=25
    )
    p.add_argument("-m", "--mask", type=str, help="mask fits file")
    p.add_argument("-s", "--sky", type=float, help="sky value. Default = 0", default=0)
    p.add_argument("-p", "--plate", type=float, help="plate scale ", default=1)
    p.add_argument("-o", "--output", type=str, help="output file", default="sb.png")
    p.add_argument(
        "-c",
        "--center",
        action="store_true",
        help=(
            "uses the center given in DS9 region file, "
            "otherwise it will found the x,y peaks within DS9 ellipse"
        ),
    )
    p.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="provide a range for x-axis: xmin - xmax ",
    )
    p.add_argument(
        "-ry",
        "--rany",
        nargs=2,
        type=float,
        help="provide a range for y-axis: ymin - ymax  ",
    )
    p.add_argument(
        "-lx", "--logx", action="store_true", help="turn the X-axis to logarithm "
    )
    p.add_argument(
        "-px", "--pix", action="store_true", help="turn the top x-axis in pixels "
    )
    p.add_argument(
        "-g", "--grid", action="store_true", help="display a grid in the plot "
    )
    p.add_argument(
        "-r", "--rad", type=float, help="value for a vertical line to add into the plot"
    )
    p.add_argument(
        "-r2",
        "--rad2",
        type=float,
        help="value for a second vertical line to add into the plot",
    )
    return p


def mainSbProf():
    """
    Calls the sbProf function based on argument parsing.
    This function serves as an example of an API.
    """

    printWelcome()

    parser = _build_parser_sbprof()
    args = parser.parse_args(argv)

    sbProf(args)

    print("Done")

    return 0


if __name__ == "__main__":
    # Optional: allow running this module directly
    sys.exit(mainMGE())  # or choose which to run
