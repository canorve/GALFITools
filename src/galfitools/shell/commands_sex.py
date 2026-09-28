import argparse

from galfitools.sex.MakeMask import makeMask
from galfitools.shell.prt import printWelcome


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

    parser.add_argument(
        "-rd",
        "--region_dir",
        type=str,
        default="kron_regions",
        help="output ellipse region ID catalog. 'None' for all ellipses",
    )

    parser.add_argument(
        "-ri",
        "--region_id",
        type=int,
        default=None,
        help="output ellipse region ID catalog. Number for specific ID ellipse or 'None' for all ellipses",
    )

    args = parser.parse_args(argv)

    makeMask(
        args.Sexfile,
        args.ImageFile,
        args.maskout,
        args.scale,
        args.satds9,
        args.region_dir,
        args.region_id,
    )
    print("Done. Mask image created ")
    return 0
