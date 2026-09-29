#!/usr/bin/env python3


import argparse
from pathlib import Path

from galfitools.sex.MakeMask import makeMask
from galfitools.sex.MakeReg import makeReg
from galfitools.sex.filter_sextractor import filter_catalog
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


def mainMakeReg(argv=None) -> int:
    printWelcome()
    parser = argparse.ArgumentParser(
        description="creates Ds9 ellipse regions from a SExtractor catalog"
    )
    parser.add_argument("Sexfile", help="SExtractor catalog file")

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

    makeReg(
        args.Sexfile,
        args.scale,
        args.region_dir,
        args.region_id,
    )
    print("Done. Mask image created ")
    return 0


def mainFilterSex(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", type=Path, help="Input ASCII or ASCII_HEAD catalog")
    parser.add_argument("output", type=Path, help="Output filtered catalog")
    parser.add_argument(
        "--flags", type=int, default=4, help="Keep FLAGS below this value (default: 4)"
    )
    parser.add_argument(
        "--mag",
        type=float,
        default=18.0,
        help="Keep MAG_BEST below this value (default: 18)",
    )
    parser.add_argument(
        "--class-star",
        type=float,
        default=0.6,
        help="Keep CLASS_STAR below this value (default: 0.6)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace existing output"
    )
    args = parser.parse_args(argv)
    try:
        total, kept = filter_catalog(
            args.input,
            args.output,
            args.flags,
            args.mag,
            args.class_star,
            args.overwrite,
        )
    except (OSError, ValueError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(f"Objects read: {total}; retained: {kept}; rejected: {total - kept}")
    print(f"Output: {args.output}")

    return 0
