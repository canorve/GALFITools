# commands_sky.py
import argparse

from galfitools.shell.prt import printWelcome
from galfitools.sky.GalfitSky import galfitSky
from galfitools.sky.SkyDs9 import SkyDs9
from galfitools.sky.SkyRing import SkyRing


# ---------- parsers (small, testable) ----------
def _build_parser_galfit_sky() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="computes the sky using GALFIT")
    p.add_argument("image", help="the image file")
    p.add_argument("mask", help="the GALFIT mask file")
    p.add_argument(
        "-s", "--scale", type=float, help="the plate scale. default = 1", default=1
    )
    p.add_argument(
        "-zp",
        "--mgzpt",
        type=float,
        help="the magnitude zero point. default=25",
        default=25,
    )
    p.add_argument(
        "-x", "--xpos", type=float, help="the x position. default=1", default=1
    )
    p.add_argument(
        "-y", "--ypos", type=float, help="the y position. default=1", default=1
    )
    p.add_argument(
        "-is",
        "--initsky",
        type=float,
        help="the initial sky value. default=0",
        default=0,
    )
    return p


def _build_parser_sky_ds9() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SkyDs9: computes sky background from a DS9 region file: Box, Ellipses, Polygons"
    )
    p.add_argument(
        "ImageFile", help="the image file where the photometry will be computed"
    )
    p.add_argument("RegFile", help="the DS9 region file")
    p.add_argument("-m", "--mask", type=str, help="the mask file")
    p.add_argument(
        "-ol",
        "--outliers",
        action="store_true",
        help="remove top 80% and bottom 20% of pixel values when computing sky",
    )
    return p


def _build_parser_sky_ring() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="SkyRing: computes the sky using gradient over rings"
    )
    p.add_argument("Image", help="FITS image of the object")
    p.add_argument("Ds9regFile", help="the DS9 ellipse region file")
    p.add_argument("-m", "--mask", type=str, help="FITS mask image")
    p.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="use the center of the ellipse; else use the (x,y) of max value",
    )
    p.add_argument(
        "-ol",
        "--outliers",
        action="store_true",
        help="remove top 80% and bottom 20% of pixel values within the ring",
    )
    p.add_argument(
        "-w",
        "--width",
        type=int,
        help="ring width for gradient method. default = 20",
        default=20,
    )
    return p


# ---------- entry points (now testable) ----------
def mainGalfitSky(argv=None) -> int:
    """Parse args and call galfitSky. Return 0 on success."""
    printWelcome()
    parser = _build_parser_galfit_sky()
    args = parser.parse_args(argv)
    galfitSky(
        args.image,
        args.mask,
        args.mgzpt,
        args.scale,
        args.xpos,
        args.ypos,
        args.initsky,
    )
    return 0


def mainSkyDs9(argv=None) -> int:
    """Parse args and call SkyDs9. Return 0 on success."""
    printWelcome()
    parser = _build_parser_sky_ds9()
    args = parser.parse_args(argv)
    mean, sig = SkyDs9(args.ImageFile, args.RegFile, args.mask, args.outliers)
    print(f"mean sky: {mean:.3f} ")
    print(f"std sky: {sig:.3f} ")
    return 0


def mainSkyRing(argv=None) -> int:
    """Parse args and call SkyRing. Return 0 on success."""
    printWelcome()
    parser = _build_parser_sky_ring()
    args = parser.parse_args(argv)

    print("Major axis of ellipse is used as initial radius.")
    mean, std, median, rad = SkyRing(
        args.Image, args.mask, args.Ds9regFile, args.width, args.center, args.outliers
    )
    print(
        f"Total sky:  mean =  {mean:.2f}; std={std:.2f}; median = {median:.2f} at radius {rad:.2f} "
    )
    return 0
