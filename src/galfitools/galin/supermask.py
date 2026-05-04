#! /usr/bin/env python3


import argparse

from pathlib import Path

from galfitools.galin.MaskDs9 import maskDs9
from galfitools.desi.get_seg_value import get_pixel_value
from galfitools.galin.imarith import imarith


def superMask(
    segmentation_file: str,
    masksky_file: str,
    ds9ellipse_file: str,
    output="supermask.fits",
    rem_masksky=False,
):
    """
    Combines SExtractor segmentation fits file, mask fits file (created with maskSky),
    and DS9 region file. It removes the central galaxy.


    Parameters
    ----------
    segmentation_file : str
        Path to the SExtractor segmentation FITS image.
    masksky_file: str
        Path to the mask FITS image.
    ds9ellipse_file: str
        Path to the DS9 region file.
    rem_masksky: bool
       If True, removes central galaxy from masksky using DS9 ellipse maskbits region file
       recommended if central galaxy has not been removed from this mask.
    output: str optional
        name of the output mask fits file. Default: megamask.fits

    Returns
    -------
    True

    Note
    -----
    DS9 region file can be a box, ellipse or polygon. The ones
    that are accepted by maskDs9 function

    """

    # reads number of the central galaxy from SExtractor segmentation file
    pix_value, x, y = get_pixel_value(segmentation_file)

    # removes central galaxy from segmentation_file
    # and optional for the masksky image. It is expected that
    # galaxy is already removed for masksky

    maskDs9(
        segmentation_file,
        ds9ellipse_file,
        0,
        pixval=pix_value,
    )

    if rem_masksky:
        maskDs9(
            masksky_file,
            ds9ellipse_file,
            0,
        )

    # FINAL MASk: segmentation file + masksky = supermask
    imarith(segmentation_file, output, masksky_file, add=0)

    return True


def mainsuperMask():
    parser = argparse.ArgumentParser(
        description=(
            "Read SExtractor segmentation image, mask created by maskSky and DS9 region file "
            "to create a super mask. It removes the central galaxy"
        )
    )
    parser.add_argument(
        "segmentation_file",
        help="Input SExtractor segmentation FITS image.",
    )

    parser.add_argument(
        "masksky_file",
        help="Input mask created with maskSky (or any other mask) file FITS image.",
    )

    parser.add_argument(
        "ds9ellipse_file",
        help="Input DS9 ellipse region file.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="supermask.fits",
        help=("output FITS file. " "Default: supermask"),
    )

    parser.add_argument(
        "--rem_masksky",
        action="store_true",
        help="remove central galaxy from masksky.",
    )

    args = parser.parse_args()

    segmentation_file = Path(args.segmentation_file)

    if not segmentation_file.exists():
        raise FileNotFoundError(f"File not found: {segmentation_file}")

    masksky_file = Path(args.masksky_file)

    if not masksky_file.exists():
        raise FileNotFoundError(f"File not found: {masksky_file}")

    ds9ellipse_file = Path(args.ds9ellipse_file)

    if not ds9ellipse_file.exists():
        raise FileNotFoundError(f"File not found: {maskbits_file}")

    superMask(
        segmentation_file,
        masksky_file,
        ds9ellipse_file,
        args.output,
        args.rem_masksky,
    )

    print(f"super mask created: {args.output}")


if __name__ == "__main__":
    mainsuperMask()
