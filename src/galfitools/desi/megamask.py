#! /usr/bin/env python3


import argparse

from pathlib import Path
import os

from galfitools.galin.MaskDs9 import maskDs9
from galfitools.desi.central_ellipse import central_ellipse
from galfitools.desi.get_seg_value import get_pixel_value
from galfitools.galin.imarith import imarith


def megaMask(
    segmentation_file: str,
    masksky_file: str,
    maskbits_file: str,
    output="megamask.fits",
    rem_masksky=False,
):
    """
    Combines SExtractor segmentation fits file, mask fits file (created with maskSky),
    and DESI maskbits fits file. It removes the central galaxy mask.


    Parameters
    ----------
    segmentation_file : str
        Path to the SExtractor segmentation FITS image.
    masksky_file: str
        Path to the mask FITS image.
    maskbits_file: str
        Path to the DESI maskbits FITS image.
    rem_masksky: bool
       If True, removes central galaxy from masksky using DS9 ellipse maskbits region file
       recommended if central galaxy has not been removed from this mask.
    output: str optional
        name of the output mask fits file. Default: megamask.fits

    Returns
    -------
    True



    """

    # it obtains the DS9  central galaxy ellipse
    ds9maskbits = "ellipse_maskbits.reg"
    central_ellipse(maskbits_file, ds9maskbits, refine=True, use_bit=True)

    maskbits_value = 4096  # DESI value for galaxies
    # reads number of the central galaxy from SExtractor segmentation file
    pix_value, x, y = get_pixel_value(segmentation_file)

    # removes central galaxy from maskbits and segmentation_file
    # and optional for the masksky image. It is expected that
    # galaxy is removed

    maskDs9(
        maskbits_file,
        ds9maskbits,
        0,
        pixval=maskbits_value,
    )

    maskDs9(
        segmentation_file,
        ds9maskbits,
        0,
        pixval=pix_value,
    )

    if rem_masksky:
        maskDs9(
            masksky_file,
            ds9maskbits,
            0,
        )

    # it proceeds to create the megamask

    tempmask = "tempmask.fits"

    # seg + masksky = tempmask
    imarith(segmentation_file, tempmask, masksky_file, add=0)

    # FINAL MASk: tempmask + maskbits = megamask
    imarith(tempmask, output, maskbits_file, add=0)

    # removing temporal mask
    tempmask_file = Path(tempmask)

    if tempmask_file.exists():
        os.remove(tempmask)

    return True


def mainmegaMask():
    parser = argparse.ArgumentParser(
        description=(
            "Read SExtractor segmentation image, mask created by maskSky and DESI maskbits "
            "to create a mega mask. It removes the central galaxy mask."
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
        "maskbits_file",
        help="Input DESI maskbits FITS image.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="megamask.fits",
        help=("output FITS file. " "Default: megamask"),
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

    maskbits_file = Path(args.maskbits_file)

    if not maskbits_file.exists():
        raise FileNotFoundError(f"File not found: {maskbits_file}")

    megaMask(
        segmentation_file,
        masksky_file,
        maskbits_file,
        args.output,
        args.rem_masksky,
    )

    print(f"mega mask created: {args.output}")


if __name__ == "__main__":
    maingetMegaMask()
