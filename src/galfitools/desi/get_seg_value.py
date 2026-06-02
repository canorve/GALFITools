#!/usr/bin/env python3

import argparse
from pathlib import Path

from astropy.io import fits


def get_pixel_value(segmentation_file, x=None, y=None, ext=0):
    """
    Return the value of a selected pixel in a segmentation image.

    By default, the selected pixel is the central pixel of the image.

    Parameters
    ----------
    segmentation_file : str
        Path to the SExtractor segmentation FITS image.
    x : int, optional
        Pixel x-coordinate using Python zero-based indexing.
        If None, the central x-coordinate is used.
    y : int, optional
        Pixel y-coordinate using Python zero-based indexing.
        If None, the central y-coordinate is used.
    ext : int, optional
        FITS extension containing the segmentation image. Default is 0.

    Returns
    -------
    value : int or float
        Value of the selected pixel.
    x : int
        Selected x pixel coordinate using Python indexing.
    y : int
        Selected y pixel coordinate using Python indexing.
    """
    with fits.open(segmentation_file) as hdul:
        data = hdul[ext].data

    if data is None:
        raise ValueError(f"No image data found in extension {ext}")

    ny, nx = data.shape

    if x is None:
        x = nx // 2

    if y is None:
        y = ny // 2

    if not (0 <= x < nx):
        raise ValueError(f"x={x} is outside the image range [0, {nx - 1}]")

    if not (0 <= y < ny):
        raise ValueError(f"y={y} is outside the image range [0, {ny - 1}]")

    value = data[y, x]

    return value, x, y


def maingetSegValue():
    parser = argparse.ArgumentParser(
        description=(
            "Read a SExtractor segmentation image and return the value "
            "of a selected pixel. By default, the central pixel is used."
        )
    )

    parser.add_argument(
        "segmentation_file",
        help="Input SExtractor segmentation FITS image.",
    )

    parser.add_argument(
        "--x",
        type=int,
        default=None,
        help=(
            "Pixel x-coordinate using Python zero-based indexing. "
            "Default: image center."
        ),
    )

    parser.add_argument(
        "--y",
        type=int,
        default=None,
        help=(
            "Pixel y-coordinate using Python zero-based indexing. "
            "Default: image center."
        ),
    )

    parser.add_argument(
        "--ext",
        type=int,
        default=0,
        help="FITS extension containing the segmentation image. Default: 0.",
    )

    args = parser.parse_args()

    segmentation_file = Path(args.segmentation_file)

    if not segmentation_file.exists():
        raise FileNotFoundError(f"File not found: {segmentation_file}")

    value, x, y = get_pixel_value(
        segmentation_file,
        x=args.x,
        y=args.y,
        ext=args.ext,
    )

    print(f"Selected pixel position: x={x}, y={y}")
    print(f"Selected segmentation value: {value}")
    print(f"Equivalent DS9 position: x={x + 1}, y={y + 1}")

    return value


if __name__ == "__main__":
    maingetSegValue()
