#!/usr/bin/env python3

import os
import sys

from astropy.io import fits


def imarith(
    ImageFile: str,
    output: str,
    image2: str,
    add: float,
    mul: float,
    div: float,
    sub: float,
) -> None:
    """makes arithmetic operations on the image

    Paremeters
    ----------
    ImageFile: str, name of the image file
    output: str, name of the output file
    image2: str, name of a second image, optional
    add: float, If a value is provided, it will be added to the image.
                Otherwise, if image2 is provided, the result will be
                the sum of image and image2.

    mul: float, If a value is provided, it will be multiplied to the image.
                Otherwise, if image2 is provided, the result will be
                the multiplication of image and image2.

    div: float, If a value is provided, it will be divided to the image.
                Otherwise, if image2 is provided, the result will be
                the division of image and image2.

    sub: float, If a value is provided, it will be substracted to the image.
                Otherwise, if image2 is provided, the result will be
                the substraction of image and image2.


    Returns
    -------
        None

    """
    if not os.path.exists(ImageFile):
        print("{} image filename does not exist!".format(ImageFile))
        sys.exit()

    if image2:
        if not os.path.exists(image2):
            print("{} image filename does not exist!".format(image2))
            sys.exit()

        hdu2 = fits.open(image2)
        dataImage2 = hdu2[0].data

    # original file
    hdu = fits.open(ImageFile)
    dataImage = hdu[0].data

    if add:
        if image2:
            newdata = dataImage + dataImage2
        else:
            newdata = dataImage + add

    if mul:
        if image2:
            newdata = dataImage * dataImage2
        else:
            newdata = dataImage * mul

    if div:
        if image2:
            newdata = dataImage / dataImage2
        else:
            newdata = dataImage / div

    if sub:
        if image2:
            newdata = dataImage - dataImage2
        else:
            newdata = dataImage - sub

    hdu[0].data = newdata
    hdu.writeto(output, overwrite=True)

    hdu.close()

    if image2:
        hdu2.close()
