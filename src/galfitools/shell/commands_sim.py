import argparse

from galfitools.shell.prt import printWelcome
from galfitools.sim.MakeSim import makeSim


def _build_parser_make_sim() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="simulates an observed galaxy from a GALFIT model"
    )
    p.add_argument("image", help="the GALFIT galaxy model")
    p.add_argument("newimage", help="the name of the new galaxy image")
    p.add_argument(
        "-s",
        "--sky",
        type=float,
        help="the sky background value. default = 0",
        default=0,
    )
    p.add_argument(
        "-std",
        "--std",
        type=float,
        help="the sky standard deviation. default = 1",
        default=1,
    )
    p.add_argument(
        "-g",
        "--gain",
        type=float,
        help="the gain value of the image. default = 1",
        default=1,
    )
    return p


def mainMakeSim(argv=None) -> int:
    """
    Parse args and call makeSim. Return 0 on success.
    Accepts an optional argv list for testing.
    """
    printWelcome()

    parser = _build_parser_make_sim()
    args = parser.parse_args(argv)

    image = args.image
    GAIN = args.gain

    skymean = args.sky
    skystd = args.std

    newimage = args.newimage

    makeSim(image, GAIN, skymean, skystd, newimage)

    return 0
