#! /usr/bin/env python3
# commands_batch.py
import argparse

from galfitools.galout.getBarSize import getMulBarSize


def mainGetMulBarSize(argv=None) -> int:
    printWelcome()
    p = argparse.ArgumentParser(
        description="getBarSize: gets the bar size from Sersic and Ferrer models"
    )
    p.add_argument("InputFile", help="file containing a list of file path GALFIT files")
    p.add_argument(
        "-d", "--dis", type=int, default=3, help="maximum distance among components"
    )
    p.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=1,
        help="number of component to be selected",
    )
    p.add_argument(
        "-o",
        "--out",
        type=str,
        default="barlength.reg",
        help="output DS9 ellipse region",
    )
    p.add_argument(
        "-co", "--output", type=str, default="barlenghts.csv", help="output csv file"
    )
    p.add_argument(
        "-p",
        "--plot",
        action="store_true",
        help="plots a file of kappa and break radius",
    )
    p.add_argument(
        "-r",
        "--red",
        action="store_true",
        help="If activated, DS9 region ellipse is red",
    )
    p.add_argument(
        "-rx",
        "--ranx",
        nargs=2,
        type=float,
        help="range of radius to search for barlength",
    )
    a = p.parse_args(argv)

    ret = getMulBarSize(
        a.InputFile, a.dis, a.numcomp, a.plot, a.ranx, a.out, a.red, a.output
    )
    print("done ")

    return 0
