#!/usr/bin/env python3

import os.path
import sys

import numpy as np


def getBoxSizeDs9(RegFile: str) -> int:

    if not os.path.exists(RegFile):
        print("%s: reg filename does not exist!" % (sys.argv[1]))
        sys.exit(1)

    v0 = []
    v1 = []
    v2 = []
    v3 = []
    v4 = []
    v5 = []

    f1 = open(RegFile, "r")

    lines = f1.readlines()

    f1.close()

    flag = False
    # reading reg file
    for line in lines:

        line = line.split("#")
        line = line[0]

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "box":

            x0 = "box"
            x2 = x1[4:]
            flag = True

        if flag is True:

            x3 = p[4]
            x4 = x3[:-2]

            v0.append(x0)
            v1.append(float(x2) - 1)
            v2.append(float(p[1]) - 1)
            v3.append(float(p[2]))
            v4.append(float(p[3]))
            v5.append(float(x4))

            flag = False

    xpos = np.array(v1)
    ypos = np.array(v2)
    rx = np.array(v3)
    ry = np.array(v4)

    # converting for galfit header:

    xmin = xpos - rx / 2
    xmax = xpos + rx / 2

    ymin = ypos - ry / 2
    ymax = ypos + ry / 2

    # [0] if there are many, only the first one is take in account
    xmin = round(xmin[0])
    xmax = round(xmax[0])

    ymin = round(ymin[0])
    ymax = round(ymax[0])

    return xmin, xmax, ymin, ymax


#############################################################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################
