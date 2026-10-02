#! /usr/bin/env python

import os
import os.path
from pathlib import Path

import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import MakeImage


def findNeigh(
    sexfile: str,
    N: str,
    scale: float,
    offset: int,
) -> None:

    """Creates a mask file from a catalog of SExtractor

    It creates a mask file for GALFIT using information from a
    SExtractor catalog. It includes masking of saturated regions.

    Parameters
    ----------
    sexfile : str
            name of the Sextractor catalog
    N: int,
            Object ID number in the catalog to search for neighbors
    scale: float,
            Scale factor by which the ellipse will be enlarged or diminished.
    offset: int,
            Value to be added to the radius by which the ellipse will be enlarged or diminished.

    Returns
    -------
    None

    """

    print(f"finding neighbors for object {N} \n")

    neigh, indexes = FindNeighbors2(sexfile, N, scale, offset)

    print("neighbors: \n")
    print(neigh)


def CheckFlag(val, check, maxx=128):
    """Check if SExtractor flag contains the check flag

    This function is useful to check if
    a object sextractor flag contains saturated region

    Parameters
    ----------
    val: SExtractor flag of the object
    check: SExtractor flag value to check
    maxx: maximum value of the flag. Default = 128

    Returns
    -------
    bool, returns True if found

    # repeated

    """

    flag = False
    mod = 1

    while mod != 0:

        res = int(val / maxx)

        if maxx == check and res == 1:

            flag = True

        mod = val % maxx

        val = mod
        maxx = maxx / 2

    return flag


def CheckOverlap(xpos, ypos, R, theta, q, xpos2, ypos2, R2, theta2, q2):
    "Check the distance of two ellipses. returns True if they overlap"

    flag = False

    theta = theta + 90  # converting from GALFIT to Sextractor position
    theta2 = theta2 + 90

    bim = q * R
    bim2 = q2 * R2

    theta = theta * np.pi / 180  # Rads!!!
    theta2 = theta2 * np.pi / 180  # Rads!!!
    dx = xpos2 - xpos
    dy = ypos2 - ypos

    dx2 = xpos - xpos2
    dy2 = ypos - ypos2

    dist = np.sqrt(dx**2 + dy**2)

    landa = np.arctan2(dy, dx)

    if landa < 0:
        landa = landa + 2 * np.pi

    landa2 = np.arctan2(dy, dx)

    if landa2 < 0:
        landa2 = landa2 + 2 * np.pi

    landa = landa - theta
    landa2 = landa2 - theta2

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)
    angle2 = np.arctan2(np.sin(landa2) / bim2, np.cos(landa2) / R2)

    xell = (
        xpos + R * np.cos(angle) * np.cos(theta) - bim * np.sin(angle) * np.sin(theta)
    )
    yell = (
        ypos + R * np.cos(angle) * np.sin(theta) + bim * np.sin(angle) * np.cos(theta)
    )

    xell2 = (
        xpos2
        + R2 * np.cos(angle2) * np.cos(theta2)
        - bim2 * np.sin(angle2) * np.sin(theta2)
    )
    yell2 = (
        ypos2
        + R2 * np.cos(angle2) * np.sin(theta2)
        + bim2 * np.sin(angle2) * np.cos(theta2)
    )

    dell = np.sqrt((xell - xpos) ** 2 + (yell - ypos) ** 2)
    dell2 = np.sqrt((xell2 - xpos2) ** 2 + (yell2 - ypos2) ** 2)

    distell = dell + dell2

    if dist <= distell:
        flag = True

    return flag


def FindNeighbors2(catfile: str, n: int, KronScale=1, offset=0):

    #  This subroutine find neighbors for every galaxy
    #  note for myself: make a tree code of this in the future

    overflag = 0

    flagcheck = 0
    flagcheck2 = 0
    flagsat = 4

    (
        Num,
        Alpha,
        Delta,
        XPos,
        YPos,
        Mag,
        Kr,
        Fluxrad,
        Isoa,
        Ai,
        E,
        Theta,
        Bkgd,
        Idx,
        Flag,
    ) = np.genfromtxt(
        catfile, delimiter="", unpack=True
    )  # sorteado

    val = np.where(n == Num)[0][0]
    idx = val.item()

    Angle = Theta - 90
    AR = 1 - E  # agregado
    RKron = KronScale * Ai * Kr + offset

    flagcheck = CheckFlag(Flag[idx], flagsat)
    A = np.empty((0))
    indexs = np.empty((0))

    for jdx, vali in enumerate(Num):

        flagcheck2 = CheckFlag(Flag[jdx], flagsat)

        if flagcheck == False and flagcheck2 == False:

            if idx != jdx:

                overflag = CheckOverlap(
                    XPos[idx],
                    YPos[idx],
                    RKron[idx],
                    Angle[idx],
                    AR[idx],
                    XPos[jdx],
                    YPos[jdx],
                    RKron[jdx],
                    Angle[jdx],
                    AR[jdx],
                )

                if overflag == True:
                    A = np.append(A, Num[jdx])
                    indexs = np.append(indexs, jdx)

    A = A.astype(int)
    indexs = indexs.astype(int)
    return A, indexs
