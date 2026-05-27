#! /usr/bin/env python

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional
import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect

from galfitools.galin.galfit import (
    Galfit,
    SelectGal,
    SelectComp,
    conver2Sersic,
    numComps,
)
from galfitools.galout.getRads import getBreak2, getKappa2
from scipy.special import gamma, gammainc, gammaincinv

# ============================================================================
# Core domain logic
# ============================================================================


def getBarSize(
    galfitFile: str,
    dis: int,
    num_comp: int,
    plot: bool,
    ranx: list,
    out: str,
    red: bool,
    scale=1.0,
    method="break_kappa",
) -> tuple[float, int, float]:
    """gets the bar size of the spiral galaxies

    It takes the average of Kappa radius (maximum curvature) and
    Break radius (maximum of double derivative) to estimate the
    bar size of the three composed model of bulge, bar, and disk.

    It assumes the bar model is the second component of the
    GALFIT file. Bar model can be a Sersic or Ferrer function.
    The rest of components must be Sersic (or related) models.

    Parameters
    ----------
    galfitFile : str
    dis : int
    num_comp : int
              Number of component where it'll obtain center
              of all components. in other words it selects
              the galaxy that contains the bar if simultaneous
              fitting of galaxies was used.
    plot : bool
            If True, it draws plots of the break and kappa radius
    ranx : list
        range of search (xmin to xmax) for the kappa radius and break
        radius. If None, it will search in a range of r=1 to 2.5*Re
        of effetive radius of the bar model.
    out : str
         Name of the output file for the DS9 ellipse region marking
         the bar.
    red : bool
            If True, draws DS9 region ellipse as red color

    scale: float
            constant to multiply the bar length. Default =1

    method: str
            indicates which method is used to measure the
            bar length. Options include 'break_kappa', 'break'
            'kappa', 're', 'disk', 'all'. Default='break_kappa'

    Returns
    -------
    rbar : float
           bar size in pixels
    N : int
        number of components of the galaxy
    theta : float
        angular position of the galactic's bar

    See also
    --------
    getBreak2 : get the break radius
    getKappa2 : get the kappa radius


    """

    galfit = Galfit(galfitFile)

    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps = SelectGal(comps, dis, num_comp)

    maskgal = comps.Active == 1
    # it assumes bar is positioned as second galfit component
    # i.e. [1]
    theta = comps.PosAng[maskgal][1]

    AxRat = comps.AxRat[maskgal][1]
    X = comps.PosX[maskgal][1]
    Y = comps.PosY[maskgal][1]

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components to compute bar size")
        print("exiting..")
        sys.exit(1)

    #########################
    # computing the slope
    #########################
    if ranx:  # pragma: no cover
        (xmin, xmax) = ranx[0], ranx[1]
    else:
        Re = comps.Rad[maskgal][
            1
        ]  # it assumes bar is positioned as second galfit component

        # it assumes bar size is in this range. Hopefully it founds the solution there:
        xmin = 1
        xmax = 2.5 * Re

        ranx = [xmin, xmax]

    rbar = 0

    options = ["break_kappa", "break", "kappa", "disk", "re", "all"]
    if not (method in options):
        print("option not found. Setting to break_kappa")
        print("options available: break_kappa, break, kappa, disk, re, all")
        method = "break_kappa"

    print(f"method used: {method}")

    if method == "break_kappa":
        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)

        rbar = scale * ((rbreak + rkappa) / 2)

    if method == "break":

        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar = scale * rbreak

    if method == "kappa":

        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar = scale * rkappa

    if method == "re":
        rbar = scale * comps.Rad[maskgal][1]

    if method == "disk":
        rbar = scale * comps.Rad[maskgal][1]

        rbar = findDisk(galfitFile, dis, theta, num_comp, plot, ranx)

    if method == "all":

        # rbar0 = EffRad

        rkappa, N2, theta2 = getKappa2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar1 = rkappa

        rbreak, N, theta = getBreak2(galfitFile, dis, theta, num_comp, plot, ranx)
        rbar2 = rbreak

        rbar3 = comps.Rad[maskgal][1]

        # rbar = (rbar0+rbar1 +rbar2+rbar3)/4
        rbar = scale * (rbar1 + rbar2 + rbar3) / 3

    # now it creates the ellipse region file
    fout = open(out, "w")

    if red:
        color = "red"
    else:
        color = "blue"

    line = "# Region file format: DS9 version 4.1 \n"
    fout.write(line)
    linea = "global color=" + color + " dashlist=8 3 width=2 "
    lineb = 'font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '
    linec = "fixed=0 edit=1 move=1 delete=1 include=1 source=1\n"
    line = linea + lineb + linec
    fout.write(line)
    line = "physical\n"
    fout.write(line)

    rbarminor = rbar * AxRat

    elline = "ellipse({:.2f}, {:.2f}, {:.2f}, {:.2f} {:.2f}) \n".format(
        X, Y, rbarminor, rbar, theta
    )
    fout.write(elline)

    fout.close()

    return rbar, N, theta


def findDisk(galfitFile, dis, num_comp, angle, plot, ranx):
    """Determines the barlength using the disk method.


    Computes the barlength when the difference between the
    galaxy surface brightness and disk is minimal. It
    is computed along the angle direction of the bar.
    It is assumed the bar is second component.


    Parameters
    ----------
    galfitFile : str
            name of the GALFIT file
    dis: float
        maximum distance among components
    num_comp: int
            Number of component from which the center of all
            components will be determined.
    angle: float
           Position angle of the major axis of the galaxy. If None
           it will take the angle of the second component. Bar is
           assumed to be the second component.
    plot: bool
          if True, it will makes plot of the second derivative
    ranx: list
           Specify the range for the plot’s x-axis: xmin to xmax
           for plotting and searching. If set to None, the default
           range is [0.1, 100].

    Returns
    -------
    rbar: float
            radius of the barlength  in pixels


    """

    galfit = Galfit(galfitFile)
    head = galfit.ReadHead()
    galcomps = galfit.ReadComps()

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    # check this later for Ferrer
    galcomps = conver2Sersic(galcomps)

    comps = SelectGal(galcomps, dis, num_comp)
    comps2 = SelectComp(galcomps, 2)  # it assumes bar is 2 component

    # taking the second component position angle

    maskgal = comps.Active == 1

    if angle:  # pragma: no cover
        theta = angle
    else:
        theta = comps.PosAng[maskgal][1]  # bar angle

    N = numComps(comps, "all")

    if N == 0:  # pragma: no cover
        print("not enough number of components  to proceed")
        print("exiting..")
        sys.exit(1)

    if plot:  # pragma: no cover

        if ranx:
            (xmin, xmax) = ranx[0], ranx[1]
        else:
            xmin = 0.1
            xmax = 100

        R = np.arange(xmin, xmax, 0.1)

        Idiff = getDiff(head, comps, comps2, R, theta)

        plt.clf()
        plt.plot(R, Idiff)
        plt.grid(True)
        plt.minorticks_on()
        plt.xlabel("Rad")
        plt.ylabel("I1 - I2")
        plt.savefig("Idiff.png")

        I1, I2 = getIs(head, comps, comps2, R, theta)

        plt.clf()
        plt.plot(R, I1)
        plt.plot(R, I2)
        plt.xscale("log")
        plt.xlabel("Rad")
        plt.ylabel("I")
        plt.grid(True)
        plt.minorticks_on()
        plt.savefig("findisk.png")

    maskgal = comps.Active == 1  # using active components only

    a = 0.1
    b = comps.Rad[maskgal][-1] * 10  # hope it doesn't crash

    try:
        rbulge = bisect(getDiffx, a, b, args=(head, comps, comps2, theta))
    except Exception:  # pragma: no cover
        print("solution not found in given range")
        rbulge = 0

    return rbulge


def getDiff(head1, comps1, comps2, R, theta):  # pragma: no cover
    """Calculates the surface brightness difference
    between two models at various radii I1 - I2.

    """
    Idiff = np.array([])

    for r in R:

        Ir1 = GetIr().Ir(head1, comps1, r, theta)
        Ir2 = GetIr().Ir(head1, comps2, r, theta)

        Ird = Ir1 - Ir2
        Idiff = np.append(Idiff, Ird)

    return Idiff


def getIs(head1, comps1, comps2, R, theta):  # pragma: no cover
    """Calculates the surface brightness of two
    models at various radii.

    """

    I1 = np.array([])
    I2 = np.array([])

    for r in R:

        Ir1 = GetIr().Ir(head1, comps1, r, theta)
        Ir2 = GetIr().Ir(head1, comps2, r, theta)

        I1 = np.append(I1, Ir1)
        I2 = np.append(I2, Ir2)

    return I1, I2


def getDiffx(r, head1, comps1, comps2, theta):
    """Calculates the surface brightness difference
    between two models at a specified radii I1 - I2.

    """

    Ir1 = GetIr().Ir(head1, comps1, r, theta)
    Ir2 = GetIr().Ir(head1, comps2, r, theta)

    # tol = .01
    tol = 0.05

    Irdx = Ir1 - Ir2 - tol

    return Irdx


class GetIr:
    """
    Class called by getDiff, getDiffx, getIs
    to obtain the surface brightness from  a set
    of Sersic functions at different radii

    the main method is Ir

    Methods
    -------
    Ir : obtains the surface brightness at R

    Itotser : method called by Ir to obtain the
              total surface brightness of the
              set components at R
    Iser : method called by Itotser to obtain
           the surface brightnness
           per component at R


    """

    def Ir(self, head, comps, R, theta):

        # comps.Rad = comps.Rad*head.scale
        comps.Flux = 10 ** ((head.mgzpt - comps.Mag) / 2.5)

        k = gammaincinv(2 * comps.Exp, 0.5)

        denom1 = (2 * np.pi * comps.Rad**2) * (np.exp(k))
        denom2 = (comps.Exp) * (k ** (-2 * comps.Exp))
        # denom3 = (gamma(2*comps.Exp))*(comps.AxRat)
        denom3 = gamma(2 * comps.Exp)

        denom = denom1 * denom2 * denom3

        comps.Ie = comps.Flux / denom

        maskgal = comps.Active == 1

        Itotr = self.Itotser(
            R,
            comps.Ie[maskgal],
            comps.Rad[maskgal],
            comps.Exp[maskgal],
            comps.AxRat[maskgal],
            comps.PosAng[maskgal],
            theta,
        )

        return Itotr

    def Itotser(
        self, R: float, Ie: list, rad: list, n: list, q: list, pa: list, theta: float
    ) -> float:

        ItotR = self.Iser(R, Ie, rad, n, q, pa, theta)

        return ItotR.sum()

    def Iser(
        self, R: float, Ie: list, Re: list, n: list, q: list, pa: list, theta: float
    ) -> float:
        """sersic flux to a determined R"""

        k = gammaincinv(2 * n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta)

        Ir = Ie * np.exp(-k * ((Rcor / Re) ** (1 / n) - 1))

        return Ir
