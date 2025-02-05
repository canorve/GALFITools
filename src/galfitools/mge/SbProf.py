#! /usr/bin/env python3

import os
import os.path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from galfitools.galin.std import GetAxis
from galfitools.galin.std import Ds9ell2Kronellv2
from galfitools.galin.std import GetInfoEllip
from galfitools.galin.std import GetPmax
from galfitools.galin.std import GetExpTime
from galfitools.mge.mge2galfit import GetSize
from matplotlib.ticker import (
    AutoLocator,
    AutoMinorLocator,
    LogLocator,
    NullFormatter,
)
from mgefit.sectors_photometry import sectors_photometry
from scipy import stats


def sbProf(args):
    """Creates a surface brightness profile

    Parameters
    ----------
    args : argument parser data class
           from library argparser

    The parameters of the args class are
    explained below:

    args.Image : str
                the FITS image
    args.Ds9Region : str
                    DS9 ellipse region file, which
                    must enclose the galaxy to be fitted.
    args.mgzpt : float
                magnitude zeropoint
    args.mask : str
                name of the mask image
    args.axrat : float
                axis ratio
    args.angle : float
                position angle
    args.sky : float
                value of the sky background
    args.plate : float
                plate scale of the image
    args.output : str
                name of the output plot
    args.center: bool
                If True, it uses the center of the DS9 ellipse.
                otherwise will use the pixel position of the
                peak
    args.ranx : list
                range of the x-axis in the plot.
    args.rany : list
                range of the y-axis in the plot.
    args.logx : bool
                if True  x axis is logarithmic
    args.pix : bool
                If true top axis indicates pixels
    args.grid : bool
                if True enables grid in plot

    args.rad : float
            if it is different than None it plots and vertical
            line at the indicated value

    args.rad2 : float
            A second vertical line. If it is different than
            None it plots and vertical line at the indicated value


    """

    image = args.Image
    ds9reg = args.Ds9Region
    mgzpt = args.mgzpt
    mask = args.mask
    axrat = args.axrat
    angle = args.angle
    sky = args.sky
    plate = args.plate
    output = args.output
    center = args.center
    ranx = args.ranx
    rany = args.rany
    logx = args.logx
    pix = args.pix
    grid = args.grid

    rad = args.rad
    rad2 = args.rad2

    # configuration

    conf = Config()

    conf.image = image
    conf.ds9reg = ds9reg
    conf.mgzpt = mgzpt
    conf.mask = mask
    conf.skylevel = sky
    conf.scale = plate
    conf.center = center
    conf.output = output
    conf.ranx = ranx
    conf.rany = rany
    conf.pix = pix
    conf.grid = grid
    conf.logx = logx
    conf.rad = rad
    conf.rad2 = rad2

    dataimg = readDataImg(conf)

    conf.exptime = GetExpTime(image)

    (ncol, nrow) = GetAxis(image)

    obj, xpos, ypos, rx, ry, angle2 = GetInfoEllip(conf.ds9reg)
    xx, yy, Rkron, theta, eps = Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle2)

    if center:
        print("center of ds9 ellipse region will be used")
        xpeak, ypeak = xpos, ypos
    else:
        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta + 90, eps, ncol, nrow)
        xpeak, ypeak = GetPmax(dataimg.img, dataimg.mask, xmin, xmax, ymin, ymax)

    # I have to switch x and y coordinates, don't ask me why

    # xpeak, ypeak = ypeak, xpeak

    conf.xc = xpeak
    conf.yc = ypeak

    if angle:
        conf.parg = angle
    else:
        conf.parg = theta

    if axrat:
        conf.qarg = axrat
        conf.eps = 1 - axrat
    else:
        conf.eps = eps
        conf.qarg = 1 - eps

    print("galaxy found at ", xpeak + 1, ypeak + 1)
    print("Axis ratio {:.2f}, Angle = {:.2f}".format(conf.qarg, conf.parg))

    print("Sky = ", sky)

    # removing background from galaxy and model images
    dataimg.img = dataimg.img - conf.skylevel
    # dataimg.model = dataimg.model - conf.skylevel

    #   numsectors=19 default
    numsectors = conf.sectors

    minlevel = conf.minlevel  # minimun value for sky

    # call to sectors_photometry for galaxy and model

    sectgalax = SectPhot(conf, dataimg, n_sectors=numsectors, minlevel=minlevel)

    print("creating plots..")

    limx, limy = EllipSectors(conf, sectgalax, n_sectors=numsectors)

    print("plot file: ", conf.output)

    ##############################################
    ##############################################
    ##############################################

    if conf.dplot:
        plt.pause(1.5)
    plt.savefig(conf.output, dpi=conf.dpival)


class DataImg:
    """Data class to save the galaxy image
    and mask image


    """

    img = np.empty([2, 2])

    mask = np.empty([2, 2])


def readDataImg(conf):
    """reads the galaxy image and
    mask image

    Parameters
    ----------
    conf :  Data class Config defined below


    Returns
    -------
    dataimg : Data class DataImg defined above which
            contains the images.


    """

    dataimg = DataImg()

    # reading galaxy and model images from file
    errmsg = "file {} does not exist".format(conf.image)

    assert os.path.isfile(conf.image), errmsg

    # hdu 1 => image   hdu 2 => model
    hdu = fits.open(conf.image)
    dataimg.img = (hdu[0].data.copy()).astype(float)
    hdu.close()

    # reading mask image from file

    if conf.mask:

        errmsg = "file {} does not exist".format(conf.mask)
        assert os.path.isfile(conf.mask), errmsg

        hdu = fits.open(conf.mask)
        mask = hdu[0].data
        dataimg.mask = np.array(mask, dtype=bool)
        hdu.close()

    else:
        dataimg.mask = np.array([])

    return dataimg


def SectPhot(conf, dataimg, n_sectors=19, minlevel=0):
    """Calls to function sectors_photometry
    for galaxy and model

    """

    maskb = dataimg.mask

    eps = 1 - conf.qarg

    # plt.clf()
    print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp = conf.xc
    xctemp = conf.yc

    angsec = conf.parg

    ###################################################
    # if fit == 'gal':
    #  galaxy:

    if conf.mask:
        sectimg = sectors_photometry(
            dataimg.img,
            eps,
            angsec,
            xctemp,
            yctemp,
            minlevel=minlevel,
            plot=conf.dplot,
            badpixels=maskb,
            n_sectors=n_sectors,
        )
    else:
        sectimg = sectors_photometry(
            dataimg.img,
            eps,
            angsec,
            xctemp,
            yctemp,
            minlevel=minlevel,
            plot=conf.dplot,
            n_sectors=n_sectors,
        )

    # if ellconf.dplot:
    plt.pause(1)  # Allow plot to appear on the screen
    plt.savefig(conf.namesec)

    return sectimg


def EllipSectors(conf, sectgalax, n_sectors=19, minlevel=0):
    """Creates the plot from data obtained by sectors_photometry"""

    # galaxy
    xradq, ysbq, ysberrq = sect2xy(sectgalax, conf, n_sectors)

    # Plotting
    limx, limy, axsec = PlotSB(xradq, ysbq, ysberrq, conf, conf.scale)

    axsec.legend(loc=1)

    return limx, limy


def sect2xy(sect, conf, n_sectors):
    """converts the sectors data to
    surface brightness data for plotting
    """

    #######################################
    #######################################
    # model
    stidx = np.argsort(sect.radius)

    mgerad = sect.radius[stidx]
    mgecount = sect.counts[stidx]
    mgeangle = sect.angle[stidx]
    mgeanrad = np.deg2rad(mgeangle)

    ab = conf.qarg

    aellab = mgerad * np.sqrt((np.sin(mgeanrad) ** 2) / ab**2 + np.cos(mgeanrad) ** 2)

    aellarc = aellab * conf.scale

    # formula according to cappellary mge manual
    mgesb = (
        conf.mgzpt
        - 2.5 * np.log10(mgecount / conf.exptime)
        + 2.5 * np.log10(conf.scale**2)
        + 0.1
    )  # - ellconf.Aext

    stidxq = np.argsort(aellarc)

    xarc = aellarc[stidxq]
    ymge = mgesb[stidxq]

    #  Function to order SB along X-axis for model

    xrad, ysb, ysberr = FindSB(xarc, ymge, n_sectors)

    return xrad, ysb, ysberr


def FindSB(xarcq, ymgeq, numsectors):
    """obtains the surface brightness
    from the divided sectors by angle

    # the xarcq array must be ordered
    # use mag instead of counts
    """
    xradq = []
    ysbq = []
    ysberrq = []
    xradq = np.array(xradq)
    ysbq = np.array(ysbq)
    ysberrq = np.array(ysberrq)

    numsave = 0
    tot = xarcq.size
    count = 0
    for i in range(tot, 0, -1):

        lima = i - numsectors
        limb = i

        if xarcq[lima:limb].size == 0:
            break
        else:
            valstd = np.std(xarcq[lima:limb])
            if valstd < 0.1:
                numsave = count
                break
            count = count + 1

    init = numsave % numsectors
    n = init

    num = np.int32((xarcq.size - init) / numsectors)
    n = xarcq.size - init
    for i in range(num, 0, -1):

        lima = n - numsectors
        limb = n

        xradq = np.append(xradq, np.mean(xarcq[lima:limb]))
        ysbq = np.append(ysbq, np.mean(ymgeq[lima:limb]))
        # ysberrq=np.append(ysberrq,np.std(ymgeq[lima:limb]))
        ysberrq = np.append(
            ysberrq, stats.sem(ymgeq[lima:limb])
        )  # standard error is more appropiate

        n = n - numsectors

    return xradq, ysbq, ysberrq


def PlotSB(xradq, ysbq, ysberrq, conf, scale):
    """Produces final best-fitting plot"""

    # subplot for arc sec axis
    plt.close("all")

    # ULISES begin
    # fig, (axsec,axred) = plt.subplots(2, sharex=True, sharey=False)
    # fig, axsec = plt.subplots(1)
    # gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
    # gs.update(hspace=0.07)
    # ULISES end

    fig, axsec = plt.subplots()

    if conf.ranx:
        (xmin, xmax) = conf.ranx[0], conf.ranx[1]

    if conf.rany:
        (ymin, ymax) = conf.rany[0], conf.rany[1]

    minrad = np.min(xradq)
    maxrad = np.max(xradq)

    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad / minrad) ** np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt / mincnt) ** np.array([+1.05, -0.05])  # inverted axis

    # ULISES begin
    # axsec = plt.subplot(gs[0])

    # axsec.set_ylabel(r"Surface Brightness $(mag\; arcsec^{-2})$")
    # ULISES end

    axsec.set_ylabel(r"Surface Brightness $(mag\; arcsec^{-2})$")

    axsec.errorbar(
        xradq,
        ysbq,
        yerr=ysberrq,
        fmt="o-",
        capsize=2,
        color="blue",
        markersize=0.7,
        label="galaxy",
        linewidth=2,
    )

    if conf.rany is True:
        axsec.set_ylim(ymax, ymin)  # inverted
    else:
        axsec.set_ylim(yran)

    if conf.logx is True:

        axsec.set_xscale("log")

        locmaj = LogLocator(base=10, numticks=12)
        axsec.xaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
        axsec.xaxis.set_minor_locator(locmin)
        axsec.xaxis.set_minor_formatter(NullFormatter())

    else:
        axsec.xaxis.set_minor_locator(AutoMinorLocator())
        axsec.xaxis.set_major_locator(AutoLocator())

    axsec.tick_params(which="both", width=2)
    axsec.tick_params(which="major", length=7)
    axsec.tick_params(which="minor", length=4, color="r")

    # begin psf fwhm
    # if ellconf.flagfwhm:
    #    xpos = ellconf.fwhm*scale
    #    axsec.axvline(x=xpos,  linestyle='--', color='k', linewidth=2)
    # end

    # ULISES begin
    # axsec.axes.xaxis.set_ticklabels([])
    # ULISES end

    # axsec.axes.xaxis.set_ticklabels([])

    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())

    # ULISES begin
    # Residual plot
    # if len(ysbq) < len(ysbm):
    #    ysbm = ysbm[len(ysbm)-len(ysbq):]

    # elif len(ysbq) > len(ysbm):
    #    ysbq = ysbq[len(ysbq)-len(ysbm):]
    #    ysberrq = ysberrq[len(ysberrq)-len(ysbq):]

    # if len(ysbq) < len(ysberrm):
    #    ysberrm = ysberrm[len(ysberrm)-len(ysbq):]
    # elif len(ysbq) > len(ysberrm):
    #    ysbq = ysbq[len(ysbq)-len(ysberrm):]

    # residual = ((ysbq-ysbm)/ysbq)*100 # (data-model)/data in percentage

    # err = ((ysbm/ysbq**2)**2) * ysberrq**2 + ((1/ysbq)**2) * ysberrm**2
    # err = np.sqrt(err)*100
    # axred = plt.subplot(gs[1])

    # if len(xradq) != len(residual):
    #    axred.errorbar(xradm, residual, yerr=err,fmt='.',capsize=2,color='k')
    # else:
    #    axred.errorbar(xradq, residual, yerr=err,fmt='.',capsize=2,color='k')

    # axred.axhline(y=0,ls='dashed', color='k')
    axsec.set_xlabel("Radius $(arcsec)$")
    # axred.set_ylabel('Residual (%)')
    # axred.set_ylim(-2,2)
    # ULISES end

    # if conf.flaglogx == True:

    #    axred.set_xscale("log")

    #    locmaj = LogLocator(base=10,numticks=12)
    #    axred.xaxis.set_major_locator(locmaj)

    #    locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
    #    axred.xaxis.set_minor_locator(locmin)
    #    axred.xaxis.set_minor_formatter(NullFormatter())
    # else:
    #    axred.xaxis.set_minor_locator(AutoMinorLocator())
    #    axred.xaxis.set_major_locator(AutoLocator())

    # axred.tick_params(which='both', width=2)
    # axred.tick_params(which='major', length=7)
    # axred.tick_params(which='minor', length=4, color='r')

    if conf.ranx:
        axsec.set_xlim(xmin, xmax)
    #  axred.set_xlim(xmin,xmax) #ulises plot
    else:
        axsec.set_xlim(xran)
    # axred.set_xlim(xran) #ulises plot

    if conf.pix is True:

        axpix = axsec.twiny()

        axpix.set_xlabel("Pixels")
        x1, x2 = axsec.get_xlim()

        axpix.set_xlim(x1 / scale, x2 / scale)

        axpix.figure.canvas.draw()

        if conf.logx is True:
            axpix.set_xscale("log")
            locmaj = LogLocator(base=10, numticks=12)
            axpix.xaxis.set_major_locator(locmaj)

            locmin = LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8), numticks=12)
            axpix.xaxis.set_minor_locator(locmin)
            axpix.xaxis.set_minor_formatter(NullFormatter())
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())
            axpix.xaxis.set_major_locator(AutoLocator())

        axpix.tick_params(which="both", width=2)
        axpix.tick_params(which="major", length=7)
        axpix.tick_params(which="minor", length=4, color="r")

        axret = axsec
    else:
        axret = axsec

    if conf.grid is True:
        # Customize the major grid
        axsec.grid(which="major", linestyle="-", linewidth="0.7", color="black")
        # Customize the minor grid
        axsec.grid(which="minor", linestyle=":", linewidth="0.5", color="black")

    # change the linewidth of the axis
    for axis in ["top", "bottom", "left", "right"]:
        axsec.spines[axis].set_linewidth(1.5)
        # axred.spines[axis].set_linewidth(1.5)

    if conf.rad:

        axsec.axvline(x=conf.rad, linestyle="--", color="m", linewidth=2)

    if conf.rad2:

        axsec.axvline(x=conf.rad2, linestyle="--", color="c", linewidth=2)

    return xran, yran, axret


class Config:
    """Data class for configuration parameters
    of SbProf
    """

    logx = False
    image = "none.fits"

    mgzpt = 25
    ds9reg = "none.reg"
    mask = "mask.fits"
    scale = 1
    center = False
    exptime = 1

    namesec = "gal.png"

    # init
    xc = 1
    yc = 1

    ranx = [0, 0]
    rany = [0, 0]

    qarg = 1
    parg = 0

    dpival = 150

    minlevel = 0
    sectors = 19

    image = "none.fits"

    # output file
    output = "out.png"

    pix = False
    grid = False

    # sky parameters:
    skylevel = 0

    dplot = True

    Aext = 0  # surface brightness correction for plots

    rad = 0
    rad2 = 0
