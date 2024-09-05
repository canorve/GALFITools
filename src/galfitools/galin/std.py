#! /usr/bin/env python


import numpy as np
from galfitools.galin.galfit import GalComps

from astropy.io import fits
import os.path


def MakeImage(newfits, sizex, sizey):
    """create a new blank FITS Image

    Parameters
    ----------
    newfits : str
              name of the new image
    sizex : int
            number of columns of the new image
    sizey : int
            number of rows of the new image

    Returns
    -------
    bool

    """

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        sp.run(
            [runcmd],
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            universal_newlines=True,
        )

    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex), dtype=np.float64)
    hdu.writeto(newfits, overwrite=True)

    return True


def GetSize(x, y, R, theta, ell, ncol, nrow):
    """Get the (x,y) coordinates that encompass the ellipse

    Parameters
    ----------
    x : float, x-center of ellipse
    y : float, y-center of ellipse
    R : float, major axis of ellipse
    theta : float, angular position of ellipse
    ell : float, ellipticity
    ncol : number of columns of the image
    nrow : number of rows of the image


    Returns
    -------
    xmin, xmax, ymin, ymax : Minimum and maximum coordinates
    that encompass the ellipse.

    """

    # k Check
    q = 1 - ell
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

    # getting size

    constx = np.sqrt(
        (R**2) * (np.cos(theta)) ** 2 + (bim**2) * (np.sin(theta)) ** 2
    )
    consty = np.sqrt(
        (R**2) * (np.sin(theta)) ** 2 + (bim**2) * (np.cos(theta)) ** 2
    )

    xmin = x - constx
    xmax = x + constx
    ymin = y - consty
    ymax = y + consty

    mask = xmin < 1
    if mask.any():
        xmin = 1

    mask = xmax > ncol
    if mask.any():
        xmax = ncol - 1

    mask = ymin < 1
    if mask.any():
        ymin = 1

    mask = ymax > nrow
    if mask.any():
        ymax = nrow - 1

    return (round(xmin), round(xmax), round(ymin), round(ymax))


def Ds9ell2Kronell(xpos, ypos, rx, ry, angle):
    """Converts DS9 ellipse parameters to geometrical parameters

    Parameters
    ----------
    obj : str
          object name
    xpos : float, x-center
    ypos : float, y-center
    rx : float, major or minor axis
    ry : float, minor or major axis
    angle : float, angular position

    Returns
    -------
    xx : float, x center position
    yy : float, y center positon
    Rkron : float, major axis
    theta: float, angular position measured from Y-axis
    e: float, ellipticity


    """

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle
        xx = xpos
        yy = ypos
    else:
        q = rx / ry
        e = 1 - q
        Rkron = ry
        theta = angle + 90
        xx = xpos
        yy = ypos

    return xx, yy, Rkron, theta, e


def GetInfoEllip(regfile: str):
    """gets ellipse information from DS9 region files

    Parameters
    ----------
    regfile : str
                DS9 region file

    Returns
    -------
    obj : str
          object name
    xpos : float, x-center
    ypos : float, y-center
    rx : float, major or minor axis
    ry : float, minor or major axis
    angle : float, angular position


    returns 0, 0, 0, 0, 0, 0 if ellipse region was not found in file


    """

    if not os.path.exists(regfile):
        print("%s: reg filename does not exist!" % (regfile))
        sys.exit()

    f1 = open(regfile, "r")

    lines = f1.readlines()

    f1.close()

    flag = False
    found = False

    # reading reg file
    for line in lines:

        b1 = line.split("(")
        p = line.split(",")

        x1 = p[0]
        if b1[0] == "ellipse":

            x0 = "ellipse"
            x2 = x1[8:]
            flag = True
            found = True

        if flag is True:
            x3 = p[4]
            x4 = x3[:-2]

            v0 = x0

            v1 = float(x2)
            v2 = float(p[1])
            v3 = float(p[2])
            v4 = float(p[3])
            v5 = float(x4)

            flag = False

    if found:
        obj = v0
        xpos = v1
        ypos = v2
        rx = v3
        ry = v4
        angle = v5

        # avoids ds9 regions with Area = 0
        if rx < 1:
            rx = 1
        if ry < 1:
            ry = 1

        return obj, xpos, ypos, rx, ry, angle

    else:
        print("ellipse region was not found in file. Exiting.. ")
        sys.exit()

    return 0, 0, 0, 0


def GetPmax(image, mask, xmin, xmax, ymin, ymax):
    """gets the peak coordinates

    Given an image, and optionally a mask, it identifies
    the (x, y) pixels where the maximum value is located.

    Parameters
    ----------
    image: 2D-array of the image
    mask: 2D-array of the mask
    xmin, xmax, ymin, ymax : int, int, int, int
                            coordinates of the image section
                            where the maximum will be obtained

    Returns
    -------
    (xpos, ypos) : (x, y) coordinates of the maximum


    """
    xmin = int(xmin)
    xmax = int(xmax)
    ymin = int(ymin)
    ymax = int(ymax)

    chuckimg = image[ymin - 1 : ymax, xmin - 1 : xmax]
    if mask.any():
        chuckmsk = mask[ymin - 1 : ymax, xmin - 1 : xmax]

        invmask = np.logical_not(chuckmsk)

        invmask = invmask * 1

        chuckimg = chuckimg * invmask
    maxy, maxx = np.where(chuckimg == np.max(chuckimg))

    xpos = maxx[0] + xmin - 1
    ypos = maxy[0] + ymin - 1

    return (xpos, ypos)


def Ds9ell2Kronellv2(xpos, ypos, rx, ry, angle):
    """Converts DS9 ellipse parameters to geometrical parameters

    Parameters
    ----------
    obj : str
          object name
    xpos : float, x-center
    ypos : float, y-center
    rx : float, major or minor axis
    ry : float, minor or major axis
    angle : float, angular position

    Returns
    -------
    xx : float, x center position
    yy : float, y center positon
    Rkron : float, major axis
    theta: float, angular position measured from Y-axis
    e: float, ellipticity


    """

    if rx >= ry:

        q = ry / rx
        e = 1 - q
        Rkron = rx
        theta = angle - 90
        xx = xpos
        yy = ypos
    else:
        q = rx / ry
        e = 1 - q
        Rkron = ry
        theta = angle  # + 90
        xx = xpos
        yy = ypos

    return xx, yy, Rkron, theta, e


def GetExpTime(Image):
    """Get exposition time
    from the image header
    # repeated

    """

    try:
        hdu = fits.open(Image)
        exptime = hdu[0].header["EXPTIME"]
        hdu.close()
    except Exception:
        exptime = 1
    return float(exptime)


def readDataImg(conf):
    """reads the galaxy image and
    mask image

    # repeated

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


def PrintHeader(
    hdl,
    A,
    B,
    C,
    D,
    E,
    F,
    G,
    xlo,
    xhi,
    ylo,
    yhi,
    convx,
    convy,
    J,
    platedx,
    platedy,
    Op,
    P,
    S,
):
    """prints GALFIT header in a file

    # repeated

    """

    # k Check
    # print to filehandle
    # the header for GALFIT

    lineZa = "=========================================="
    lineZb = "========================================================\n"
    lineZ = lineZa + lineZb
    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}       # Input Data image (FITS file) \n".format(A)
    lineB = "B) {}   # Output data image block  \n".format(B)
    lineC = 'C) {}   # Sigma image name (made from data if blank or "none")  \n'.format(
        C
    )
    lineD = "D) {}   # Input PSF image and (optional) diffusion kernel\n".format(D)
    lineE = "E) {}   # PSF fine sampling factor relative to data  \n".format(E)
    lineF = "F) {}     # Bad pixel mask (FITS image or ASCII coord list) \n".format(F)
    lineG = "G) {}     # File with parameter constraints (ASCII file)  \n".format(G)
    lineH = "H) {} {} {} {}    # Image region to fit (xmin xmax ymin ymax)\n".format(
        xlo, xhi, ylo, yhi
    )
    lineI = "I) {} {}    # Size of the convolution box (x y)  \n".format(convx, convy)
    lineJ = "J) {}       # Magnitude photometric zeropoint  \n".format(J)
    lineK = "K) {} {}    # Plate scale (dx dy). [arcsec per pixel] \n".format(
        platedx, platedy
    )
    lineO = "O) {}    # Display type (regular, curses, both)   \n".format(Op)
    lineP = "P) {}    # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps\n".format(P)
    lineS = "S) {}    # Modify/create objects interactively?  \n".format(S)
    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS          \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:\n"
    line3 = "#    nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n"
    line4 = "#    ferrer, powsersic, sky, and isophote.   \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified:\n"
    line7 = "#  C0 (diskyness/boxyness),   \n"
    line8 = "#  Fn (n=integer, Azimuthal Fourier Modes),  \n"
    line9 = "#  R0-R10 (PA rotation, for creating spiral structures).  \n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number\n"
    line12 = "# column 2:     \n"
    line13 = "#  -- Parameter 0: the allowed functions are: sersic, nuker, expdisk\n"
    line14 = "#        edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky \n"
    line15 = "#        -- Parameter 1-10: value of the initial parameters \n"
    line16 = "#        -- Parameter C0:   For diskiness/boxiness \n"
    line17 = "#              <0 = disky   \n"
    line18 = "#              >0 = boxy    \n"
    line19 = "#       -- Parameter Z: Outputting image options, the options are:\n"
    line20 = "#  0 = normal, i.e. subtract final model from the data to create \n"
    line21 = "#            the residual image  \n"
    line22 = "#            1 = Leave in the model -- do not subtract from the data\n"
    line23 = "#            \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0)\n"
    line25 = "# column 4: comment  \n"
    line26 = " \n"

    line27a = "========================================="
    line27b = "=========================================================\n"
    line27 = line27a + line27b

    hdl.write(lineZ)
    hdl.write(lineX)
    hdl.write(lineA)
    hdl.write(lineB)
    hdl.write(lineC)
    hdl.write(lineD)
    hdl.write(lineE)
    hdl.write(lineF)
    hdl.write(lineG)
    hdl.write(lineH)
    hdl.write(lineI)
    hdl.write(lineJ)
    hdl.write(lineK)
    hdl.write(lineO)
    hdl.write(lineP)
    hdl.write(lineS)
    hdl.write(lineY)

    hdl.write(line0)
    hdl.write(line1)
    hdl.write(line2)
    hdl.write(line3)
    hdl.write(line4)
    hdl.write(line5)
    hdl.write(line6)
    hdl.write(line7)
    hdl.write(line8)
    hdl.write(line9)
    hdl.write(line10)

    hdl.write(line11)
    hdl.write(line12)
    hdl.write(line13)
    hdl.write(line14)
    hdl.write(line15)
    hdl.write(line16)
    hdl.write(line17)
    hdl.write(line18)
    hdl.write(line19)
    hdl.write(line20)
    hdl.write(line21)
    hdl.write(line22)
    hdl.write(line23)
    hdl.write(line24)
    hdl.write(line25)
    hdl.write(line26)
    hdl.write(line27)

    return True


def PrintSky(hdl, ncomp, sky, Z, fit):
    """Print GALFIT sky function to filehandle

    # repeated
    """

    # k Check

    line00 = "# Object number: {}   \n".format(ncomp)
    line01 = " 0)      sky       #    Object type  \n"
    line02 = " 1) {}     {}  # sky background  \n".format(sky, fit)
    line03 = " 2) 0.000      0    # dsky/dx (sky gradient in x) \n"
    line04 = " 3) 0.000      0     # dsky/dy (sky gradient in y)\n"
    line05 = " Z) {}        # Skip this model in output image?  \n".format(Z)
    line06 = "\n"
    line07a = "================================="
    line07b = "===============================================\n"
    line07 = line07a + line07b

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)

    return True


def PrintSersic(
    hdl, ncomp, xpos, ypos, magser, reser, nser, axratser, angleser, Z, fit
):
    """print GALFIT Sersic function to filehandle

    # not used
    # repeated

    """
    # k Check

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}  \n".format(ncomp)
    line01 = " 0)     sersic      #  Object type    \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  #  position x, y   \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}  #  total magnitude     \n".format(magser, fit)
    line04 = " 4) {:.2f}  {}   #  R_e  [Pixels]\n".format(reser, fit)
    line05 = " 5) {}  {}    #  Sersic exponent (deVauc=4, expdisk=1)\n".format(
        nser, fit
    )
    line06 = " 9) {:.2f}  {}  #  axis ratio (b/a)  \n".format(axratser, fit)
    line07 = "10) {:.2f}   {}  #  position angle (PA)  \n".format(angleser, fit)
    line08 = " Z) {}       #  Skip this model in output image? \n".format(Z)
    line09 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)
    hdl.write(line09)

    return True


def PrintGauss(hdl, ncomp, xpos, ypos, magass, fwhm, axratgass, anglegass, Z, fit):
    """print GALFIT GAUSS function to filehandle

    # not used
    # repeated

    """

    # print to filehandle
    # a sersic function given
    # by the parameters

    line00 = "# Object number: {}  \n".format(ncomp)
    line01 = " 0)     gaussian     #  Object type  \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {} #  position x, y \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}  #  total magnitude  \n".format(magass, fit)
    line04 = " 4) {:.2f}   {}  #  FWHM  [Pixels]  \n".format(fwhm, fit)
    line05 = " 9) {:.2f}   {}  #  axis ratio (b/a)  \n".format(axratgass, fit)
    line06 = "10) {:.2f}   {}   #  position angle (PA)  \n".format(anglegass, fit)
    line07 = " Z) {}    #  Skip this model in output image? \n".format(Z)
    line08 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)

    return True


def PrintExp(hdl, ncomp, xpos, ypos, magexp, rsexp, axratexp, angleexp, Z, fit):
    """print GALFIT exponential function to filehandle

    # repeated
    # not used

    """

    # k Check

    # print to filehandle
    # a exponential function given
    # by the parameters

    line00 = "# Object number: $ncomp    \n"
    line01 = " 0)     expdisk     # Object type       \n"
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  # position x, y  [pixel] \n".format(
        xpos, ypos, fit, fit
    )
    line03 = " 3) {:.2f}   {}   # total magnitude \n".format(magexp, fit)
    line04 = " 4) {:.2f}   {} # Rs  [Pixels]  \n".format(rsexp, fit)
    line05 = " 9) {:.2f}   {} # axis ratio (b/a) \n".format(axratexp, fit)
    line06 = "10) {:.2f}   {}  # position angle (PA)  \n".format(angleexp, fit)
    line07 = " Z) {} # Skip this model in output image? \n".format(Z)
    line08 = "\n"

    hdl.write(line00)
    hdl.write(line01)
    hdl.write(line02)
    hdl.write(line03)
    hdl.write(line04)
    hdl.write(line05)
    hdl.write(line06)
    hdl.write(line07)
    hdl.write(line08)

    return True


def GetAxis(Image):
    """Get number of rows and columns from the image"""

    i = 0  # index indicated where the data is located

    # if checkCompHDU(Image):
    #    i = 1

    hdu = fits.open(Image)

    ncol = hdu[i].header["NAXIS1"]
    nrow = hdu[i].header["NAXIS2"]

    hdu.close()

    return ncol, nrow
