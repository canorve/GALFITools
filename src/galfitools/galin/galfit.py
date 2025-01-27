import copy
import os

import numpy as np
from astropy.io import fits
from scipy.optimize import bisect
from scipy.special import gammainc, gammaincinv


class GalHead:
    """
    Data class to store the header info of the galfit file


    """

    inputimage = "none.fits"
    outimage = "none-out.fits"
    sigimage = "none"
    psfimage = "none"
    psfsamp = 1
    maskimage = "none"
    constraints = "none"
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    convx = 1
    convy = 1
    mgzpt = 25
    scale = 1
    scaley = 1
    display = "regular"
    P = 0

    # extra info
    imgidx = "sci"
    flagidx = False
    num = 1
    flagnum = False
    exptime = 1
    tempmask = "tempmask.fits"


class GalComps:
    """
    Data class to store the galfit components of the galfit file


    """

    N = np.array([])
    NameComp = np.array([])  # 0)
    PosX = np.array([])  # 1)
    PosY = np.array([])  # 2)
    Mag = np.array([])  # 3)
    Rad = np.array([])  # 4)
    Exp = np.array([])  # 5)
    Exp2 = np.array([])  # 6)  for moffat
    Exp3 = np.array([])  # 7)  for moffat
    # 8)  There is No 8 in any galfit model
    AxRat = np.array([])  # 9)  AxisRatio
    PosAng = np.array([])  # 10) position angle
    skip = np.array([])  # z)  skip model

    Active = np.array([])  # activate component  for galaxy

    # store the flags related to parameters
    PosXFree = np.array([])  # 1)
    PosYFree = np.array([])  # 2)
    MagFree = np.array([])  # 3)
    RadFree = np.array([])  # 4)
    ExpFree = np.array([])  # 5)
    Exp2Free = np.array([])  # 6)  for moffat
    Exp3Free = np.array([])  # 7)  for moffat
    # 8)  There is No 8 in any galfit model
    AxRatFree = np.array([])  # 9)  AxisRatio
    PosAngFree = np.array([])  # 10) position angle

    # computed parameters:
    Rad50 = np.array([])
    SerInd = np.array([])
    Rad50kpc = np.array([])
    Rad50sec = np.array([])
    Rad90 = np.array([])
    AbsMag = np.array([])
    Lum = np.array([])
    Flux = np.array([])
    PerLight = np.array([])
    me = np.array([])
    mme = np.array([])
    kser = np.array([])

    KronRad = np.array([])
    PetRad = np.array([])


class GalSky:
    """
    Data class to store the sky component of the GALFIT file


    """

    sky = 0
    dskyx = 0
    dskyy = 0
    skip = 0

    skyfree = 1
    dskyxfree = 0
    dskyyfree = 0


class DataImg:
    """
    Data class to store the galaxy, model and residual images
    of the output file of GALFIT


    """

    img = np.array([[1, 1], [1, 1]])
    model = np.array([[1, 1], [1, 1]])
    imres = np.array([[1, 1], [1, 1]])

    mask = np.array([[1, 1], [1, 1]])
    sigma = np.array([[1, 1], [1, 1]])
    imsnr = np.array([[1, 1], [1, 1]])
    imchi = np.array([[1, 1], [1, 1]])
    impsf = np.array([[1, 1], [1, 1]])


class Galfit:
    """
    Class to read GALFIT parameters file

    Attributes
    ----------
    File: str
         name of the GALFIT file

    Methods
    -------
    ReadHead : GalHead
        reads the header of the GALFIT file

    ReadComps : GalComps
        reads the galfit components of the GALFIT file

    ReadSky : GalSky
        reads the galfit sky components of the GALFIT file

    """

    def __init__(self, File: str):
        self.File = File

    def ReadHead(self) -> GalHead:

        inputf = self.File

        galhead = GalHead()  # class for header

        GalfitFile = open(inputf, "r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split("#", 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        while index < len(lines):

            # ================================================================
            # IMAGE and GALFIT CONTROL PARAMETERS
            # A) tempfits/A2399-3-2.fits      # Input data image (FITS file)
            # B) A2399-215615.96-m073822.7-337-out.fits      # Output data image block
            # C) tempfits/none-3-2.fits      # Sigma image name
            # D) psfs/PSF-1309-721.fits          # Input PSF image
            # E) 1                   # PSF fine sampling factor relative to data
            # F) mask-337            # Bad pixel mask (FITS image or ASCII coord list)
            # G) constraints         # File with parameter constraints (ASCII file)
            # H) 129  809  265  809  # Image region to fit (xmin xmax ymin ymax)
            # I) 60     60           # Size of the convolution box (x y)
            # J) 21.672              # Magnitude photometric zeropoint
            # K) 0.680  0.680        # Plate scale (dx dy)   [arcsec per pixel]
            # O) regular             # Display type (regular, curses, both)
            # P) 0                   # Choose: 0=optimize,1=model,2=imgblock,3=subcomps

            line = lines[index]
            (tmp) = line.split()

            if len(tmp) > 1:  # avoids empty options

                if tmp[0] == "A)":  # input image
                    galhead.inputimage = tmp[1]

                if tmp[0] == "B)":  # out image
                    galhead.outimage = tmp[1]

                if tmp[0] == "C)":  # sigma
                    galhead.sigimage = tmp[1]

                if tmp[0] == "D)":  # psf file
                    galhead.psfimage = tmp[1]

                if tmp[0] == "E)":  # psf sampling
                    galhead.psfsamp = int(tmp[1])

                if tmp[0] == "F)":  # mask image
                    try:
                        galhead.maskimage = tmp[1]
                    except IndexError:
                        galhead.maskimage = "None"

                if tmp[0] == "G)":  # psf file
                    galhead.constraints = tmp[1]

                if tmp[0] == "H)":  # region fit box
                    galhead.xmin = int(tmp[1])
                    galhead.xmax = int(tmp[2])
                    galhead.ymin = int(tmp[3])
                    galhead.ymax = int(tmp[4])

                if tmp[0] == "I)":  # convolution size
                    galhead.convx = int(tmp[1])
                    try:
                        galhead.convy = int(tmp[2])
                    except Exception:
                        galhead.convy = galhead.convx

                if tmp[0] == "J)":  # mgzpt
                    galhead.mgzpt = float(tmp[1])

                if tmp[0] == "K)":  # plate scale
                    galhead.scale = float(tmp[1])
                    try:
                        galhead.scaley = float(tmp[2])
                    except Exception:
                        galhead.scaley = galhead.scale

                if tmp[0] == "O)":  # display
                    galhead.display = tmp[1]

                if tmp[0] == "P)":  # optional output
                    galhead.P = int(tmp[1])

            index += 1

        #  check for extension in name
        chars = set("[]")
        numbers = set("1234567890")
        if any((c in chars) for c in galhead.inputimage):
            print("Ext Found")
            galhead.flagidx = True
            (filename, imgidxc) = galhead.inputimage.split("[")
            (imgidx, trash) = imgidxc.split("]")

            if any((n in numbers) for n in imgidx):
                galhead.flagnum = True
                (imgidx, num) = imgidx.split(",")
                num = int(num)

            galhead.inputimage = filename
            galhead.imgidx = imgidx
            galhead.num = num

        ####################

        return galhead

    def ReadComps(self) -> GalComps:
        File = self.File

        galcomps = GalComps()

        GalfitFile = open(File, "r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split("#", 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        N = 0

        while index < len(lines):

            line = lines[index]
            (tmp) = line.split()

            # init values
            NameComp = "none"
            PosX = 0
            PosY = 0
            Mag = 99
            Rad = 0
            Exp = 0
            Exp2 = 0
            Exp3 = 0
            AxRat = 1
            PosAng = 0
            skip = 0

            flagcomp = False

            Xfree = 1
            Yfree = 1
            Magfree = 1
            Radfree = 1
            Expfree = 1
            Exp2free = 1
            Exp3free = 1
            AxRatfree = 1
            PosAngfree = 1

            if (tmp[0] == "0)") and (tmp[1] != "sky"):

                namec = tmp[1]
                N = N + 1
                NameComp = namec
                while tmp[0] != "Z)":

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":  # center
                        PosX = float(tmp[1])
                        PosY = float(tmp[2])
                        Xfree = int(tmp[3])
                        Yfree = int(tmp[4])
                    if tmp[0] == "3)":  # mag
                        Mag = float(tmp[1])
                        Magfree = int(tmp[2])
                    if tmp[0] == "4)":
                        Rad = float(tmp[1])
                        Radfree = int(tmp[2])
                    if tmp[0] == "5)":
                        Exp = float(tmp[1])
                        Expfree = int(tmp[2])
                    if tmp[0] == "6)":
                        Exp2 = float(tmp[1])
                        try:
                            Exp2free = int(tmp[2])
                        except Exception:
                            Exp2free = 0
                    if tmp[0] == "7)":
                        Exp3 = float(tmp[1])
                        try:
                            Exp3free = int(tmp[2])
                        except Exception:
                            Exp3free = 0
                    if tmp[0] == "9)":
                        AxRat = float(tmp[1])
                        try:
                            AxRatfree = int(tmp[2])
                        except Exception:
                            AxRatfree = 0
                    if tmp[0] == "10)":
                        PosAng = float(tmp[1])
                        PosAngfree = int(tmp[2])
                    if tmp[0] == "Z)":
                        skip = int(tmp[1])

                galcomps.PosX = np.append(galcomps.PosX, PosX)
                galcomps.PosY = np.append(galcomps.PosY, PosY)
                galcomps.NameComp = np.append(galcomps.NameComp, NameComp)
                galcomps.N = np.append(galcomps.N, N)

                galcomps.Mag = np.append(galcomps.Mag, Mag)
                galcomps.Rad = np.append(galcomps.Rad, Rad)
                galcomps.Exp = np.append(galcomps.Exp, Exp)
                galcomps.Exp2 = np.append(galcomps.Exp2, Exp2)
                galcomps.Exp3 = np.append(galcomps.Exp3, Exp3)
                galcomps.AxRat = np.append(galcomps.AxRat, AxRat)
                galcomps.PosAng = np.append(galcomps.PosAng, PosAng)
                galcomps.skip = np.append(galcomps.skip, skip)
                galcomps.Active = np.append(galcomps.Active, flagcomp)

                galcomps.PosXFree = np.append(galcomps.PosXFree, Xfree)
                galcomps.PosYFree = np.append(galcomps.PosYFree, Yfree)
                galcomps.MagFree = np.append(galcomps.MagFree, Magfree)
                galcomps.RadFree = np.append(galcomps.RadFree, Radfree)
                galcomps.ExpFree = np.append(galcomps.ExpFree, Expfree)
                galcomps.Exp2Free = np.append(galcomps.Exp2Free, Exp2free)
                galcomps.Exp3Free = np.append(galcomps.Exp3Free, Exp3free)
                galcomps.AxRatFree = np.append(galcomps.AxRatFree, AxRatfree)
                galcomps.PosAngFree = np.append(galcomps.PosAngFree, PosAngfree)

            index += 1

        GalfitFile.close()

        # filling the rest of arrays with zeros:
        arsiz = len(galcomps.N)

        galcomps.Rad50 = np.zeros(arsiz)
        galcomps.SerInd = np.zeros(arsiz)
        galcomps.Rad50kpc = np.zeros(arsiz)
        galcomps.Rad50sec = np.zeros(arsiz)
        galcomps.Rad90 = np.zeros(arsiz)
        galcomps.AbsMag = np.zeros(arsiz)
        galcomps.Lum = np.zeros(arsiz)
        galcomps.Flux = np.zeros(arsiz)
        galcomps.PerLight = np.zeros(arsiz)
        galcomps.me = np.zeros(arsiz)
        galcomps.mme = np.zeros(arsiz)
        galcomps.kser = np.zeros(arsiz)
        galcomps.KronRad = np.zeros(arsiz)
        galcomps.PetRad = np.zeros(arsiz)

        return galcomps

    def ReadSky(self) -> GalSky:

        File = self.File

        galsky = GalSky()

        GalfitFile = open(File, "r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split("#", 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        while index < len(lines):

            line = lines[index]
            (tmp) = line.split()

            # init values
            NameComp = "sky"
            sky = 0
            dskyx = 0
            dskyy = 0
            skip = 0

            skyfree = 1
            dskyxfree = 0
            dskyyfree = 0

            if (tmp[0] == "0)") and (tmp[1] == NameComp):

                while tmp[0] != "Z)":

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":  # center
                        sky = float(tmp[1])
                        skyfree = float(tmp[2])
                    if tmp[0] == "2)":  # mag
                        dskyx = float(tmp[1])
                        dskyxfree = int(tmp[2])
                    if tmp[0] == "3)":
                        dskyy = float(tmp[1])
                        dskyyfree = int(tmp[2])
                    if tmp[0] == "Z)":
                        skip = int(tmp[1])

                # note if this work in this way or using append
                galsky.NameComp = NameComp

                galsky.sky = sky
                galsky.dskyx = dskyx
                galsky.dskyy = dskyy
                galsky.skip = skip

                galsky.skyfree = skyfree
                galsky.dskyxfree = dskyxfree
                galsky.dskyyfree = dskyyfree

            index += 1

        GalfitFile.close()

        return galsky


def numParFree(galcomps: GalComps) -> int:
    """Obtains the number of free parameters

    Count the number of free parameters of the surface brightness model
    of the GALFIT file.

    Parameters
    ----------
    galcomps : GalComps Data Class defined above


    Returns
    -------
    pt : int
        Number of free parameters


    See Also
    --------
    numParSkyFree: Obtains the number of free parameters of the sky
                    component


    Notes
    -----
    It only counts the components with galcomps.Active == True
    This prevents to count components of other galaxies that
    were simultanously fitted.

    Sky component is ignored.


    """

    p1 = 0
    p2 = 0
    p3 = 0
    p4 = 0
    p5 = 0
    p6 = 0
    p7 = 0
    p8 = 0
    p9 = 0

    parmask1 = (galcomps.Active == 1) & (galcomps.PosXFree == 1)
    parmask2 = (galcomps.Active == 1) & (galcomps.PosYFree == 1)
    parmask3 = (galcomps.Active == 1) & (galcomps.MagFree == 1)
    parmask4 = (galcomps.Active == 1) & (galcomps.RadFree == 1)
    parmask5 = (galcomps.Active == 1) & (galcomps.ExpFree == 1)
    parmask6 = (galcomps.Active == 1) & (galcomps.Exp2Free == 1)
    parmask7 = (galcomps.Active == 1) & (galcomps.Exp3Free == 1)
    parmask8 = (galcomps.Active == 1) & (galcomps.AxRatFree == 1)
    parmask9 = (galcomps.Active == 1) & (galcomps.PosAngFree == 1)

    if parmask1.any():
        p1 = np.sum(galcomps.PosXFree[parmask1])
    if parmask2.any():
        p2 = np.sum(galcomps.PosYFree[parmask2])
    if parmask3.any():
        p3 = np.sum(galcomps.MagFree[parmask3])
    if parmask4.any():
        p4 = np.sum(galcomps.RadFree[parmask4])
    if parmask5.any():
        p5 = np.sum(galcomps.ExpFree[parmask5])
    if parmask6.any():
        p6 = np.sum(galcomps.Exp2Free[parmask6])
    if parmask7.any():
        p7 = np.sum(galcomps.Exp3Free[parmask7])
    if parmask8.any():
        p8 = np.sum(galcomps.AxRatFree[parmask8])
    if parmask9.any():
        p9 = np.sum(galcomps.PosAngFree[parmask9])

    pt = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9

    return int(pt)


def numParSkyFree(galsky: GalSky) -> int:
    """Obtains the number of free parameters of the sky component

    Count the number of free parameters of the sky component
    of the GALFIT file.

    Parameters
    ----------
    galsky: GalSky data Class defined above


    Returns
    -------
    pt : int
        Number of free parameters


    See Also
    --------
    numParFree : Obtains the number of free parameters


    """

    p1 = 0
    p2 = 0
    p3 = 0

    parmask1 = galsky.skyfree == 1
    parmask2 = galsky.dskyxfree == 1
    parmask3 = galsky.dskyyfree == 1

    if parmask1:
        p1 = 1
    if parmask2:
        p2 = 1
    if parmask3:
        p3 = 1

    pt = p1 + p2 + p3

    return int(pt)


def numComps(galcomps: GalComps, name: str) -> int:
    """gets the number of components for a galaxy.

    Given the number of components of a GALFIT file,
    it discerns which components belong to the galaxy. This
    is useful when simultaneous fitting was used.

    Parameters
    ----------
    galcomps : GalComps data class defined above
    name :  str
            indicates which model components to count:
            all, sersic, expdisk, gaussian, devauc

    Returns
    -------
    N : Number of components of the galaxy


    """

    if name == "all":
        nummask = (galcomps.Active == 1) & (
            (galcomps.NameComp == "sersic")
            | (galcomps.NameComp == "expdisk")
            | (galcomps.NameComp == "gaussian")
            | (galcomps.NameComp == "devauc")
        )

    else:
        nummask = (galcomps.Active == 1) & (galcomps.NameComp == name)

    N = galcomps.Active[nummask].size

    return N


def SelectGal(galcomps: GalComps, distmax: float, n_comp: int) -> GalComps:
    """Selects which components belong to the galaxy.

    From the data class GalComps, it changes the flag Active to
    True for those components where distance is less than distmax
    from component n_comp.

    This is useful when there are simultaneous
    fitting of galaxies in a single GALFIT file.


    Parameters
    ----------
    galcomps : GalComps data class defined above
    distmax : float
              maximum distance among components
    n_comp : int
            The center of the galaxy will be determined from this
            component (n_comp). The distances of the remaining
            components will be calculated relative to this center.

    Returns
    -------
    galcomps : GalComps data class with flag Active = True for
              those components which belong to the galaxy selected
              by n_comp.

    Notes
    ------
    If a galaxy is made up of three components, it does not matter which
    of the three you select as a center. As long these belong to the
    same galaxy.


    """

    galcomps.Active.fill(False)

    idx = np.where(galcomps.N == n_comp)

    assert idx[0].size != 0, "component not found"

    n_idx = idx[0].item(0)

    if distmax > 10:
        print("Warning: maximum distance among components is greater than 10")
        print(
            "Equations' solutions only apply for components that share the same center"
        )

    galcomps.Active[n_idx] = True  # main component

    posx = galcomps.PosX[n_idx]
    posy = galcomps.PosY[n_idx]

    dx = galcomps.PosX - posx
    dy = galcomps.PosY - posy

    dist = np.sqrt((dx) ** 2 + (dy) ** 2)

    maskdist = dist <= distmax

    galcomps.Active[maskdist] = True

    return galcomps


def conver2Sersic(galcomps: GalComps) -> GalComps:
    """Converts the exponential, gaussian or De Vaucouleurs to Sersic

    Using the format of GALFIT with the GalComps data class,
    it converts the parameters of the functions exponential,
    gaussian or De Vaucouleurs to parameters of Sersic function.


    Parameters
    -----------
    galcomps : GalComps data class defined above


    Returns
    -------
    galcomps : GalComps data class with the functions converted to Sersic
              parameters

    """

    comps = copy.deepcopy(galcomps)

    maskdev = comps.NameComp == "devauc"
    maskexp = comps.NameComp == "expdisk"
    maskgas = comps.NameComp == "gaussian"

    K_GAUSS = 0.6931471805599455  # constant k for gaussian
    K_EXP = 1.6783469900166612  # constant k for expdisk
    SQ2 = np.sqrt(2)

    # for gaussian functions
    if maskgas.any():
        comps.Exp[maskgas] = 0.5
        comps.Rad[maskgas] = comps.Rad[maskgas] / 2.354  # converting to sigma
        comps.Rad[maskgas] = (
            SQ2 * (K_GAUSS**0.5) * comps.Rad[maskgas]
        )  # converting to Re

    # for de vaucouleurs
    if maskdev.any():
        comps.Exp[maskdev] = 4

    # for exponential disks
    if maskexp.any():
        comps.Exp[maskexp] = 1
        comps.Rad[maskexp] = K_EXP * comps.Rad[maskexp]  # converting to Re

    return comps


def GetRadAng(R: float, q: list, pa: list, theta: float) -> float:
    """Obtains the radius along the specified angular direction.

    Given the angular position and axis ratio of the ellipse
    it gives the radius of the ellipse in the theta direction.

    If theta is in the major axis direction it returns the
    major axis. On the other hand, if theta is in the direction
    of the minor axis it returns the minor axis.


    Parameters
    ----------
    R : float
        Radius along the major axis
    q : list
        Axis ratio of ellipse
    pa : list
        Angular position of ellipse measured from Y-axis (same as GALFIT)
    theta : float
            angle that defines the direction

    Returns
    -------
    aell : float
            radius along the angle specified by theta

    """

    # changing measured angle from y-axis to x-axis
    # and changing to rads:
    newpa = (pa + 90) * np.pi / 180  # angle of every component
    theta = (theta + 90) * np.pi / 180  # angle of direction of R

    # bim = q * R

    ecc = np.sqrt(1 - q**2)

    alpha = theta - newpa  # this is the direction

    bell = R * np.sqrt(1 - (ecc * np.cos(alpha)) ** 2)

    aell = bell / q  # rad to evalue for every component

    return aell


def galfitLastFit(directory: str) -> str:
    """determine the last fit completed by GALFIT

    Every time GALFIT produces a successful fit it produces
    a file called galfit.XX, where the XX represent a number
    that increases every time where GALFIT

    Parameters
    ----------
    directory : str
                directory containing the galfit files

    Returns
    -------
    max_file : str
                last fitting model file produced by GALFIT

    """
    # Directory containing the files
    # directory = "." #actual directory

    # List all files in the directory
    files = os.listdir(directory)

    # Filter files that match the prefix 'galfit.'
    galfit_files = [f for f in files if f.startswith("galfit.") and f[7:].isdigit()]

    # Extract the numerical suffix and find the maximum
    max_file = max(galfit_files, key=lambda x: int(x[7:]))

    return max_file


def galPrintHeader(hdl: str, galhead: GalHead) -> bool:
    """prints GALFIT header to a file

    Given a file handler, prints the GALFIT header parameters of
    to a file


    Parameters
    ----------
    hdl : str
        file handler where the header information will be print
    galhead : GalHead data class defined above

    Returns
    -------
        bool

    """

    A = galhead.inputimage
    B = galhead.outimage
    C = galhead.sigimage
    D = galhead.psfimage
    E = galhead.psfsamp
    F = galhead.maskimage
    G = galhead.constraints
    xlo = galhead.xmin
    xhi = galhead.xmax
    ylo = galhead.ymin
    yhi = galhead.ymax
    convx = galhead.convx
    convy = galhead.convy
    J = galhead.mgzpt
    platedx = galhead.scale
    platedy = galhead.scaley
    varO = galhead.display
    P = galhead.P
    # S = 0

    # k Check
    # print to filehandle
    # the header for GALFIT

    lineZa = "================================================="
    lineZb = "=================================================\n"
    lineZ = lineZa + lineZb
    lineX = "# IMAGE PARAMETERS \n"
    lineA = "A) {}    # Input Data image (FITS file)            \n".format(A)
    lineB = "B) {}    # Output data image block               \n".format(B)
    lineC = "C) {}    # Sigma image name   \n".format(C)
    lineD = "D) {}    # Input PSF image and (optional) diffusion kernel \n".format(D)
    lineE = "E) {}    # PSF fine sampling factor relative to data    \n".format(E)
    lineF = "F) {}    # Bad pixel mask (FITS image or ASCII coord list) \n".format(F)
    lineG = "G) {}    # File with parameter constraints (ASCII file)     \n".format(G)
    lineH = "H) {} {} {} {}   # Image region to fit (xmin xmax ymin ymax)  \n".format(
        xlo, xhi, ylo, yhi
    )
    lineI = "I) {} {}  # Size of the convolution box (x y)      \n".format(convx, convy)
    lineJ = "J) {}     # Magnitude photometric zeropoint          \n".format(J)
    lineK = "K) {} {}  # Plate scale (dx dy). [arcsec per pixel]\n".format(
        platedx, platedy
    )
    lineO = "O) {}     # Display type (regular, curses, both) \n".format(varO)
    lineP = "P) {}     # Choose 0=optimize, 1=model, 2=imgblock, 3=subcomps \n".format(
        P
    )

    lineY = " \n"

    line0 = "# INITIAL FITTING PARAMETERS                                \n"
    line1 = "# \n"
    line2 = "#   For object type, allowed functions are:                     \n"
    line3 = "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,\n"
    line4 = "#       ferrer, powsersic, sky, and isophote.               \n"
    line5 = "# \n"
    line6 = "#  Hidden parameters will only appear when they're specified: \n"
    line7 = "#      C0 (diskyness/boxyness),                             \n"
    line8 = "#      Fn (n=integer, Azimuthal Fourier Modes),             \n"
    line9 = "#      R0-R10 (PA rotation, for creating spiral structures). \n"
    line10 = "# \n"

    line11 = "# column 1:  Parameter number\n"
    line12 = "# column 2:           \n"
    line13 = "#      -- Parameter 0: the allowed functions are: sersic,nuker,expdisk\n"
    line14 = "#      edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky\n"
    line15 = "#      -- Parameter 1-10: value of the initial parameters \n"
    line16 = "#      -- Parameter C0:   For diskiness/boxiness\n"
    line17 = "#                     <0 = disky      \n"
    line18 = "#                     >0 = boxy    \n"
    line19 = "#      -- Parameter Z:  the options are:  \n"
    line20 = "#                     0 = normal,  \n"
    line21 = "#                     the residual image\n"
    line22 = "#                     1 = Leave in the model \n"
    line23 = "#                       \n"
    line24 = "# column 3: allow parameter to vary (yes = 1, no = 0) \n"
    line25 = "# column 4: comment     \n"
    line26 = " \n"

    line27a = "==================================================="
    line27b = "===============================================\n"
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
    # hdl.write(lineS)
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


def galPrintComp(hdl: str, ncomp: int, idx: int, galcomps: GalComps) -> bool:
    """prints GALFIT component parameters to a file

    Given a file handler, prints the GALFIT component parameters of
    one model the surface brightness models to a GALFIT file


    Parameters
    ----------
    hdl : str
        file handler where the component parameter information will be print
    ncomp : int
        prints the number of component in the GALFIT file
    idx : int
        index to indicate which component of the galcomps data class will be print

    galcomps : GalComps data class defined above

    Returns
    -------
        bool

    """

    # print to filehandle

    line00 = "# Object number: {}   \n".format(ncomp)
    line01 = " 0)   {}   #  Object type      \n".format(galcomps.NameComp[idx])
    line02 = " 1) {:.2f}  {:.2f}  {}  {}  #  position x, y  [pixel] \n".format(
        galcomps.PosX[idx],
        galcomps.PosY[idx],
        galcomps.PosXFree[idx],
        galcomps.PosYFree[idx],
    )
    line03 = " 3) {:.2f}  {}    #  total magnitude  \n".format(
        galcomps.Mag[idx], galcomps.MagFree[idx]
    )
    line04 = " 4) {:.2f}  {}    #  R_e     [Pixels] \n".format(
        galcomps.Rad[idx], galcomps.RadFree[idx]
    )
    line05 = " 5) {}     {}   #  Sersic exponent (deVauc=4, expdisk=1) \n".format(
        galcomps.Exp[idx], galcomps.ExpFree[idx]
    )
    line06 = " 6)  {}  {}      #  ---------------- \n".format(
        galcomps.Exp2[idx], galcomps.Exp2Free[idx]
    )
    line07 = " 7)  {}  {}      #  ---------------- \n".format(
        galcomps.Exp3[idx], galcomps.Exp3Free[idx]
    )
    line08 = (
        " 8)  0.0000       0   #  ----------------                                \n"
    )
    line09 = " 9) {:.2f}   {}   #  axis ratio (b/a)  \n".format(
        galcomps.AxRat[idx], galcomps.AxRatFree[idx]
    )
    line10 = "10) {:.2f}    {}  #  position angle (PA)   \n".format(
        galcomps.PosAng[idx], galcomps.PosAngFree[idx]
    )
    lineZ = " Z) {}         #  Skip this model in output image?  \n".format(
        galcomps.skip[idx]
    )
    line11 = "\n"

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
    hdl.write(line10)
    hdl.write(lineZ)
    hdl.write(line11)

    return True


def galPrintSky(hdl: str, ncomp: int, galsky: GalSky) -> bool:
    """prints GALFIT sky component parameters to a file

    Given a file handler, prints the GALFIT sky component parameters
    to a GALFIT file

    Parameters
    ----------
    hdl : str
        file handler where the sky component parameter information will be print
    ncomp : int
        prints the number of component in the GALFIT file
    galsky : GalSky data class defined above

    Returns
    -------
        bool

    """

    line00 = "# Object number: {}                 \n".format(ncomp)
    line01 = " 0)      sky            #    Object type   \n"
    line02 = " 1) {:.2f}   {}   # sky background        [ADU counts]  \n".format(
        galsky.sky, galsky.skyfree
    )
    line03 = " 2) {:.2f}  {}    # dsky/dx (sky gradient in x) \n".format(
        galsky.dskyx, galsky.dskyxfree
    )
    line04 = " 3) {:.2f}  {}    # dsky/dy (sky gradient in y) \n".format(
        galsky.dskyy, galsky.dskyyfree
    )
    line05 = " Z) {}  # Skip this model in output image?  (yes=1, no=0) \n".format(
        galsky.skip
    )
    line06 = "\n"
    line07a = "================================================"
    line07b = "================================\n"
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
