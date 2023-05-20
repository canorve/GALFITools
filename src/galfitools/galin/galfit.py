#    EllipSect: An analysis tool for GALFIT output 
#    Copyright (C) 2022  Christopher Añorve 

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



from ellipsect.lib.libs import *

from ellipsect import *


import copy
from scipy.special import gamma, gammainc, gammaincinv

class GalHead():
    '''store the header of galfit file'''

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
    convx = 0
    convy = 0
    mgzpt = 25
    scale = 1
    scaley = 1
    display = "regular"
    P = 0

    #extra info  
    imgidx = "sci"
    flagidx = False
    num = 1
    flagnum = False
    exptime = 1
    tempmask = "tempmask.fits"

### class for Galfit components
class GalComps:
    '''stores the components of galfit file'''

    N = np.array([])
    NameComp = np.array([])        #0)
    PosX = np.array([])            #1)   
    PosY = np.array([])            #2)   
    Mag = np.array([])             #3)
    Rad = np.array([])             #4)
    Exp = np.array([])             #5)
    Exp2 = np.array([])            #6)  for moffat
    Exp3 = np.array([])            #7)  for moffat
                                   #8)  There is No 8 in any galfit model
    AxRat = np.array([])           #9)  AxisRatio
    PosAng = np.array([])          #10) position angle
    skip = np.array([])            #z)  skip model

    Active = np.array([])            #activate component  for galaxy

    # store the flags related to parameters
    PosXFree = np.array([])            #1)   
    PosYFree = np.array([])            #2)   
    MagFree = np.array([])             #3)
    RadFree = np.array([])             #4)
    ExpFree = np.array([])             #5)
    Exp2Free = np.array([])            #6)  for moffat
    Exp3Free = np.array([])            #7)  for moffat
                                   #8)  There is No 8 in any galfit model
    AxRatFree = np.array([])           #9)  AxisRatio
    PosAngFree = np.array([])          #10) position angle

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
    '''stores the value of the GALFIT file'''

    sky = 0 
    dskyx = 0
    dskyy = 0
    skip = 0 

    skyfree = 1 
    dskyxfree = 0 
    dskyyfree = 0 
 

class DataImg:

    img = np.array([[1,1],[1,1]])
    model = np.array([[1,1],[1,1]])
    imres = np.array([[1,1],[1,1]])

    mask = np.array([[1,1],[1,1]])
    sigma = np.array([[1,1],[1,1]])
    imsnr = np.array([[1,1],[1,1]])
    imchi = np.array([[1,1],[1,1]])
    impsf = np.array([[1,1],[1,1]])


def readDataImg(ellconf, galhead):

    dataimg = DataImg()

    #reading galaxy and model images from file
    errmsg="file {} does not exist".format(galhead.outimage)

    assert os.path.isfile(galhead.outimage), errmsg

    if ellconf.flagmodel == False:
        # hdu 1 => image   hdu 2 => model
        hdu = fits.open(galhead.outimage)
        dataimg.img = (hdu[1].data.copy()).astype(float)
        dataimg.model = (hdu[2].data.copy()).astype(float)
        dataimg.imres = (hdu[3].data.copy()).astype(float)
        hdu.close()

    else:
        hdu = fits.open(galhead.inputimage)

        if galhead.flagidx:
            if galhead.flagnum:
                dataimg.img = (hdu[galhead.imgidx, galhead.num].data).astype(float)
            else:    
                dataimg.img = (hdu[galhead.imgidx].data).astype(float)
        else:
            dataimg.img = (hdu[0].data).astype(float)

        hdu.close()

        hdu = fits.open(ellconf.inputmodel)
        dataimg.model = (hdu[0].data).astype(float)
        hdu.close()

        dataimg.imres = galpar.img - galpar.model


    ### reading mask image from file

    if galhead.tempmask is not None:

        errmsg="file {} does not exist".format(galhead.tempmask)
        assert os.path.isfile(galhead.tempmask), errmsg

        if ellconf.flagmodel == False:
            hdu = fits.open(galhead.tempmask)
            mask = hdu[0].data
            dataimg.mask=np.array(mask,dtype=bool)
            hdu.close()
        else:
            hdu = fits.open(galhead.maskimage)
            mask = hdu[0].data
            dataimg.mask=np.array(mask,dtype=bool)
            hdu.close()

    else:
        dataimg.mask = None
 

    return dataimg



class Galfit():

    def __init__(self, File: str):
        self.File = File


    def ReadHead(self) -> GalHead:
        '''reads header of galfit file'''
        inputf = self.File 

        galhead = GalHead() # class for header

        maskimage = ""
        #    skylevel=0

        GalfitFile = open(inputf,"r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        while index < len(lines):

            #================================================================================
            # IMAGE and GALFIT CONTROL PARAMETERS
            #A) tempfits/A2399-3-2.fits      # Input data image (FITS file)
            #B) A2399-215615.96-m073822.7-337-out.fits      # Output data image block
            #C) tempfits/none-3-2.fits      # Sigma image name (made from data if blank or "none")
            #D) psfs/PSF-1309-721.fits          # Input PSF image and (optional) diffusion kernel
            #E) 1                   # PSF fine sampling factor relative to data
            #F) mask-337            # Bad pixel mask (FITS image or ASCII coord list)
            #G) constraints         # File with parameter constraints (ASCII file)
            #H) 129  809  265  809  # Image region to fit (xmin xmax ymin ymax)
            #I) 60     60           # Size of the convolution box (x y)
            #J) 21.672              # Magnitude photometric zeropoint
            #K) 0.680  0.680        # Plate scale (dx dy)   [arcsec per pixel]
            #O) regular             # Display type (regular, curses, both)
            #P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

            line = lines[index]
            (tmp) = line.split()

            if len(tmp) > 1: #avoids empty options

                if tmp[0] == "A)":     # input image
                    galhead.inputimage=tmp[1]

                if tmp[0] == "B)":     # out image
                    galhead.outimage = tmp[1]

                if tmp[0] == "C)":   # sigma 
                    galhead.sigimage = tmp[1]

                if tmp[0] == "D)":  # psf file 
                    galhead.psfimage = tmp[1]

                if tmp[0] == "E)":  #psf sampling 
                    galhead.psfsamp = int(tmp[1])

                if tmp[0] == "F)":     # mask image
                    try:
                        galhead.maskimage = tmp[1]
                    except IndexError:
                        galhead.maskimage = "None"

                if tmp[0] == "G)":  # psf file 
                    galhead.constraints = tmp[1]


                if tmp[0] == "H)":     # region fit box
                    galhead.xmin = int(tmp[1])
                    galhead.xmax = int(tmp[2])
                    galhead.ymin = int(tmp[3])
                    galhead.ymax = int(tmp[4])

                if tmp[0] == "I)":     # convolution size 
                    galhead.convx = int(tmp[1])
                    galhead.convy = int(tmp[2])

                if tmp[0] == "J)":     # mgzpt
                    galhead.mgzpt = float(tmp[1])

                if tmp[0] == "K)":     # plate scale
                    galhead.scale = float(tmp[1])
                    galhead.scaley = float(tmp[2])

                if tmp[0] == "O)":     # display 
                    galhead.display = tmp[1]

                if tmp[0] == "P)":     # optional output 
                    galhead.P = int(tmp[1])


            index += 1


        ###### check for extension in name #####
        chars = set('[]') 
        numbers=set('1234567890')
        if any((c in chars) for c in galhead.inputimage): 
            print("Ext Found") 
            galhead.flagidx=True
            (filename,imgidxc) = galhead.inputimage.split("[")
            (imgidx,trash)=imgidxc.split("]")

            if any((n in numbers) for n in imgidx):
                galhead.flagnum=True
                (imgidx,num)=imgidx.split(",")
                num=int(num)

            galhead.inputimage=filename
            galhead.imgidx=imgidx
            galhead.num=num
      
        ####################

        return galhead


    def ReadComps(self) -> GalComps:
        '''reads all the components in the galfit file'''

        File = self.File

        galcomps = GalComps()
         
        GalfitFile = open(File,"r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        N = 0


        while index < len(lines):

            line = lines[index]
            (tmp) = line.split()

            #init values
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



            if (tmp[0] == "0)") and (tmp[1] != 'sky'):

                namec = tmp[1] 
                N = N + 1 
                NameComp = namec
                while (tmp[0] != "Z)"):

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":   # center
                        PosX = float(tmp[1])
                        PosY = float(tmp[2])
                        Xfree = int(tmp[3])
                        Yfree = int(tmp[4])
                    if tmp[0] == "3)" :    # mag 
                        Mag = float(tmp[1])
                        Magfree = int(tmp[2])
                    if tmp[0] == "4)" :
                        Rad = float(tmp[1])
                        Radfree = int(tmp[2])
                    if tmp[0] == "5)" : 
                        Exp = float(tmp[1])
                        Expfree = int(tmp[2])
                    if tmp[0] == "6)": 
                        Exp2 = float(tmp[1])
                        Exp2free = int(tmp[2])
                    if tmp[0] == "7)":  
                        Exp3 = float(tmp[1])
                        Exp3free = int(tmp[2])
                    if tmp[0] == "9)":
                        AxRat = float(tmp[1])
                        AxRatfree = int(tmp[2])
                    if tmp[0] == "10)":
                        PosAng = float(tmp[1])
                        PosAngfree = int(tmp[2])
                    if tmp[0] == "Z)":
                        skip=int(tmp[1])

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

        #filling the rest of arrays with zeros: 
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
        '''reads the sky value of the galfit file'''

        File = self.File    

        galsky = GalSky()

        GalfitFile = open(File,"r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in GalfitFile)
        lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        # remove lines containing only comments
        lines = (line.rstrip() for line in lines)
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)
        index = 0

        N = 0


        while index < len(lines):

            line = lines[index]
            (tmp) = line.split()

            #init values
            NameComp = "sky"
            sky = 0
            dskyx = 0
            dskyy = 0
            skip = 0

            skyfree = 1
            dskyxfree = 0
            dskyyfree = 0


            if (tmp[0] == "0)") and (tmp[1] == NameComp):

                namec = tmp[1] 

                while (tmp[0] != "Z)"):

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":   # center
                        sky = float(tmp[1])
                        skyfree = float(tmp[2])
                    if tmp[0] == "2)" :    # mag 
                        dskyx = float(tmp[1])
                        dskyxfree = int(tmp[2])
                    if tmp[0] == "3)" :
                        dskyy = float(tmp[1])
                        dskyyfree = int(tmp[2])
                    if tmp[0] == "Z)":
                        skip=int(tmp[1])

                #note if this work in this way or using append
                galsky.NameComp =  NameComp
                    
                galsky.sky = sky
                galsky.dskyx =  dskyx
                galsky.dskyy =  dskyy
                galskyskip =  skip

                galsky.skyfree =  skyfree
                galsky.dskyxfree = dskyxfree
                galsky.dskyyfree =  dskyyfree
     

            index += 1

        GalfitFile.close()

     
        return galsky


def numParFree(galcomps: GalComps) -> int:
    '''obtains the number of free parameters. This function does NOT 
    count the sky as a free param'''



    p1 = 0
    p2 = 0
    p3 = 0
    p4 = 0
    p5 = 0
    p6 = 0
    p7 = 0
    p8 = 0
    p9 = 0


    parmask1 = (galcomps.Active == True) & (galcomps.PosXFree == 1) 
    parmask2 = (galcomps.Active == True) & (galcomps.PosYFree == 1) 
    parmask3 = (galcomps.Active == True) & (galcomps.MagFree == 1) 
    parmask4 = (galcomps.Active == True) & (galcomps.RadFree == 1) 
    parmask5 = (galcomps.Active == True) & (galcomps.ExpFree == 1) 
    parmask6 = (galcomps.Active == True) & (galcomps.Exp2Free == 1) 
    parmask7 = (galcomps.Active == True) & (galcomps.Exp3Free == 1) 
    parmask8 = (galcomps.Active == True) & (galcomps.AxRatFree == 1) 
    parmask9 = (galcomps.Active == True) & (galcomps.PosAngFree == 1) 

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


    return pt




def numComps(galcomps: GalComps, name: str = 'none') -> int:
    '''obtains the number of components'''

    if name == 'none':
        nummask = (galcomps.Active == True) 
        N = galcomps.Active[nummask].size

    else:
        nummask = (galcomps.Active == True) & (galcomps.NameComp == name)
        N = galcomps.Active[nummask].size

    return N



def SelectGal(galcomps: GalComps, distmax: float, n_comp: int) -> GalComps:
    '''changes Flag to true for those components who belongs
        to the same galaxy of n_comp'''

    galcomps.Active.fill(False)

    idx = np.where(galcomps.N ==  n_comp)

    assert idx[0].size !=0, 'component not found'

    n_idx = idx[0].item(0)


    galcomps.Active[n_idx] = True #main component

    posx = galcomps.PosX[n_idx] 
    posy = galcomps.PosY[n_idx]      

    dx = galcomps.PosX - posx
    dy = galcomps.PosY - posy

    dist = np.sqrt((dx)**2 + (dy)**2)

    maskdist = (dist <= distmax ) 

    galcomps.Active[maskdist] = True

    return galcomps



### Sersic components 
class GetReff:
    '''class to obtain the effective radius for the whole galaxy'''

    def GetReSer(self, galhead: GalHead, comps: GalComps, eff: float, theta: float) -> float:


        maskgal = (comps.Active == True) 

        comps.Flux = 10**((galhead.mgzpt - comps.Mag)/2.5)

        totFlux = comps.Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + galhead.mgzpt

        a = 0.1
        b = comps.Rad[maskgal][-1] * 1000  # hope it doesn't crash

        Reff = self.solveSerRe(a, b, comps.Flux[maskgal], comps.Rad[maskgal], 
                comps.Exp[maskgal], comps.AxRat[maskgal], comps.PosAng[maskgal], totFlux, eff, theta)


        return Reff, totmag



    def solveSerRe(self, a: float, b: float, flux: list, rad: list, n: list, q: list, pa: list, totFlux: float, eff: float, theta: float) -> float:
        "return the Re of a set of Sersic functions. It uses Bisection"


        Re = bisect(self.funReSer, a, b, args=(flux, rad, n, q, pa, totFlux, eff, theta))

        return Re


    def funReSer(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, totFlux: float, eff: float, theta: float) -> float:
        

        fun = self.Ftotser(R, flux, rad, n, q, pa, theta) - totFlux*eff

        return fun
     
    def Ftotser(self, R: float, flux: list, rad: list, n: list, q: list, pa: list, theta: float) -> float:

        ftotR = self.Fser(R, flux, rad, n, q, pa, theta) 

        return ftotR.sum()



    def Fser(self, R: float, Flux: list, Re: list, n: list, q: list, pa: list, theta: float) -> float:
        '''sersic flux to a determined R'''
        
        k = gammaincinv(2*n, 0.5)

        Rcor = GetRadAng(R, q, pa, theta) 

        X = k*(Rcor/Re)**(1/n) 

        Fr = Flux*gammainc(2*n, X) ##esta funcion esta mal 
        
        return Fr



def conver2Sersic(galcomps: GalComps) -> GalComps:
    ''' function to convert exponential, gaussian params to Sersic params'''

    comps =  copy.deepcopy(galcomps)

    maskdev = comps.NameComp == "devauc"
    maskexp = comps.NameComp == "expdisk"
    maskgas = comps.NameComp == "gaussian"


    K_GAUSS = 0.6931471805599455 #constant k for gaussian
    K_EXP = 1.6783469900166612 # constant k for expdisk
    SQ2 = np.sqrt(2) 

    #for gaussian functions
    if maskgas.any():
        comps.Exp[maskgas] = 0.5 
        comps.Rad[maskgas] = comps.Rad[maskgas]/2.354 #converting to sigma 
        comps.Rad[maskgas] = SQ2*(K_GAUSS**0.5)*comps.Rad[maskgas] #converting to Re 


    #for de vaucouleurs
    if maskdev.any():
        comps.Exp[maskdev] = 4

    #for exponential disks
    if maskexp.any():
        comps.Exp[maskexp] = 1
        comps.Rad[maskexp] = K_EXP*comps.Rad[maskexp] #converting to Re


    return comps



def GetRadAng(R: float, q: list, pa: list, theta: float) -> float:
    '''Given an ellipse and an angle it returns the radius in angle direction. 
    Theta are the values for the galaxy and the others for every component'''



    #changing measured angle from y-axis to x-axis
    # and changing to rads:
    newpa = (pa + 90)*np.pi/180 #angle of every component
    theta = (theta + 90)*np.pi/180 #angle of direction of R 

    #bim = q * R

    ecc = np.sqrt(1 - q**2)

    alpha = theta - newpa #this is the direction 


    bell =  R*np.sqrt(1 - (ecc*np.cos(alpha))**2)


    aell = bell/q  #rad to evalue for every component



    return aell 







