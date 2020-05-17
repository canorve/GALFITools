#! /usr/bin/env python3

import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path  
import scipy
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mimetypes

from astropy.io.votable import parse
from mgefit.sectors_photometry import sectors_photometry

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,NullFormatter,
                               AutoMinorLocator,LogLocator,LinearLocator,AutoLocator)

### Dictionary for Absolute mag of the Sun taken from Willmer 2018
SunMag = {
        "U":5.61,"B":5.44,"V": 4.81,"R":4.43,"I":4.1,
        "J": 3.67,"H": 3.32,"K": 3.27,
        "u":5.49,"g":5.23,"r": 4.53,"i":4.19,"z":4.01,
        "L":3.26
        } 

####



def main():

    if (len(sys.argv[1:]) == 0):
        print ('Missing arguments')
        print ("Usage:\n %s [GALFITOutputFile] [-options] " % (sys.argv[0]))
        print ("use help to display more information about 'options' arguments: ")
        print ("%s -help " % (sys.argv[0]))

        sys.exit()

    #class for saving user's parameters
    params=InputParams()

    OptionHandleList = ['-logx', '-q', '-pa','-comp','-pix','-ranx','-rany','-grid','-dpi','-sbout','-noplot','-minlevel','-sectors','-phot','-object','-filter','-snr','-help','-checkimg','-noned','-distmod','-magcor','-scalekpc','-sbdim']
    options = {}
    for OptionHandle in OptionHandleList:
        options[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)] if OptionHandle in sys.argv else None
    if options['logx'] != None:
        params.flaglogx=True
        print("X axis is logarithm")
    if options['q'] != None:
        params.flagq=True
    if options['pa'] != None:
        params.flagpa=True
    if options['pix'] != None:
        params.flagpix=True
    if options['ranx'] != None:
        params.flagranx[0]=True
    if options['rany'] != None:
        params.flagrany[0]=True
    if options['grid'] != None:
        params.flagrid=True
    if options['dpi'] != None:
        params.flagdpi=True
    if options['comp'] != None:
        params.flagcomp=True
        print("Plotting subcomponents ")
    if options['noplot'] != None:
        params.flagnoplot=True
        print("images will not be displayed")
    if options['sbout'] != None:
        params.flagsbout=True
        print("surface brightness output file will be created")
    if options['phot'] != None:
        params.flagphot=True
        print("output photometry file will be created")
    if options['minlevel'] != None:
        params.flagminlevel=True
    if options['sectors'] != None:
        params.flagsectors=True
    if options['object'] != None:
        params.flagobj=True
    if options['filter'] != None:
        params.flagband=True
    if options['snr'] != None:
        params.flagsnr=True
    if options['checkimg'] != None:
        params.flagcheck=True
    if options['noned'] != None:
        params.flagned=True
    if options['distmod'] != None:
        params.flagmod=True
    if options['magcor'] != None:
        params.flagmag=True
    if options['scalekpc'] != None:
        params.flagscale=True
    if options['sbdim'] != None:
        params.flagdim=True



    if options['help'] != None:
        Help()




    ################## search arguments after the option:
    if params.flagpa == True:
        opt={}
        OptionHandle="-pa"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.parg=np.float(opt['pa'])

    if params.flagq == True:
        opt={}
        OptionHandle="-q"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.qarg=np.float(opt['q'])

    if params.flagranx[0] == True:
        opt={}
        OptionHandle="-ranx"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]

        params.rangex=opt["ranx"]
        if "-" in params.rangex:
            params.flagranx[1] = True
            params.ranx=opt['ranx']
        else:
            params.ranx=np.float(opt['ranx'])

    if params.flagrany[0]== True:
        opt={}
        OptionHandle="-rany"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]

        params.rangey=opt["rany"]
        if "-" in params.rangey:
            params.flagrany[1] = True
            params.rany=opt['rany']
        else:
            params.rany=np.float(opt['rany'])

    if params.flagdpi == True:
        opt={}
        OptionHandle="-dpi"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.dpival=np.int(opt['dpi'])

    if params.flagnoplot == True:
        params.dplot=False

    if params.flagminlevel == True:
        opt={}
        OptionHandle="-minlevel"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.minlevel=np.int(opt['minlevel'])

    if params.flagsectors == True:
        opt={}
        OptionHandle="-sectors"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.minlevel=np.int(opt['sectors'])

    if params.flagobj == True:
        opt={}
        OptionHandle="-object"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.objname=np.str(opt['object'])

    if params.flagband == True:
        opt={}
        OptionHandle="-filter"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.band=np.str(opt['filter'])


    if params.flagmod == True:
        opt={}
        OptionHandle="-distmod"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.InDistMod=np.float(opt['distmod'])


    if params.flagmag == True:
        opt={}
        OptionHandle="-magcor"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.InMagCor=np.float(opt['magcor'])


    if params.flagscale == True:
        opt={}
        OptionHandle="-scalekpc"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.InScale=np.float(opt['scalekpc'])


    if params.flagdim == True:
        opt={}
        OptionHandle="-sbdim"
        opt[OptionHandle[1:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        params.InSbDim=np.float(opt['sbdim'])



    params.galfile= sys.argv[1]

    print("angle in multi-plot is measured from the galaxy's major axis ")

    ########  End of parameter reading #############
    ################################################
    ################################################


    #class for GALFIT's parameters
    galpar=GalfitParams()

    #class for GALFIT's components
    galcomps=GalfitComps()


    ######################################
    ####### Read Galfit File #############
    #  xc,yc,q,ang,skylevel,scale,file,mgzpt,exptime,mask=ReadGALFITout(params.galfile,galpars)
    ReadGALFITout(params.galfile,galpar)

    ######################################
    ######################################

    if params.flagq == True:
        galpar.q=params.qarg

    if params.flagpa == True:
        galpar.ang=params.parg


    str = "q = {} is used ".format(galpar.q)
    print(str)

    str = "pa = {} is used ".format(galpar.ang)
    print(str)

    str = "sky = {} ".format(galpar.skylevel)
    print(str)

    ##
    str = "dpi = {} for plots ".format(params.dpival)
    print(str)
    ##

    ##
    str = "minlevel = {} ".format(params.minlevel)
    print(str)

    ##
    str = "number of sectors = {}  ".format(params.sectors)
    print(str)



    (tmp)=galpar.outimage.split(".")

    params.namefile=tmp[0]

    # names for the different png

    params.namepng=params.namefile + ".png"
    params.namesec=params.namefile + "-gal.png"
    params.namemod=params.namefile + "-mod.png"
    params.namemul=params.namefile + "-mul.png"
    params.namesub=params.namefile + "-sub.fits"

    params.namesig=params.namefile + "-sig.fits"


    params.sboutput=params.namefile + "-sbout"
    params.output=params.namefile + "-out.txt"

    params.namened=params.namefile + "-ned.xml"



    params.namesnr=params.namefile + "-snr.fits"

    params.namecheck=params.namefile + "-check.fits"



    if params.flagsbout == True: 
        msg="surface brightness output file: {} ".format(params.sboutput+".txt")
        print(msg)

    if params.flagphot == True: 
        msg="output photometry file: {} ".format(params.output)
        print(msg)


    # read all the object components of the model in galfit.XX
    ReadNComp(params.galfile,galpar.xc,galpar.yc,galcomps)
    print("Number of components = ",len(galcomps.N))




    #reading galaxy and model images from file
    errmsg="file {} does not exist".format(galpar.outimage)

    assert os.path.isfile(galpar.outimage), errmsg

    # hdu 1 => image   hdu 2 => model
    hdu = fits.open(galpar.outimage)
    galpar.img = (hdu[1].data).astype(float)
    galpar.model = (hdu[2].data).astype(float)
    galpar.imres = (hdu[3].data).astype(float)
    hdu.close()


    # removing background from galaxy and model images 
    galpar.img = galpar.img - galpar.skylevel
    galpar.model = galpar.model - galpar.skylevel


    ### reading mask image from file

    if galpar.tempmask is not None:

        errmsg="file {} does not exist".format(galpar.tempmask)
        assert os.path.isfile(galpar.tempmask), errmsg

        hdu = fits.open(galpar.tempmask)
        mask = hdu[0].data
        galpar.mask=np.array(mask,dtype=bool)
        hdu.close()
    else:
        galpar.mask=None
    ####


    #   numsectors=19
    #   numsectors=15
    numsectors=params.sectors

    # minlevel=-100  # minimun value for sky
    # minlevel=15  # minimun value for sky
    minlevel=params.minlevel  # minimun value for sky

    # initial values for image matrixes 
    sectgalax=sectmodel=sectcomps=[]

    #call to sectors_photometry for galaxy and model
    sectgalax,sectmodel=SectPhot(galpar, params, n_sectors=numsectors, minlevel=minlevel)

    
    if params.flagcomp:
        #change to minlevel =0 ?
        sectcomps=SectPhotComp(galpar, params, galcomps, n_sectors=numsectors, minlevel=minlevel)



    print("creating plots..")

    limx,limy=EllipSectors(params, galpar, galcomps, sectgalax,sectmodel, sectcomps,n_sectors=numsectors)


    ##############################################
    ##############################################
    ##############################################

    if params.dplot:
        plt.pause(1.5)
    plt.savefig(params.namepng,dpi=params.dpival)
    plt.close()


    ########################################################
    ################ Multiplots: ###########################
    ########################################################

    print("creating multi-plots..")

    #MulEllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectmulcomps)

    MulEllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps)


    if params.dplot:
        plt.pause(1.5)

    plt.savefig(params.namemul,dpi=params.dpival)
    plt.close()


    ########################################################
    ############ Computing output photometry ###############
    ########################################################


    if params.flagphot:
        print("Computing output photometry ... ")
        OutPhot(params, galpar, galcomps, sectgalax, sectmodel, sectcomps)


    if galpar.tempmask != None:
        os.remove(galpar.tempmask) # removing temp mask file


    ##############       #############
    ##############  END  #############
    ##############       #############

    #     ______________________________________________________________________
    #    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
    #   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
    #   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/


### class for parameters
class InputParams:

    #flags
    flaglogx=False
    flagq=False
    flagpa=False
    flagcomp=False
    flagpix=False
    flagranx=[False,False]
    flagrany=[False,False]
    flagnoplot=False
    flagrid=False
    flagdpi=False
    flagsbout=False
    flagphot=False
    flagminlevel=False
    flagsectors=False
    flagobj=False
    flagband=False
    flagweb=False # to check connection to ned
    flagsnr = False

    flagcheck=False
    flagned=False

    flagmod=False
    flagmag=False
    flagscale=False
    flagdim=False
   

    #init
    qarg=1
    parg=0
    ranx=1
    rany=1
    dplot=True

    dpival=100

    minlevel=0
    sectors=15

    #input file
    galfile= "galfit.01"

    #sb output file
    sboutput = "sbout.txt"

    #output file
    output = "out.txt"

    # object name to search in NED

    objname="none"
    band="R"

    InDistMod=0
    InMagCor=0
    InScale=1
    InSbDim=0
   

    namefile="none"
    namepng="none.png"
    namesec="none-gal.png"
    namemod="none-mod.png"
    namemul="none-mul.png"
    namesub="none-sub.fits"
    namesig="none-sig.fits"
    namesnr="none-snr.fits"
    namened="none-ned.xml"
    namecheck="none-check.fits"



### class for Galfit parameters
class GalfitParams:

    xc=1        #for sectors_photometry
    yc=1        #for sectors_photometry
    q=1         #for sectors_photometry
    ang=0       #for sectors_photometry
    skylevel=0
    scale=1
    inputimage="galaxy.fits"
    outimage="galaxy-out.fits"
    maskimage="galaxy-mask.fits"
    mgzpt=25
    exptime=1
    tempmask="mask.fits"
    xmin=1
    xmax=2
    ymin=1
    ymax=2

    band="R"

    img = np.array([[1,1],[1,1]])
    model = np.array([[1,1],[1,1]])
    imres = np.array([[1,1],[1,1]])

    mask = np.array([[1,1],[1,1]])
    sigma = np.array([[1,1],[1,1]])
    imsnr = np.array([[1,1],[1,1]])


### class for Galfit components
class GalfitComps:

    # init sub values
    Comps=np.array([])  #remove this?
    N=np.array([])

    NameComp=np.array([])  #0)
    PosX=np.array([])            #1)   
    PosY=np.array([])            #2)   
    Mag=np.array([])             #3)
    Rad=np.array([])             #4)
    Exp=np.array([])             #5)
    Exp2=np.array([])            #6)  for moffat
    Exp3=np.array([])            #7)  for moffat
                                  #8)  There is No 8 in any galfit model
    AxRat=np.array([])           #9)  AxisRatio
    PosAng =np.array([])         #10) position angle
    skip=np.array([])            #z)  skip model
    freepar=np.array([])            # Number of free params

    # computed parameters:
    Rad50=np.array([])
    SerInd=np.array([])
    Rad50kpc=np.array([])
    Rad50sec=np.array([])
    Rad90=np.array([])
    AbsMag=np.array([])
    Lum=np.array([])
    Flux=np.array([])
    PerLight=np.array([])
    me=np.array([])
    mme=np.array([])
    kser = np.array([])


##### end of classes

def Help():

    print ("Usage:\n %s [GALFITOutputFile] [-logx] [-q AxisRatio] [-pa PositionAngle] [-comp] [-pix] [-ranx/y Value] [-grid] [-dpi Value] [-noplot] [-phot] " % (sys.argv[0]))
    print ("More options: [-sbout] [-noplot] [-minlevel Value] [-sectors Value] [-object Name] [-filter Name] [-snr] [-help] [-checkimg] [-noned] [-distmod Value] [-magcor Value] [-scalekpc Value][-sbdim Value] ") 

    print ("GALFITOutputFile: GALFIT output file ")
    print ("logx: activates X-axis as logarithm ")
    print ("q: introduce axis ratio ")
    print ("pa: introduce position angle (same as GALFIT) ")
    print ("comp: plots individual components ")
    print ("pix: plot the top x-axis in pixels ")
    print ("ranx: constant that multiplies the range of the x axis or xmin-xmax range")
    print ("rany: constant that multiplies the range of the y axis or ymin-ymax range")
    print ("grid: display a grid in the plot ")
    print ("dpi: dots per inch used for images files ")
    print ("noplot: avoid displaying windows and directly creates images")
    print ("sbout: creates output file containing the surface brightness profiles")


    print ("                OUTPUT               ")
    print ("phot: Compute photometry. Check the created output file")
    print ("the below options are used only if 'phot' is enabled ")    
    print ("      snr: Creates Signal to Noise image ")
    print ("      object: used for 'phot' to search in NED  ")
    print ("      filter: used for 'phot' to indicate band for NED ")
    print ("      noned: avoid to connect to NED")
    print ("any of the following options disabled the connection of NED ")    
    print ("      distmod: Introduce Distance Modulus ")
    print ("      magcor: Introduce Galactic Extinction ")
    print ("      scalekpc: Introduce equivalence of ''/kiloparsec ")
    print ("      sbdim: Introduce surface brightness dimming")
    print ("                ADVANCED               ")
    print ("minlevel: parameter given directly to sectors_photometry.")
    print ("                      It stops when it founds this value")

    print ("sectors: parameter given directly to sectors_photometry. Divide elipse in 'sectors' ")
    print ("                      Check sectors_photometry manual")
    print ("checkimg: save the images used for sectors_photometry in individual components")


    print ("help: This menu ")



    print ("Example:\n %s galfit.01 -logx" % (sys.argv[0]))
    print ("or Example:\n %s galfit.02 -q 0.35 -pa 60 -comp -ranx 2 -out " % (sys.argv[0]))
    print ("or Example:\n %s galfit.02 -q 0.35 -pa 60 -comp -ranx 1-20" % (sys.argv[0]))
    print ("see https://github.com/canorve/GALFITools/blob/master/docs/Ellipse.md  for more examples")



    sys.exit()


    return True



def SectPhot(galpar, params, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for galaxy and model """


    maskb=galpar.mask


    eps=1-galpar.q

    if params.dplot:
        plt.clf()
        print("")

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp=galpar.xc
    xctemp=galpar.yc
    # and angle is different as well:
    angsec=90-galpar.ang
    #    angsec=ang


    ###############################
    #  galaxy:

    sectgalax = sectors_photometry(galpar.img, eps, angsec, xctemp, yctemp, minlevel=minlevel,
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namesec)

    ###################################################
    
    #  model:
    sectmodel = sectors_photometry(galpar.model, eps, angsec, xctemp, yctemp,minlevel=minlevel,
            plot=params.dplot, badpixels=maskb, n_sectors=n_sectors)


    if params.dplot:
        plt.pause(1)  # Allow plot to appear on the screen
        plt.savefig(params.namemod)



    return sectgalax,sectmodel



def SectPhotComp(galpar, params, galcomps, n_sectors=19, minlevel=0):
    """ calls to function sectors_photometry for subcomponents """

    if (params.flagphot) and (not(os.path.isfile(params.namesig))):

        print("running galfit to create sigma image and individual model images...")

        rungal = "galfit -o3 -outsig {}".format(params.galfile)
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        # changing name to subcomponents
        runchg = "mv subcomps.fits {}".format(params.namesub)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(params.namesub)
        assert os.path.isfile(params.namesub), errmsg

        # changing name to sigma image

        runchg = "mv sigma.fits {}".format(params.namesig)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(params.namesig)
        assert os.path.isfile(params.namesig), errmsg

    else: 

        print("running galfit to create individual model images...")

        rungal = "galfit -o3 {}".format(params.galfile)
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        # changing name to subcomponents
        runchg = "mv subcomps.fits {}".format(params.namesub)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(params.namesub)
        assert os.path.isfile(params.namesub), errmsg

        if (os.path.isfile(params.namesig) and (params.flagphot)):
            print("using existing sigma image")

    ##

    hdu = fits.open(params.namesub)

    subimgs=[]

    cnt=0  # image =0 do not count
    while(cnt<len(galcomps.Comps)):
        if galcomps.Comps[cnt] == True:
            img = hdu[cnt+2].data.astype(float)
            subimgs.append(img)
        cnt=cnt+1
    hdu.close()


    maskb=galpar.mask


    eps=1-galpar.q

    ############
    # I have to switch x and y values because they are different axes for
    # numpy:
    yctemp=galpar.xc
    xctemp=galpar.yc
    # and angle is different as well:
    angsec=90-galpar.ang

    epsmul=eps
    angsecmul=angsec
    #print("eps,angle mul ",epsmul,angsecmul)
    #############
    sectcomps=[]
    #sectmulcomps=[]
    n=0

    while(n<len(galcomps.N)):

        subim=subimgs[n]

        eps=1-galcomps.AxRat[n]
        angsec=90-galcomps.PosAng[n]


        if params.flagcheck:
            scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,plot=1, badpixels=maskb, n_sectors=n_sectors)
            plt.savefig("Comp"+str(n)+".png")
        else:
            scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,plot=0, badpixels=maskb, n_sectors=n_sectors)


        #scmpmul = sectors_photometry(subim, epsmul, angsecmul, xctemp, yctemp,minlevel=minlevel,plot=0, badpixels=maskb, n_sectors=n_sectors)
        #plt.savefig("Cmul"+str(n)+".png")

        sectcomps.append(scmp)

        #sectmulcomps.append(scmpmul)

        n=n+1


    return sectcomps


def EllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps,n_sectors=19, minlevel=0):


    xradm = []
    ysbm = []
    ysberrm = []


    # galax 

    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]
    mgecount=sectgalax.counts[stidxg]
    mgeangle=sectgalax.angle[stidxg]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

    #changing to arc sec
    aellarcg=aellabg*galpar.scale


    ####
    # formula according to cappellary mge manual:
    mgesbg= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ####

    stidxq = np.argsort(aellarcg)


    xarcg = aellarcg[stidxq]
    ymgeg = mgesbg[stidxq]

    #############  Function to order SB along X-axis for galaxy

    xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)

    ################

    # model
    stidxm = np.argsort(sectmodel.radius)

    mgerad=sectmodel.radius[stidxm]
    mgecount=sectmodel.counts[stidxm]
    mgeangle=sectmodel.angle[stidxm]
    mgeanrad=np.deg2rad(mgeangle)

    ab=galpar.q

    aellabm= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarcm=aellabm*galpar.scale

    # formula according to cappellary mge manual
    mgesbm= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    ##

    stidxq = np.argsort(aellarcm)

    xarcm = aellarcm[stidxq]
    ymgem = mgesbm[stidxq]




    ######  Function to order SB along X-axis for model

    xradm, ysbm, ysberrm    = FindSB(xarcm, ymgem, n_sectors)

    ################ Plotting

    limx,limy,axsec=PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,galpar.scale)


    ### surface brightness output file

    if params.flagsbout == True: 

        #print to file    
        PrintEllFilesGax(params,galpar,xradq,ysbq,ysberrq,xradm,ysbm,ysberrm)


    #### Creating Subcomponents images with Galfit


    if params.flagcomp:


        xradq,ysbq,n=SubComp(params, galpar, galcomps, sectcomps, axsec, n_sectors=n_sectors)



    axsec.legend(loc=1)

    return limx,limy


def PrintEllFilesGax(params,galpar,xradq,ysbq,ysberrq,xradm,ysbm,ysberrm):
    "print surface brightness of galaxy and model to file"

    # output for galaxy
    OUTFH = open (params.sboutput+".gal.txt","w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#            Galaxy                                   \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB        SBerr       \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec) (error)   \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradq))):
        lineout= "{0:.3f} {1:.3f} {2:.3f} \n".format(xradq[idx],ysbq[idx],ysberrq[idx])
        OUTFH.write(lineout)

    OUTFH.close()

    # output for model 
    OUTFH = open (params.sboutput+".mod.txt","w")

    lineout= "#        sectors_photometry used with q={} and pa={} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#           Surface Brightness   Model                \n"
    OUTFH.write(lineout)

    lineout= "#    rad       SB        SBerr \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec) (error)  \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradm))):
            
        lineout= "{0:.3f} {1:.3f} {2:.3f} \n".format(xradm[idx],ysbm[idx],ysberrm[idx])

        OUTFH.write(lineout)

    OUTFH.close()








def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,params,scale):
    """  Produces final best-fitting plot  """

    # subplot for arc sec axis
    fig, axsec = plt.subplots()


    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)


    minrad = np.min(xradq)
    if params.flagranx[1] == False:
        maxrad = np.max(xradq) * params.ranx
    else:
        maxrad = np.max(xradq)

    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt/mincnt)**np.array([-0.05, +1.05])

    if params.flagrany[1] == False:
        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran2 #inverted axis
        yran[1] = yran1


    axsec.set_xlabel("radius ('')")
    axsec.set_ylabel("mag/''")

    axsec.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy",linewidth=2)
    axsec.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model",linewidth=2)
    if params.flagranx[1] == False:
        axsec.set_xlim(xran)
    else:
        axsec.set_xlim(xmin,xmax)


    if params.flagrany[1] == False:
        axsec.set_ylim(yran)
    else:
        axsec.set_ylim(ymax,ymin) #inverted


    if params.flaglogx == True:

        axsec.set_xscale("log")

        locmaj = LogLocator(base=10,numticks=12)
        axsec.xaxis.set_major_locator(locmaj)

        locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
        axsec.xaxis.set_minor_locator(locmin)
        axsec.xaxis.set_minor_formatter(NullFormatter())

    else:
        axsec.xaxis.set_minor_locator(AutoMinorLocator())
        axsec.xaxis.set_major_locator(AutoLocator())


    axsec.tick_params(which='both', width=2)
    axsec.tick_params(which='major', length=7)
    axsec.tick_params(which='minor', length=4, color='r')

    axsec.yaxis.set_minor_locator(AutoMinorLocator())
    axsec.yaxis.set_major_locator(AutoLocator())

    if params.flagpix == True:

        axpix = axsec.twiny()

        axpix.set_xlabel("(pixels)")
        x1, x2 = axsec.get_xlim()

        axpix.set_xlim(x1/scale, x2/scale)

        axpix.figure.canvas.draw()


        if params.flaglogx == True:
            axpix.set_xscale("log")
            locmaj = LogLocator(base=10,numticks=12)
            axpix.xaxis.set_major_locator(locmaj)

            locmin = LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
            axpix.xaxis.set_minor_locator(locmin)
            axpix.xaxis.set_minor_formatter(NullFormatter())
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())
            axpix.xaxis.set_major_locator(AutoLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        axret=axsec
    else:
        axret=axsec

    if params.flagrid == True:
        # Customize the major grid
        axsec.grid(which='major', linestyle='-', linewidth='0.7', color='black')
        # Customize the minor grid
        axsec.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

    return xran,yran,axret


def SubComp(params, galpar, galcomps, sectcomps, axsec, n_sectors=19):

    N=len(galcomps.N)

    #color value
    values = range(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    ####################
    ab=galpar.q
    n=0

    while(n<N):

        
        namec=galcomps.NameComp[n]

        scmp = sectcomps[n] 


        ###################################################

        stidxg = np.argsort(scmp.radius)

        mgerad=scmp.radius[stidxg]
        mgecount=scmp.counts[stidxg]
        mgeangle=scmp.angle[stidxg]
        mgeanrad=np.deg2rad(mgeangle)


        aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

        #converting to arcsec

        aellarcg=aellabg*galpar.scale


        # formula according to cappellary mge manual

        mgesbg= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1


        stidxq = np.argsort(aellarcg)

        xarcg = aellarcg[stidxq]
        ymgeg = mgesbg[stidxq]

        xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)

        colorVal = scalarMap.to_rgba(values[n])
        PlotSub(xradq,ysbq,n,axsec,namec,colorVal)


        if params.flagsbout == True:
            ncomp=n+1
            ncomp=str(ncomp)

            PrintEllFilesComps(params,galpar,namec,ncomp,xradq,ysbq,ysberrq)


        n=n+1


    return  xradq,ysbq,n

def PrintEllFilesComps(params,galpar,namecomp,ncomp,xradq,ysbq,ysberrq):
    "Print surface brigthness of components to file "
    #subcomponent model 

    OUTFH = open (params.sboutput+".comp-"+ncomp+".txt","w")

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model  {}   component {}            \n".format(namecomp,ncomp)
    OUTFH.write(lineout)

    lineout= "#     rad      SB   SBerror  \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) (error) \n"
    OUTFH.write(lineout)

    for idx, item in reversed(list(enumerate(xradq))):

        lineout= "{0:.3f} {1:.3f} {2:.3f}  \n".format(xradq[idx],ysbq[idx],ysberrq[idx])

        OUTFH.write(lineout)

    OUTFH.close()




def PlotSub(xradq,ysbq,nsub,axsec,namec,colorval):
    """
    Produces subcomponent plot

    """

    substr=namec+" "+np.str(nsub+1)

    # axsec.plot(xradq, ysbq,'--',color='skyblue',linewidth=4,markersize=0.7,label=substr)
    axsec.plot(xradq, ysbq,'--',color=colorval,linewidth=1.5,markersize=0.7,label=substr)



def MulEllipSectors(params, galpar, galcomps, sectgalax, sectmodel, sectcomps):

  
    fignum=1


    eps=1-galpar.q

    if params.flagranx[1] == True:
        (xmin,xmax)=params.ranx.split("-")
        xmin=np.float(xmin)
        xmax=np.float(xmax)

    if params.flagrany[1] == True:
        (ymin,ymax)=params.rany.split("-")
        ymin=np.float(ymin)
        ymax=np.float(ymax)


    yctemp=galpar.xc
    xctemp=galpar.yc


    angsec=90-galpar.ang

    ######################

    sg = sectgalax

    sm = sectmodel
    ###################################################

    stidx = np.argsort(sg.radius)

    #   galaxy
    mgerad=sg.radius[stidx]

    mgecount=sg.counts[stidx]
    mgeangle=sg.angle[stidx]
    mgeanrad=np.deg2rad(mgeangle)


    # model

    stidx = np.argsort(sm.radius)

    mgemodrad=sm.radius[stidx]

    mgemodcount=sm.counts[stidx]
    mgemodangle=sm.angle[stidx]
    mgemodanrad=np.deg2rad(mgemodangle)


    # converting to pixels

    mgerad=mgerad*galpar.scale
    mgemodrad=mgemodrad*galpar.scale


    # formula according to cappellary mge manual
    # galaxy:
    mgesb= galpar.mgzpt - 2.5*np.log10(mgecount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
    # Model:
    mgemodsb= galpar.mgzpt - 2.5*np.log10(mgemodcount/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1


    if params.flagcomp:

        wtemp=[]
        mgesbsub=[]
        mgeradsub=[]
        mgeanglesub=[]
        sectorsub=[]

        ###############################
        
        ab=galpar.q
        ni=0
        while(ni<len(galcomps.N)):

            subcmp = sectcomps[ni]


            subidx = np.argsort(subcmp.radius)

            temprad=subcmp.radius[subidx]

            #converting to arcsec
            temprad=temprad*galpar.scale

            mgecountsub=subcmp.counts[subidx]

            tempangle=subcmp.angle[subidx]
            mgeanradsub=np.deg2rad(tempangle)


            # formula according to cappellary mge manual
            tempmge= galpar.mgzpt - 2.5*np.log10(mgecountsub/galpar.exptime) + 2.5*np.log10(galpar.scale**2) + 0.1
            

            tempsectorsub = np.unique(tempangle)

            sectorsub.append(tempsectorsub)
            mgesbsub.append(tempmge)
            mgeradsub.append(temprad)
            mgeanglesub.append(tempangle)

            ni+=1


    minrad = np.min(mgerad)

    if params.flagranx[1] == False:
        maxrad = np.max(mgerad) * params.ranx
    else:
        maxrad = np.max(mgerad)

    minsb = np.min(mgesb)
    maxsb = np.max(mgesb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = minsb * (maxsb/minsb)**np.array([+1.05,-0.05])


    if params.flagrany[1] == False:

        yran1=yran[0]
        yran2=yran[1]

        lyran= yran2 - yran1

        yranmid= yran1 + lyran/2

        lyran=lyran*params.rany

        yran1 = yranmid - lyran/2
        yran2 = yranmid + lyran/2

        yran[0] = yran1
        yran[1] = yran2


    sectors = np.unique(mgeangle)
    n = sectors.size
    dn = int(round(n/6.))
    nrows = (n-1)//dn + 1 # integer division

    plt.clf()

    fig, axsec = plt.subplots(nrows, 2, sharex=True, sharey='col', num=fignum)
    fig.subplots_adjust(hspace=0.01)


    if params.flagpix:
        axpix = axsec[0,0].twiny()
        axpix2 = axsec[0,1].twiny()

    fig.text(0.04, 0.5, 'Surface brightness', va='center', rotation='vertical')
    fig.text(0.96, 0.5, 'error (%)', va='center', rotation='vertical')

    axsec[-1, 0].set_xlabel("radius ('')")
    axsec[-1, 1].set_xlabel("radius ('')")

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        if params.flagpix:
            axpix.set_xscale("log")
            axpix2.set_xscale("log")

    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())

    axsec[-1, 0].tick_params(which='both', width=2)
    axsec[-1, 0].tick_params(which='major', length=7)
    axsec[-1, 0].tick_params(which='minor', length=4, color='r')

    if params.flaglogx == True:
        axsec[-1, 0].xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
    else:
        axsec[-1, 0].xaxis.set_minor_locator(AutoMinorLocator())
    axsec[-1, 1].tick_params(which='both', width=2)
    axsec[-1, 1].tick_params(which='major', length=7)
    axsec[-1, 1].tick_params(which='minor', length=4, color='r')


    #    row = 7 # old values
    row = nrows -1
    for j in range(0, n, dn):
        w = np.nonzero(mgeangle == sectors[j])[0]

        w = w[np.argsort(mgerad[w])]
        r = mgerad[w]

        wmod = np.nonzero(mgemodangle == sectors[j])[0]
        wmod = wmod[np.argsort(mgemodrad[wmod])]

        if (len(mgemodrad) < len(mgerad)):
            r2 = mgemodrad[wmod]
        else:
            wmod=w
            r2 = mgemodrad[wmod]


        txtang= sectors[j]
        txt = "$%.f^\circ$" % txtang

        if params.flagranx[1] == False:
            axsec[row, 0].set_xlim(xran)
        else:
            axsec[row, 0].set_xlim(xmin,xmax)

        if params.flagrany[1] == False:
            axsec[row, 0].set_ylim(yran)
        else:
            axsec[row, 0].set_ylim(ymax,ymin) #inverted

        if params.flaglogx == False:
            axsec[row, 0].plot(r, mgesb[w], 'C3o')

            axsec[row, 0].plot(r2, mgemodsb[wmod], 'C0-', linewidth=2)

        else:
            axsec[row, 0].semilogx(r, mgesb[w], 'C3o')

            axsec[row, 0].semilogx(r2, mgemodsb[wmod], 'C0-', linewidth=2)

        if params.flagsbout == True: 

            rtxtang=np.int(np.round(txtang)) 

            PrintFilesGax(params,galpar,rtxtang,r,mgesb,w,r2,mgemodsb,wmod)


        if params.flagrid == True:
            # Customize the major grid
            axsec[row,0].grid(which='major', linestyle='-', linewidth='0.7', color='black')
            # Customize the minor grid
            axsec[row,0].grid(which='minor', linestyle=':', linewidth='0.5', color='black')

            #  axsec[row,0].grid(True)

        if params.flagcomp == True:
            ii=0
                #color value
            values = range(len(galcomps.N))
            jet = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

            while(ii<len(galcomps.N)):

                #wtemp = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
        
######################## Patch for angle :############  
                alpha = sectors[j]
                angsec2= 90-galcomps.PosAng[ii]
                if angsec < 0:
                    angsec = 360 + angsec
                if angsec2 < 0:
                    angsec2 = 360 + angsec2

                alpha2 = alpha  + angsec - angsec2 

                if alpha2 < 0:
                    alpha2 = 360 + alpha2

                if alpha2 > 90 and alpha2 <=270:
                    alpha2 = np.abs(180-alpha2)

                if alpha2 > 270 and alpha2 <=360:
                    alpha2 = np.abs(360-alpha2)


                # search for the nearest angle for subcomponent:
                jj=(np.abs(sectorsub[ii]-alpha2)).argmin()  

                diffangle =  sectorsub[ii][jj] - alpha2

                #print("Check C{}: Axrat: {:.3f}, alpha: {:.3f} angsec: {:.3f} ; theta2: {:.3f} sector {:.3f}; alpha2 {:.3f} ".format(ii,galcomps.AxRat[ii],alpha,angsec,90-galcomps.PosAng[ii],sectorsub[ii][jj],alpha2))



###############################################
                
                #wtemp = np.nonzero(mgeanglesub[ii] == sectorsub[ii][j])[0]
                wtemp = np.nonzero(mgeanglesub[ii] == sectorsub[ii][jj])[0]
                wtemp = wtemp[np.argsort(mgeradsub[ii][wtemp])]

      
                rtemp = mgeradsub[ii][wtemp]

                colorval = scalarMap.to_rgba(values[ii])
                if params.flaglogx == False:
                #    axsec[row, 0].plot(rtemp, mgesbsub[ii][wtemp],'--',color='skyblue', linewidth=2)
                    axsec[row, 0].plot(rtemp, mgesbsub[ii][wtemp],'--',color=colorval, linewidth=1.5)
                else:
                    axsec[row, 0].semilogx(rtemp, mgesbsub[ii][wtemp], '--',color=colorval, linewidth=1.5)

                #introduce output 
                if params.flagsbout == True:
                    ncomp=ii+1
                    ncomp=str(ncomp)

                    PrintFilesComps(params,galpar,galcomps,rtxtang,ncomp,diffangle,rtemp,mgesbsub,ii,wtemp)

                ii+=1

        axsec[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=axsec[row, 0].transAxes)

        if (len(mgemodrad) < len(mgerad)):
            sberr=1-mgemodsb[wmod]/mgesb[wmod]
            axsec[row, 1].plot(r2, sberr*100, 'C0o')
        else:
            sberr=1-mgemodsb[w]/mgesb[w]
            axsec[row, 1].plot(r, sberr*100, 'C0o')


        axsec[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        axsec[row, 1].yaxis.tick_right()
        axsec[row, 1].yaxis.set_label_position("right")
        axsec[row, 1].set_ylim([-19.5, 20])
        # axsec[row, 1].set_ylim([-20, 20])


        if params.flagranx[1] == False:
            axsec[row, 1].set_xlim(xran)
        else:
            axsec[row, 1].set_xlim(xmin,xmax)


        axsec[row, 0].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 0].tick_params(which='both', width=2)
        axsec[row, 0].tick_params(which='major', length=7)
        axsec[row, 0].tick_params(which='minor', length=4, color='r')

        axsec[row, 1].yaxis.set_minor_locator(AutoMinorLocator())
        axsec[row, 1].tick_params(which='both', width=2)
        axsec[row, 1].tick_params(which='major', length=7)
        axsec[row, 1].tick_params(which='minor', length=4, color='r')


        row -= 1

    if params.flagpix == True:
        axpix.set_xlabel("(pixels)")
     
        #x1, x2 = axsec[7,0].get_xlim() ## buggy for some data have to change it for code below:
        
        if params.flagranx[1] == False:
            x1= xran[0]
            x2= xran[1]
        else:
            x1=xmin
            x2=xmax

        axpix.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix.figure.canvas.draw()

        axpix2.set_xlabel("(pixels)")
        axpix2.set_xlim(x1/galpar.scale, x2/galpar.scale)
        axpix2.figure.canvas.draw()

        ##
        if params.flaglogx == True:
            axpix.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix.xaxis.set_minor_locator(AutoMinorLocator())

        axpix.tick_params(which='both', width=2)
        axpix.tick_params(which='major', length=7)
        axpix.tick_params(which='minor', length=4, color='r')

        if params.flaglogx == True:
            axpix2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=15))
        else:
            axpix2.xaxis.set_minor_locator(AutoMinorLocator())
        axpix2.tick_params(which='both', width=2)
        axpix2.tick_params(which='major', length=7)
        axpix2.tick_params(which='minor', length=4, color='r')


def PrintFilesGax(params,galpar,rtxtang,r,mgesb,w,r2,mgemodsb,wmod):
    "Print surface parameters of galaxy and model to outfile "

    # galaxy
    OUTFH = open (params.sboutput+"-"+str(rtxtang)+".gal.txt","w")

    lineout= "# Values along radius with ang = {} from major axis \n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel]\n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#            Galaxy                                   \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB              \n"
    OUTFH.write(lineout)

    lineout= "# (arcsec) (mag/arcsec)    \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(r):

        lineout= "{0:.3f} {1:.3f} \n".format(r[idx],mgesb[w][idx])

        OUTFH.write(lineout)

    OUTFH.close()

    #model
    OUTFH = open (params.sboutput+"-"+str(rtxtang)+".mod.txt","w")

    lineout= "# Values along radius with ang = {} from major axis \n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model                \n"
    OUTFH.write(lineout)

    lineout= "#     rad      SB     \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(r2):

        lineout= "{0:.3f} {1:.3f}  \n".format(r2[idx],mgemodsb[wmod][idx])

        OUTFH.write(lineout)

    OUTFH.close()


def PrintFilesComps(params,galpar,galcomps,rtxtang,ncomp,diffangle,rtemp,mgesbsub,ii,wtemp):

    #subcomponent model 
    OUTFH = open (params.sboutput+"-"+str(rtxtang)+".sub-"+ncomp+".txt","w")

    lineout= "# Values along radius with ang = {} from major axis \n".format(rtxtang)
    OUTFH.write(lineout)

    lineout= "# sectors_photometry used with  q = {} and pa = {} (same as GALFIT) \n".format(galcomps.AxRat[ii],90-galcomps.PosAng[ii])
    OUTFH.write(lineout)

    lineout= "# In the multiplot, there is a difference of {:.3f} for the one indicated in the top right corner.\n".format(diffangle)
    OUTFH.write(lineout)

    lineout= "# The above is due to differences in the sectors_photometry for the galaxy and individual components.\n"
    OUTFH.write(lineout)

    lineout= "#  OutImage = {}  magzpt = {}  exptime = {}  plate scale = {} [arcsec per pixel] \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTFH.write(lineout)

    lineout= "#  xc = {}  yc = {}  sky = {}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTFH.write(lineout)

    lineout= "#        Model   component {}            \n".format(ncomp)
    OUTFH.write(lineout)

    lineout= "#     rad      SB     \n"
    OUTFH.write(lineout)

    lineout= "#   (arcsec) (mag/arcsec) \n"
    OUTFH.write(lineout)

    for idx, item in enumerate(rtemp):

        lineout= "{0:.3f} {1:.3f}  \n".format(rtemp[idx], mgesbsub[ii][wtemp][idx])

        OUTFH.write(lineout)

    OUTFH.close()




def FindSB(xarcq, ymgeq, numsectors):
    # the xarcq array must be ordered
    # use mag instead of counts

    xradq=[]
    ysbq=[]
    ysberrq=[]
    xradq=np.array(xradq)
    ysbq=np.array(ysbq)
    ysberrq=np.array(ysberrq)

    numsave=0
    tot=xarcq.size
    count=0
    for i in range(tot,0,-1):

        lima=i-numsectors
        limb=i

        valstd=np.std(xarcq[lima:limb])
        if valstd < 0.1:
            numsave=count
            break
        count=count+1
    init=numsave%numsectors

    n=init

    num=np.int((xarcq.size-init)/numsectors)
    n=xarcq.size-init
    for i in range(num,0,-1):

        lima=n-numsectors
        limb=n

        xradq=np.append(xradq,np.mean(xarcq[lima:limb]))
        ysbq=np.append(ysbq,np.mean(ymgeq[lima:limb]))
        ysberrq=np.append(ysberrq,np.std(ymgeq[lima:limb]))

        n=n-numsectors

    return xradq, ysbq, ysberrq


def ReadGALFITout(inputf,galpar):

    flagfirst = True

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
        if tmp[0] == "A)":     # input image
            galpar.inputimage=tmp[1]

        if tmp[0] == "B)":     # out image
            galpar.outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            galpar.maskimage=tmp[1]

        if tmp[0] == "H)":     # region fit box
            galpar.xmin=int(tmp[1])
            galpar.xmax=int(tmp[2])
            galpar.ymin=int(tmp[3])
            galpar.ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            galpar.mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            galpar.scale=float(tmp[1])


        # first component
        if tmp[0] == "0)" and flagfirst == True:     # plate scale

            flagfirst=False

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    galpar.xc=float(tmp[1])
                    galpar.yc=float(tmp[2])

                if tmp[0] == "9)":    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)": # position angle
                    galpar.ang=float(tmp[1])

        # sersic component
        if tmp[0] == "0)" and tmp[1] == "sersic":     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcser=float(tmp[1])
                    ycser=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcser)**2+(galpar.yc-ycser)**2)



                if tmp[0] == "9)" and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if tmp[0] == "10)" and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])

        # second component exponential model
        if tmp[0] == "0)" and tmp[1] == "expdisk" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcexp=float(tmp[1])
                    ycexp=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcexp)**2+(galpar.yc-ycexp)**2)


                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])


        if tmp[0] == "0)" and tmp[1] == "gaussian":

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgauss=float(tmp[1])
                    ycgauss=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((galpar.xc-xcgauss)**2+(galpar.yc-ycgauss)**2)

                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    galpar.q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
                    galpar.ang=float(tmp[1])



        if tmp[0] == "0)" and tmp[1] == "sky" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":    # axis ratio
                    galpar.skylevel=float(tmp[1])

        index += 1

    GalfitFile.close()


    errmsg="file {} does not exist".format(galpar.inputimage)
    assert os.path.isfile(galpar.inputimage), errmsg

    galpar.exptime=GetExpTime(galpar.inputimage)


    #errmsg="xc and yc are unknown "
    #assert ("xc" in locals()) and ("yc" in locals())  , errmsg

    print("center is at xc, yc = ",galpar.xc,galpar.yc)

    # correcting coordinates
    galpar.xc=galpar.xc-galpar.xmin+1
    galpar.yc=galpar.yc-galpar.ymin+1


    if os.path.isfile(galpar.maskimage):
        mime=mimetypes.guess_type(galpar.maskimage)

        flagbm = not(mime[0] == "text/plain")

        errmsg="Sorry the mask file: {}  must be binary, not ASCII ".format(maskimage)
        assert flagbm is True, errmsg


        galpar.tempmask="tempmask.fits"
        GetFits(galpar.maskimage, galpar.tempmask, galpar.xmin, galpar.xmax, galpar.ymin, galpar.ymax)

    else:
        errmsg="Mask file does not exist"
        print(errmsg)
        galpar.tempmask=None

    #return xc,yc,q,pa,skylevel,scale,outimage,mgzpt,exptime,mask


def ReadNComp(inputf,X,Y,galcomps):
    ## search and count model components

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    n=0

    distmin=3 # minimum distance for among component centers

    while index < len(lines):

        line = lines[index]
        (tmp) = line.split()

        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

            X=X+xmin-1
            Y=Y+ymin-1

        #init values
        Comps=False
        N=0
        NameComp="none"
        PosX=0
        PosY=0
        Mag=99
        Rad=0
        Exp=0
        Exp2=0
        Exp3=0
        AxRat=1
        PosAng=0
        skip=1
        flagcomp=False  # detect components
        freepar=0
        if tmp[0] == "0)":

            namec=tmp[1] 
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xc=float(tmp[1])
                    yc=float(tmp[2])

                    dist = np.sqrt((xc-X)**2+(yc-Y)**2)
                    if (dist < distmin and namec != "sky"):
                        n=n+1
                        PosX=xc
                        PosY=yc
                        Comps=True
                        NameComp=namec
                        N= n 
                        flagcomp=True
                        par1=int(tmp[3])
                        par2=int(tmp[4])
                        # count the number of free params
                        freepar=freepar+par1+par2  
                    else:
                        Comps=False

                if tmp[0] == "3)" and flagcomp == True:    # axis ratio
                    Mag=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "4)" and flagcomp == True:    # axis ratio
                    Rad=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "5)" and flagcomp == True:    # axis ratio
                    Exp=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "6)" and flagcomp == True:    # axis ratio
                    Exp2=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "7)" and flagcomp == True:    # axis ratio
                    Exp3=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "9)" and flagcomp == True:    # axis ratio
                    AxRat=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "10)" and flagcomp == True:    # axis ratio
                    PosAng=float(tmp[1])
                    par=int(tmp[2])
                    freepar+=par
                if tmp[0] == "z)" and flagcomp == True:    # axis ratio
                    skip=int(tmp[1])
            if (flagcomp == True):

                galcomps.PosX=np.append(galcomps.PosX,PosX)
                galcomps.PosY=np.append(galcomps.PosY,PosY)
                galcomps.Comps=np.append(galcomps.Comps,Comps)
                galcomps.NameComp=np.append(galcomps.NameComp,NameComp)
                galcomps.N=np.append(galcomps.N, N)
                
                galcomps.Mag=np.append(galcomps.Mag,Mag)
                galcomps.Rad=np.append(galcomps.Rad,Rad)
                galcomps.Exp=np.append(galcomps.Exp,Exp)
                galcomps.Exp2=np.append(galcomps.Exp2,Exp2)
                galcomps.Exp3=np.append(galcomps.Exp3,Exp3)
                galcomps.AxRat=np.append(galcomps.AxRat,AxRat)
                galcomps.PosAng=np.append(galcomps.PosAng,PosAng)
                galcomps.skip=np.append(galcomps.skip,skip)
                galcomps.freepar=np.append(galcomps.freepar,freepar)

        index += 1

    GalfitFile.close()

    tot=galcomps.Comps.size

    # computed parameters:
    galcomps.Rad50=np.array([0.0]*tot)
    galcomps.SerInd=np.array([0.0]*tot)
    galcomps.Rad50kpc=np.array([0.0]*tot)
    galcomps.Rad50sec=np.array([0.0]*tot)
    galcomps.Rad90=np.array([0.0]*tot)
    galcomps.AbsMag=np.array([99.0]*tot)
    galcomps.Lum=np.array([0.0]*tot)
    galcomps.Flux=np.array([0.0]*tot)
    galcomps.PerLight=np.array([0.0]*tot)
    galcomps.me=np.array([99.0]*tot)
    galcomps.mme=np.array([99.0]*tot)
    galcomps.kser = np.array([0.0]*tot)


    galcomps.N=galcomps.N.astype(int)

    return True



def GetExpTime(Image):
    "Get exposition time from the image"

    hdu = fits.open(Image)
    exptime = hdu[0].header.get("EXPTIME",1) # return 1 if not found

    hdu.close()
    return exptime


def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"

    if os.path.isfile(Imageout):
        print("{} deleted; a new one is created \n".format(Imageout))
        runcmd = "rm {}".format(Imageout)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)

    hdu = fits.open(Image)
    dat = hdu[0].data[ylo - 1:yhi, xlo - 1:xhi]
    hdu[0].data = dat
    hdu.writeto(Imageout, overwrite=True)
    hdu.close()


def OutPhot(params, galpar, galcomps, sectgalax, sectmodel, sectcomps):
    """ Output photometry for further analysis """


    # masks to identify components: 

    maskmag=(galcomps.NameComp != "ferrer") & (galcomps.NameComp != "nuker") & (galcomps.NameComp != "edgedisk") & (galcomps.NameComp != "king") 
    maskgalax = (galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | (galcomps.NameComp == "expdisk")  | (galcomps.NameComp == "gaussian") 


    maskdisk = (galcomps.NameComp == "expdisk") # ignore this at the moment: #| (galcomps.NameComp == "edgedisk") 
    maskbulge = (galcomps.NameComp == "sersic") | (galcomps.NameComp == "devauc") | (galcomps.NameComp == "moffat") | (galcomps.NameComp == "ferrer") | (galcomps.NameComp == "king") | (galcomps.NameComp == "gaussian") | (galcomps.NameComp == "psf")
    masksersic = (galcomps.NameComp == "sersic") 
    maskexp= (galcomps.NameComp == "expdisk")
    maskdevauc=(galcomps.NameComp == "devauc")
    maskgauss=(galcomps.NameComp == "gaussian")


    if maskmag.any():

        galcomps.Flux[maskmag]=10**((galpar.mgzpt - galcomps.Mag[maskmag])/2.5)
        totFlux=galcomps.Flux[maskmag].sum()
        totMag=-2.5*np.log10(totFlux) + galpar.mgzpt
        #print("total magnitud = ",totMag)
        #print("total Flux = ",totFlux)
    #else:
    #    print("Total magnitud can not be computed with the actual galfit functions")



    if maskdisk.any():
        BulgeFlux = 10**((galpar.mgzpt -  galcomps.Mag[maskbulge])/2.5)
        DiskFlux  = 10**((galpar.mgzpt -  galcomps.Mag[maskdisk])/2.5)
        totBulgeF = BulgeFlux.sum()
        totDiskF =  DiskFlux.sum()
        BulgeToTotal= totBulgeF / (totBulgeF + totDiskF)
        #print("BulgeToTotal = ",BulgeToTotal)
    else:
        BulgeToTotal = 1
        #print("BulgeToTotal = 1")



    ################## Computing component variables ###########

    # sersic index
    galcomps.SerInd[masksersic] = galcomps.Exp[masksersic] 
    galcomps.SerInd[maskdevauc] = 4
    galcomps.SerInd[maskexp]    = 1
    galcomps.SerInd[maskgauss]  = 0.5


    # effective radius
    galcomps.Rad50[masksersic] = galcomps.Rad[masksersic] 
    galcomps.Rad50[maskdevauc] = galcomps.Rad[maskdevauc] 

    galcomps.Rad50[maskexp] = 1.678*galcomps.Rad[maskexp] 
    galcomps.Rad50[maskgauss] = 0.5*galcomps.Rad[maskgauss] 

    ################
    # computing Rad 90% light with aproximation taken from my Thesis

    galcomps.Rad50sec[maskgalax] = galcomps.Rad50[maskgalax] * galpar.scale 

    galcomps.Rad90[maskgalax] = galcomps.Rad50[maskgalax] * (1.53 + 0.73 * galcomps.SerInd[maskgalax] + 0.07 * galcomps.SerInd[maskgalax]**2) 

    ######

    galcomps.kser[maskgalax]  = GetKAprox(galcomps.SerInd[maskgalax])

    # computing meanme and me 
    galcomps.mme[maskgalax]  = galcomps.Mag[maskgalax]  + 2.5 * np.log10(2 * np.pi * galcomps.AxRat[maskgalax] * galcomps.Rad50sec[maskgalax]**2 )


    fn = (( galcomps.AxRat[maskgalax] * galcomps.SerInd[maskgalax] * np.exp( galcomps.kser[maskgalax])) / (galcomps.kser[maskgalax] ** (2 * galcomps.SerInd[maskgalax] )) ) * ( np.exp(scipy.special.gammaln(2*galcomps.SerInd[maskgalax])) )

    galcomps.me[maskgalax] = galcomps.mme[maskgalax] +  2.5 * np.log10( fn )



    Num=len(galcomps.Flux[maskmag])

    galcomps.PerLight[maskmag]= (galcomps.Flux[maskmag] / totFlux )*100
    namecomp=galcomps.NameComp[maskmag]
    N=galcomps.N[maskmag]

    N=N.astype(int)

    #n=0
    #while(n<Num):

    #    print("Num, componente, %perlight ",N[n],namecomp[n],galcomps.PerLight[n])

    #    n+=1

   
    #print("Number of free params: ",int(galcomps.freepar.sum()))


    # create sigma image
    if(not(os.path.isfile(params.namesig))):

        print("running galfit to create sigma image...")

        rungal = "galfit -o1 -outsig {} ".format(params.galfile) 
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        runchg = "mv sigma.fits {}".format(params.namesig)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE, stderr=sp.PIPE,
            universal_newlines=True)

        errmsg="file {} does not exist".format(params.namesig)
        assert os.path.isfile(params.namesig), errmsg
    else:
        if not(params.flagcomp):
            print("using existing sigma image ")



    # call to Tidal
    (tidal,objchinu,bump,snr,stdsnr,totsnr,rss,ndof)=Tidal(params, galpar, galcomps, sectgalax, 2)

    
    #print("Tidal = ",tidal)  
    #print("Local Chinu = ",objchinu)  
    #print("Bumpiness = ",bump)  
    #print("mean SNR = ",snr)  
    #print("std SNR = ",stdsnr)  
    #print("total SNR sum over area = ",totsnr)  
    #print("Residual sum of squares =  ",rss)  
    #print("degrees of freedom = ",ndof)  
   
    header = fits.getheader(galpar.inputimage)


    if not(params.flagobj):
        if "OBJECT" in header: 
            params.objname=header["OBJECT"] 
            params.flagobj=True
            print("using object name: {} to search in NED ".format(params.objname))            

        elif "TARGNAME": 
            params.objname=header["TARGNAME"] 
            params.flagobj=True
            print("using object name: {} to search in NED ".format(params.objname))            
        else:
            print("WARNING: name for object not found in header nor it was provided by user") 
            print("Luminosity and absolute magnitude will not be computed") 
    else:
        print("using object name: {} to search in NED ".format(params.objname))            


    if not(params.flagband):
        if "BAND" in header: 
            params.flagband = True
            params.band=header["BAND"] 
        elif "FILTER" in header: 
            params.flagband = True
            params.band=header["FILTER"] 
        elif "FILTNAM" in header: 
            params.flagband = True
            params.band=header["FILTNAM"] 
        elif "FILTNAM1" in header: 
            params.flagband = True
            params.band=header["FILTNAM1"] 
        else:
            print("WARNING: filter not found. I will use default filter: ",params.band) 
            print("use -filter option to change band") 
    else:
        print("using {} band to correct for galactic extinction ".format(params.band)) 

    if (params.flagobj and not(params.flagned)): 
        (GalExt,DistMod,DistMod2,Scalekpc,SbDim)=NED(params, galpar, galcomps)
    else:
        if params.flagned: 
            print("No search in NED. Lum and abs Mag will not be computed")

    if maskmag.any():

        if (params.flagweb and params.flagobj and not(params.flagned)):

            CorMag = totMag - GalExt # corrected by galactic extinction 
            
            AbsMag=CorMag - DistMod # No K correction applied

            AbsMag2=CorMag - DistMod2 # No K correction applied

            # per component: 
            CompCorMag = galcomps.Mag[maskmag] - GalExt # corrected by galactic extinction 
            galcomps.AbsMag[maskmag] = CompCorMag - DistMod # No K correction applied
        

            if params.band in SunMag:
                MSun = SunMag[params.band]

                Lum = 10**((MSun - AbsMag)/2.5)
                # per component 
                galcomps.Lum[maskmag]= 10**((MSun - galcomps.AbsMag[maskmag])/2.5)
            else:
                print("Absolute Magnitude for Band {} was not found. Check filter name ".format(params.band))
                print("Luminosity will not be computed.")

                Lum = 0

            #print("Magnitud Absoluta",AbsMag)
            #print("Magnitud Absoluta using Distance Modulus independen of z ",AbsMag2)
            #print("check references in ",params.namened)


    if maskgalax.any():

        if (params.flagweb and params.flagobj and not(params.flagned)):

            galcomps.Rad50kpc[maskgalax] = galcomps.Rad50[maskgalax] * galpar.scale * Scalekpc

            galcomps.mme[maskgalax] = galcomps.mme[maskgalax] - GalExt - SbDim
            galcomps.me[maskgalax] = galcomps.me[maskgalax] - GalExt - SbDim
 
        else:
            print("mean surface brightness at Re is not corrected for galactic extintion nor surface brightness dimming ")





    ################  INFORMATION CRITERIA ####################


    freepar=int(galcomps.freepar.sum())

    npix = ndof + freepar

    #    ;  AKAIKE INFORMATION CRITERION
    #    ; AIC = χ2 + 2k

    AICrit = objchinu * ndof + 2*freepar
    

    #    ; BAYESIAN INFORMATION CRITERION
    #    ; BIC = χ2 + k * ln(n)

    BICrit = objchinu * ndof + freepar * np.log(npix)

    ## for output only: 
    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]

    aell = mgerad.max() 
    bell = mgerad.max() * galpar.q



    #######  file output:  ######

    print("Creating output photometry file: ",params.output)

    OUTPHOT = open (params.output,"w")

    lineout= "#    Output photometry for {} \n".format(galpar.outimage)
    OUTPHOT.write(lineout)

    lineout= "#  sector_photometry used with q={} and pa={} (same as GALFIT) \n".format(galpar.q,galpar.ang)
    OUTPHOT.write(lineout)

    lineout= "# Total Magnitude and many other variables does not include \n"
    OUTPHOT.write(lineout)

    lineout= "# the following components: ferrer, nuker, edgedisk and king  \n"
    OUTPHOT.write(lineout)


    lineout= "#  OutImage = {}  Mgzpt = {}  exptime = {}  plate scale = {} ''/pix \n".format(galpar.outimage,galpar.mgzpt,galpar.exptime,galpar.scale)
    OUTPHOT.write(lineout)


    lineout= "#  xc = {}  yc = {}  sky = {:.2f}   \n".format(galpar.xc, galpar.yc, galpar.skylevel)
    OUTPHOT.write(lineout)


    lineout = "# most of the photometry is computed within an ellipse defined by sectors_photometry funciton \n"
    OUTPHOT.write(lineout)

    lineout = "# This box size is computed using an ellipse with a = {:.2f} and b = {:.2f} centered at xc,yc \n".format(aell,bell)
    OUTPHOT.write(lineout)

    lineout = "# This ellipse is computed using minlevel = {} (use -minlevel option to change it) \n".format(params.minlevel)
    OUTPHOT.write(lineout)

    lineout = "# All the photometric quantities are computed for band {} (use -filter option to change it) \n".format(params.band)
    OUTPHOT.write(lineout)

    lineout = "# correction constants applied here: \n"
    OUTPHOT.write(lineout)

    lineout = "# Galactic Extinction= {}; Distance Modulus= {}; Distance Modulus independent of z = {} \n".format(GalExt,DistMod,DistMod2)
    OUTPHOT.write(lineout)

    lineout = "# cosmology corrected scale kpc/arcsec =  {}; Surface brightness dimming (mag) = {} \n".format(Scalekpc,SbDim)
    OUTPHOT.write(lineout)

    lineout = "# Magnitudes are not corrected by K-Correction \n"
    OUTPHOT.write(lineout)

    lineout = "# Check references in NED file {} \n\n".format(params.namened)
    OUTPHOT.write(lineout)


    if maskmag.any():
        #Flux=10**((galpar.mgzpt -  galcomps.Mag[maskmag])/2.5)
        #totFlux=Flux.sum()
        #totMag=-2.5*np.log10(totFlux) + galpar.mgzpt

        lineout = "total apparent magnitude (without corrections) = {:.3f}  \n".format(totMag)
        OUTPHOT.write(lineout)

        lineout = "total flux (without corrections) = {:.3f}  \n".format(totFlux)
        OUTPHOT.write(lineout)
    else:
        lineout="Total magnitude can not be computed with the actual galfit functions \n"
        OUTPHOT.write(lineout)


    lineout="Bulge To Total Ratio = {:.3f} \n".format(BulgeToTotal)
    OUTPHOT.write(lineout)

    lineout="Tidal = {:.3f} \n".format(tidal)
    OUTPHOT.write(lineout)

    lineout="Local Chinu = {:.3f} \n".format(objchinu)
    OUTPHOT.write(lineout)

    lineout="Bumpiness = {:.3f} \n".format(bump)
    OUTPHOT.write(lineout)

    lineout = "mean SNR = {:.3f} \n".format(snr)
    OUTPHOT.write(lineout)

    lineout = "std SNR = {:.3f} \n".format(stdsnr)
    OUTPHOT.write(lineout)

    lineout = "total SNR sum = {:.3f} \n".format(totsnr)  
    OUTPHOT.write(lineout)

    lineout = "Residual sum of squares (RSS) = {:.3f}  \n".format(rss)  
    OUTPHOT.write(lineout)

    lineout= "degrees of freedom = {} \n".format(ndof)  
    OUTPHOT.write(lineout)


    lineout="Number of free params = {} \n".format(int(galcomps.freepar.sum()))
    OUTPHOT.write(lineout)
   

    if maskmag.any():

        if (params.flagweb and params.flagobj and not(params.flagned)):

            #CorMag = totMag - GalExt # corrected by galactic extinction 
            
            #AbsMag=CorMag - DistMod # No K correction applied

            #AbsMag2=CorMag - DistMod2 # No K correction applied

            lineout="Absolute Magnitud = {:.3f} \n".format(AbsMag)
            OUTPHOT.write(lineout)

            lineout="Absolute Magnitud using Dist Mod independent of z = {:.3f} \n".format(AbsMag2)
            OUTPHOT.write(lineout)

            lineout="Luminosity = {:.3f} (solar lum) \n".format(Lum)
            OUTPHOT.write(lineout)



    lineout= "Akaike Information Criterion = {:.3f} \n".format(AICrit)  
    OUTPHOT.write(lineout)

    lineout = "Bayesian Information Criterion = {:.3f} \n".format(BICrit)  
    OUTPHOT.write(lineout)


    lineout = "#########################################\n"  
    OUTPHOT.write(lineout)

    lineout = "# Photometric properties per component: #\n"  
    OUTPHOT.write(lineout)

    lineout = "########## Columns: #####################\n"  
    OUTPHOT.write(lineout)

    lineout = "# Number Component %PerLight me(mag) <me>(mag) Flux AbsMag Luminosity(SolarLum) Rad90(pix) Re(kpc)   \n"  
    OUTPHOT.write(lineout)


    for idx, item in enumerate(galcomps.N) :
        lineout= "{0} {1} {2:^6.3f} {3:^6.3f} {4:^6.3f} {5:^10.3f} {6:^6.3f} {7:^17.3f} {8:^10.3f} {9:^10.3f} \n".format(galcomps.N[idx],galcomps.NameComp[idx],galcomps.PerLight[idx],galcomps.me[idx],galcomps.mme[idx],galcomps.Flux[idx],galcomps.AbsMag[idx],galcomps.Lum[idx],galcomps.Rad90[idx],galcomps.Rad50kpc[idx])
        OUTPHOT.write(lineout)

    OUTPHOT.close()




def Tidal(params, galpar, galcomps, sectgalax, rmin):
    "Computes Tidal  values as defined in Tal et al. 2009 AJ. It algo computes Bumpiness"
    "(Blakeslee 2006 ApJ) value defined between rmin and radius"

    # rmin = minimum radius to compute Tidal and Bumpiness (to avoid PSF Mismatch)

    #init values

    pixcount=0
    pixcountsnr=0
    pixcountchi=0
    sumtidal=0
    sumsig=0
    sigma=1
    flux=0
    sumflux=0
    meanflux=0
    resbump=0
    sumres=0
    sumchinu=0
    chinu=0
    varchi=1
    sflux=0
    meanres=0
    varres=0
    numbump=0
    ndof=0

    rss=0

    tidal=0
    objchinu=0
    bump=0
    snr=0

    totsnr=0
    stdsnr=0
    ################

  
    stidxg = np.argsort(sectgalax.radius)

    mgerad=sectgalax.radius[stidxg]
    mgecount=sectgalax.counts[stidxg]
    mgeangle=sectgalax.angle[stidxg]
    mgeanrad=np.deg2rad(mgeangle)



    ab=galpar.q

    ell=1-ab

    aell = mgerad.max() 

    bell = mgerad.max() * ab

    #changing to arc sec
    aellarc=aell*galpar.scale

    #print("major axis, minor axis (pix) ",aell,bell)

    NCol=len(galpar.img[0])
    NRow=len(galpar.img)

    #print("max size ",NCol,NRow)

    #Obj.Angle = Obj.Theta - 90

    #angle computed from y-axis to x-axis  
    Theta=galpar.ang + 90


    (xlo, xhi, ylo, yhi) = GetSize(galpar.xc, galpar.yc, aell, Theta, ab, NCol, NRow)

    #print("size box ",xmin, xmax, ymin, ymax)


    xser=galpar.xc
    yser=galpar.yc


    imgal = galpar.img


    immodel = galpar.model


    imres = galpar.imres


    immask = galpar.mask


    hdu = fits.open(params.namesig)
    header=hdu[0].header  
    galpar.sigma=hdu[0].data
    #hdu.close()

    imsigma = galpar.sigma.astype(float)


    # creates a new image for snr 
    #NCol=len(galpar.img[0])
    #NRow=len(galpar.img)
    #MakeImage(params.namesnr, NCol, NRow):


    header['TypeIMG'] = ('SNR', 'Signal to Noise Ratio image')
    hdu[0].header  =header
    galpar.imsnr=imgal/imsigma
    hdu[0].data = galpar.imsnr

    if params.flagsnr:
        hdu.writeto(params.namesnr, overwrite=True)
        print("SNR image created.. ",params.namesnr)

    hdu.close()


    imell=immask.copy()
    imell.fill(False)


    #    for objchinu, Tidal and SNR
    #    maskm = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates
    #maskm =immask[ylo - 1:yhi, xlo - 1:xhi] == False  # big image coordinates
    maskm =immask == False  # big image coordinates


    #maskm =immask == False  # big image coordinates

    #   mask including rmin for Bumpiness only
    #    maskbum = dat[ylo - 1:yhi, xlo - 1:xhi] == num  # big image coordinates
    #maskbum = immask[ylo - 1:yhi, xlo - 1:xhi] == False  # big image coordinates
    maskbum = immask == False  # big image coordinates


    #############

    theta = 0

    ypos, xpos = np.mgrid[ylo - 1:yhi, xlo - 1:xhi]

    dx = xpos - xser
    dy = ypos - yser

    dist = np.sqrt(dx**2 + dy**2)

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / rmin, np.cos(landa) / rmin)

    xell = xser + rmin * np.cos(angle)
    yell = yser + rmin * np.sin(angle)

    dell = np.sqrt((xell - xser)**2 + (yell - yser)**2)

    mask = dist < dell

    #  correcting for rmin
    maskbum[ylo - 1:yhi, xlo - 1:xhi][mask] = False


    ## identifying area to compute photometry: 
    imell = ExtractEllip(imell, True, xser, yser, aell, Theta, ell, xlo, xhi, ylo, yhi)
    ## xlo, xhi, ylo, yhi makes sure that ellipse will not be outside of this range


    #maskm=maskm*imell
    #maskbum=maskbum*imell

    #maskm=maskm[ylo - 1:yhi, xlo - 1:xhi]*imell
    maskm=maskm*imell

    #maskbum=maskbum[ylo - 1:yhi, xlo - 1:xhi]*imell
    maskbum=maskbum*imell



    hdu = fits.open(galpar.tempmask)
    header=hdu[0].header  

    header['TypeIMG'] = ('Check', 'Check area where output photometry was computed')
    hdu[0].header  =header


    hdu[0].data = (~maskm).astype("int")*1000
    hdu.writeto(params.namecheck, overwrite=True)
    hdu.close()




    if maskbum.any():

    # for Bumpiness

        #galfluxbum  = imgal[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        #modfluxbum  = immodel[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        #sigfluxbum  = imsigma[ylo - 1:yhi, xlo - 1:xhi][maskbum]
        galfluxbum  = imgal[maskbum]
        modfluxbum  = immodel[maskbum]
        sigfluxbum  = imsigma[maskbum]

    ####

    # for Tidal, SNR, and objchinu

    #        sumflux  = np.sum(imgal[ylo - 1:yhi, xlo - 1:xhi][maskm] - sky)
        #sumflux  = np.sum(imgal[ylo - 1:yhi, xlo - 1:xhi][maskm])
        #sumsig   = np.sum(imsigma[ylo - 1:yhi, xlo - 1:xhi][maskm])

        #galflux  = imgal[ylo - 1:yhi, xlo - 1:xhi][maskm]
        #modflux  = immodel[ylo - 1:yhi, xlo - 1:xhi][maskm]

        #sigflux  = imsigma[ylo - 1:yhi, xlo - 1:xhi][maskm]


        sumflux  = np.sum(imgal[maskm])
        sumsig   = np.sum(imsigma[maskm])

        galflux  = imgal[maskm]
        modflux  = immodel[maskm]

        sigflux  = imsigma[maskm]




        resflux = (galflux - modflux)**2

        #rss2=(resflux).sum()

        #print("Residual sum squares ",rss2)
        
        #  local chinu

        varchi = sigflux**2
        chinu  = np.sum(resflux/varchi)
        

        #pixcountchi = np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskm])
        pixcountchi = np.size(immodel[maskm])


        if(pixcountchi > 11):

            ndof=pixcountchi - int(galcomps.freepar.sum())
            objchinu= chinu / ndof
        else:
            objchinu=-1
            ndof=-1


        # snr
        #if(np.size(galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm]) > 0):

         #   totsnr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].sum()

          #  snr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].mean()
           # stdsnr=galpar.imsnr[ylo - 1:yhi, xlo - 1:xhi][maskm].std()

        if(np.size(galpar.imsnr[maskm]) > 0):

            totsnr=galpar.imsnr[maskm].sum()

            snr=galpar.imsnr[maskm].mean()
            stdsnr=galpar.imsnr[maskm].std()


        else:
            print("I can't compute SNR")
            snr=-1
            stdsnr=0
            totsnr=0
        # Tidal parameter

        tgal = np.abs((galflux)/(modflux) - 1)
        sumtidal=np.sum(tgal)
        #pixcountid=np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskm])
        pixcountid=np.size(immodel[maskm])


        if pixcountid > 0:
            tidal = (sumtidal / pixcountid)
        else:
            tidal=-1

        # bumpiness

        resbump  = (galfluxbum - modfluxbum)**2
        # sflux    = np.sum(np.abs(modfluxbum - sky))
        sflux    = np.sum(np.abs(modfluxbum))
        varres   = sigfluxbum * sigfluxbum
        numbump  = np.sum(resbump - varres)

        #pixcountbum=np.size(immodel[ylo - 1:yhi, xlo - 1:xhi][maskbum])
        pixcountbum=np.size(immodel[maskbum])


        # Bumpiness
        if pixcountbum > 0:

            meansflux = sflux / pixcountbum
            meanres   = numbump / pixcountbum

            if (meanres < 0):
                meanres=0

            bump      = (np.sqrt(meanres)) / meansflux

        else:
            bump=-1


    # computing RSS: 
    #rss=(imres[ylo - 1:yhi, xlo - 1:xhi][maskm]**2).sum()
    rss=(imres[maskm]**2).sum()
    
  
    return (tidal,objchinu,bump,snr,stdsnr,totsnr,rss,ndof)



def ExtractEllip(imagemat, idn, x, y, R, theta, ell, xmin, xmax, ymin, ymax):
    "This subroutine extract an ellipse within an box to compute photometric parameters "
    "It returns the area of ellipse with idn values. "


    xmin = np.int(xmin)
    xmax = np.int(xmax)
    ymin = np.int(ymin)
    ymax = np.int(ymax)

    q = (1 - ell)
    bim = q * R

    theta = theta * np.pi / 180  # Rads!!!

    ypos, xpos = np.mgrid[ymin - 1:ymax, xmin - 1:xmax]

    dx = xpos - x
    dy = ypos - y

    landa = np.arctan2(dy, dx)

    mask = landa < 0
    if mask.any():
        landa[mask] = landa[mask] + 2 * np.pi

    landa = landa - theta

    angle = np.arctan2(np.sin(landa) / bim, np.cos(landa) / R)

    xell = x + R * np.cos(angle) * np.cos(theta) - bim * \
        np.sin(angle) * np.sin(theta)
    yell = y + R * np.cos(angle) * np.sin(theta) + bim * \
        np.sin(angle) * np.cos(theta)

    dell = np.sqrt((xell - x)**2 + (yell - y)**2)
    dist = np.sqrt(dx**2 + dy**2)

    mask = dist < dell
    imagemat[ypos[mask], xpos[mask]] = idn

    return imagemat



def GetSize(x, y, R, theta, q, ncol, nrow):
    "this subroutine get the maximun"
    "and minimim pixels for Kron and sky ellipse"
    
    #theta is measured from x-axis    
    #q = (1 - ell)
    bim = q * R

    theta = theta * (np.pi / 180)  # rads!!

# getting size

    xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                       (bim**2) * (np.sin(theta))**2)

    ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                       (bim**2) * (np.cos(theta))**2)

    mask = xmin < 1
    if mask.any():
        if isinstance(xmin,np.ndarray):
            xmin[mask] = 1
        else:
            xmin = 1

    mask = xmax > ncol

    if mask.any():
        if isinstance(xmax,np.ndarray):
            xmax[mask] = ncol
        else:
            xmax = ncol

    mask = ymin < 1
    if mask.any():
        if isinstance(ymin,np.ndarray):
            ymin[mask] = 1
        else:
            ymin = 1

    mask = ymax > nrow
    if mask.any():
        if isinstance(ymax,np.ndarray):
            ymax[mask] = nrow
        else:
            ymax = nrow


    return (int(xmin), int(xmax), int(ymin), int(ymax))

def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex))
    hdu.writeto(newfits, overwrite=True)

    return True


def NED(params, galpar, galcomps):
    "connect to NED database to obtain Gal Extinction and other variables"
    
    objname=params.objname
    band=params.band

    params.flagweb=True

    objname=params.objname

    nedweb="https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname="

    filened=params.namened

    if params.flagmod or params.flagmag or params.flagscale or params.flagdim:

        GalExt=params.InMagCor
        DistMod=params.InDistMod
        DistMod2=params.InDistMod
        Scalekpc=params.InScale
        SbDim=params.InSbDim
    else:
        #checar si el archivo existe para no hacer conexion a internet
        if(not(os.path.isfile(filened))):

            # command for wget
            #wget -O NED_51.xml "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=m+51"
            wgetcmd = 'wget -O {} "{}{}"'.format(filened,nedweb,objname)

            print("Running: ",wgetcmd)

            errwg = sp.run([wgetcmd], shell=True, stdout=sp.PIPE,stderr=sp.PIPE, universal_newlines=True)

            if errwg.returncode != 0:
                print("can't connect to NED webserver. Is your internet connection working? ")
                print("Luminosity and absolute magnitude will not be computed") 
                params.flagweb=False
        else:
            print("using existing {} file ".format(filened))


        print("reading ",filened)
        votable=parse(filened,pedantic=False) 

        try: 
            table=votable.get_table_by_index(0) 

        except: 
            print("I can't read file or object name can be found in file")
            print("check object name or delete NED file for a new web search")
            print("luminosity and absolute magnitude will not be computed")
            params.flagweb=False 


        if params.flagweb==True:
            # si flag == True
            tablephot=votable.get_table_by_id("NED_DerivedValuesTable") 

            dataphot=tablephot.array

            lumdist=dataphot["luminosity_distance"].data[0] # units in Mpc

            #DistMod= 5 * np.log10(lumdist/10) 

            DistMod=dataphot["luminosity_distance_moduli"].data[0] # units in magnitudes

            tablext=votable.get_table_by_id("NED_BasicDataTable") 
            dataext=tablext.array

            extband="gal_extinc_" + band 
            try: 
                GalExt=dataext[extband].data[0] 
            except: 
                print("can't found {} in {} check filter name. GalExt=0 ".format(extband,filened))           
                GalExt=0


            print("Luminosity distance: (Mpc) ",lumdist)

            print("Module Distance (mag): ",DistMod)


            print("Galactic Extinction for band {} : {}".format(band,GalExt))


            Scalekpc=dataphot["cosmology_corrected_scale_kpc/arcsec"].data[0] # to convert to kpc

            print("Scale kpc/arcsec",Scalekpc)


            SbDim=dataphot["surface_brightness_dimming_mag"].data[0] 

            print("SB dimming in mag",SbDim)


            ## modulo de distancia calculado en forma independiente del redshift: 
            tabledist=votable.get_table_by_id("Redshift_IndependentDistances") 

            datadist=tabledist.array

            DistMod2=datadist["DistanceModulus"].data[0] 

            DistMod2=float(DistMod2)

            print("Distance Modulus (z independent) ",DistMod2)
        else:
            GalExt=0
            DistMod=0
            DistMod2=0
            Scalekpc=0
            SbDim=0



    return (GalExt,DistMod,DistMod2,Scalekpc,SbDim)


def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

## solve the Sersic equation
# to get the dependence of K over
# Sersic index

    count = 1

    #limits
    lima=0
    limb=100

#fx is the function to solve
    fxa = fx(n,lima)
    fxb = fx(n,limb)

    resk= (lima + limb)/2

    fxres=fx(n,resk)


    if(fxa * fxb < 0):

        while(np.abs(fxres) > 0.00000001):

            if(fxa * fxres > 0):
                lima=resk
            elif(fxa * fxres < 0):
                limb=resk
            elif(fxres==0):
                break
            resk= (lima + limb)/2
            fxres=fx(n,resk)

            count+=1

            if (count >= 10000):
                break

    else:
        print("no solution in the range: ({},{})\n".format(lima,limb))

    return (resk)


def fx(n,k):
    "function to solve to get the relation between Sersic index and K"


    func = np.exp(scipy.special.gammaln(2*n)) - 2 * np.exp(scipy.special.gammaln(2*n)) * scipy.special.gammainc(2*n,k)


    return(func)


def GetKAprox(n):
    "Aproximation to solve the dependence of K on the Sersic index"


    K = 2 * n - 1/3 + 4/(405*n) + 46 / (25515*n**2) + 131 / (1148175 * n**3) - 2194697 / (30690717750*n**4)


    return (K)



##############################################
##############################################
## Functions no longer used in this script: ##
##############################################
##############################################
##       deprecated Functions               ##

def ReadGauss(xpos,ypos,inputf):

    flagauss = False
    maskimage = ""

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    ngauss=1

    magas=[]
    fwhmgas=[]
    qgas=[]
    pagas=[]

    magas=np.array(magas)
    fwhmgas=np.array(fwhmgas)
    qgas=np.array(qgas)
    pagas=np.array(pagas)

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
        if tmp[0] == "A)":     # input image
            inputimage=tmp[1]
            exptime=GetExpTime(inputimage)

        if tmp[0] == "B)":     # out image
            outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            maskimage=tmp[1]

        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            scale=float(tmp[1])

        if tmp[0] == "0)" and tmp[1] == "gaussian":     # plate scale
        #  flagexp == True forces to take the value of exp instead of gauss

            flagauss=True
            xcgas=0
            ycgas=0
            q=0
            pa=0
            flagdist=False
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgas=float(tmp[1])
                    ycgas=float(tmp[2])

                    # correcting coordinates
                    xcgas=xcgas-xmin+1
                    ycgas=ycgas-ymin+1


                    dist = np.sqrt((xpos-xcgas)**2+(ypos-ycgas)**2)
                    if (dist < 5):
                        flagdist=True

                if (tmp[0] == "3)") and (flagdist==True):    # axis ratio
                    magauss=float(tmp[1])

                if (tmp[0] == "4)") and (flagdist==True):    # axis ratio
                    fwhmgauss=float(tmp[1])

                if (tmp[0] == "9)") and (flagdist==True):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (flagdist==True): # position angle
                    pa=float(tmp[1])

            if (flagdist == True):


            # saving data
                magas=np.append(magas,magauss)
                fwhmgas=np.append(fwhmgas,fwhmgauss)
                qgas=np.append(qgas,q)
                pagas=np.append(pagas,pa)

        index += 1

    GalfitFile.close()

    return magas,fwhmgas,qgas,pagas

def ReadSersic(xpos,ypos,inputf):

    flagser = False

    #   init values
    maskimage = ""

    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    ngauss=1

    magser=[]
    reser=[]
    nser=[]
    qser=[]
    paser=[]

    magser=np.array(magser)
    reser=np.array(reser)
    nser=np.array(nser)
    qser=np.array(qser)
    paser=np.array(paser)

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
        if tmp[0] == "A)":     # input image
            inputimage=tmp[1]
            exptime=GetExpTime(inputimage)

        if tmp[0] == "B)":     # out image
            outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            maskimage=tmp[1]

        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            scale=float(tmp[1])

        if tmp[0] == "0)" and tmp[1] == "sersic":     # plate scale

            flagser=True
            xcgas=0
            ycgas=0
            q=0
            pa=0
            flagdist=False
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcser=float(tmp[1])
                    ycser=float(tmp[2])

                    # correcting coordinates
                    xcser=xcser-xmin+1
                    ycser=ycser-ymin+1


                    dist = np.sqrt((xpos-xcser)**2+(ypos-ycser)**2)
                    if (dist < 5):
                        flagdist=True

                if (tmp[0] == "3)") and (flagdist==True):    # axis ratio
                    mag=float(tmp[1])

                if (tmp[0] == "4)") and (flagdist==True):    # axis ratio
                    re=float(tmp[1])

                if (tmp[0] == "5)") and (flagdist==True):    # axis ratio
                    n=float(tmp[1])

                if (tmp[0] == "9)") and (flagdist==True):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (flagdist==True): # position angle
                    pa=float(tmp[1])

            if (flagdist == True):


            # saving data
                magser=np.append(magser,mag)
                reser=np.append(reser,re)
                nser=np.append(nser,n)
                qser=np.append(qser,q)
                paser=np.append(paser,pa)


        index += 1

    GalfitFile.close()

    return magser,reser,nser,qser,paser


def ReadExp(xpos,ypos,inputf):

    flagexp = False
    maskimage = ""


    GalfitFile = open(inputf,"r")

    # All lines including the blank ones
    lines = (line.rstrip() for line in GalfitFile)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)

    index = 0
    ngauss=1

    magexp=[]
    rsexp=[]
    qexp=[]
    paexp=[]

    magexp=np.array(magexp)
    rsexp=np.array(rsexp)
    qexp=np.array(qexp)
    paexp=np.array(paexp)

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
        if tmp[0] == "A)":     # input image
            inputimage=tmp[1]
            exptime=GetExpTime(inputimage)

        if tmp[0] == "B)":     # out image
            outimage=tmp[1]

        if tmp[0] == "F)":     # mask image
            maskimage=tmp[1]

        if tmp[0] == "H)":     # region fit box
            xmin=int(tmp[1])
            xmax=int(tmp[2])
            ymin=int(tmp[3])
            ymax=int(tmp[4])

        if tmp[0] == "J)":     # mgzpt
            mgzpt=float(tmp[1])

        if tmp[0] == "K)":     # plate scale
            scale=float(tmp[1])

        if tmp[0] == "0)" and tmp[1] == "expdisk":     # plate scale
        #  flagexp == True forces to take the value of exp instead of gauss

            flagauss=True
            xcgas=0
            ycgas=0
            q=0
            pa=0
            flagdist=False
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcexp=float(tmp[1])
                    ycexp=float(tmp[2])

                    # correcting coordinates
                    xcexp=xcexp-xmin+1
                    ycexp=ycexp-ymin+1

                    dist = np.sqrt((xpos-xcexp)**2+(ypos-ycexp)**2)
                    if (dist < 5):
                        flagdist=True

                if (tmp[0] == "3)") and (flagdist==True):    # axis ratio
                    mag=float(tmp[1])

                if (tmp[0] == "4)") and (flagdist==True):    # axis ratio
                    rs=float(tmp[1])

                if (tmp[0] == "9)") and (flagdist==True):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (flagdist==True): # position angle
                    pa=float(tmp[1])

            if (flagdist == True):


            # saving data
                magexp=np.append(magexp,mag)
                rsexp=np.append(rsexp,rs)
                qexp=np.append(qexp,q)
                paexp=np.append(paexp,pa)


        index += 1

    GalfitFile.close()

    return magexp,rsexp,qexp,paexp



def PloTotal(magas,fwhmgas,qgas,pagas,magser,reser,nser,qser,paser,magexp,rsexp,qexp,paexp,patot,xlim):


    xgauss =  np.arange(xlim[0],xlim[1],0.1)


    radx,ytot=GalTotal(magas,fwhmgas,qgas,pagas,magser,reser,nser,qser,paser,magexp,rsexp,qexp,paexp,patot,xgauss)

    strtot="init Model"
    plt.plot(radx, ytot,'--',color='blue',markersize=0.7,label=strtot)



def PlotMulTotal(magas,fwhmgas,qgas,pagas,magser,reser,nser,qser,paser,magexp,rsexp,qexp,paexp,patot,ax,row,xlim,flaglogx):


    xgauss =  np.arange(xlim[0],xlim[1],0.1)


    radx,ytot=GalTotal(magas,fwhmgas,qgas,pagas,magser,reser,nser,qser,paser,magexp,rsexp,qexp,paexp,patot,xgauss)

    strtot="init Model"

    if flaglogx == False:
        ax[row,0].plot(radx, ytot,'--',color='blue',markersize=0.7,label=strtot)
    else:
        ax[row,0].semilogx(radx, ytot,'--',color='blue',markersize=0.7,label=strtot)


def PlotMulGauss(magas,fwhmgas,qgas,pagas,angle,ax,row,xlim,flaglogx):


    xgauss =  np.arange(xlim[0],xlim[1],0.1)

    ngauss=1

    alpha= angle - pagas

    for idx, item in enumerate(magas):


        radx,ygauss=GalGauss(magas[idx],fwhmgas[idx],qgas[idx],pagas[idx],alpha[idx],xgauss)
        strgas="gauss " + str(ngauss)
        if flaglogx == False:
            ax[row, 0].plot(radx, ygauss,'--',color='green',markersize=0.7,label=strgas)
        else:
            ax[row, 0].semilogx(radx, ygauss,'--',color='green',markersize=0.7,label=strgas)

        ngauss=ngauss+1




def PlotGauss(magas,fwhmgas,qgas,pagas,angle,xlim):


    xgauss =  np.arange(xlim[0],xlim[1],0.1)

    ngauss=1

    alpha= angle - pagas


    for idx, item in enumerate(magas):

        radx,ygauss=GalGauss(magas[idx],fwhmgas[idx],qgas[idx],pagas[idx],alpha[idx],xgauss)


        strgas="gauss " + str(ngauss)
        plt.plot(radx, ygauss,'--',color='green',markersize=0.7,label=strgas)
        ngauss=ngauss+1


def PlotSersic(magser,reser,nser,qser,paser,angle,xlim):

    xser =  np.arange(xlim[0],xlim[1],0.1)

    alpha= angle - paser

    num=1

    for idx, item in enumerate(magser):

        radx,yser=GalSersic(magser[idx],reser[idx],nser[idx],qser[idx],paser[idx],alpha[idx],xser)
        strgas="sersic " + str(num)
        plt.plot(radx, yser,'--',color='red',markersize=0.7,label=strgas)
        num=num+1



def PlotMulSersic(magser,reser,nser,qser,paser,angle,ax,row,xlim,flaglogx):

    xser =  np.arange(xlim[0],xlim[1],0.1)

    num=1

    alpha= angle - paser

    for idx, item in enumerate(magser):
        radx,yser=GalSersic(magser[idx],reser[idx],nser[idx],qser[idx],paser[idx],alpha[idx],xser)
        strgas="sersic " + str(num)
        if flaglogx == False:
            ax[row, 0].plot(radx, yser,'--',color='red',markersize=0.7,label=strgas)
        else:
            ax[row, 0].semilogx(radx, yser,'--',color='red',markersize=0.7,label=strgas)

        num=num+1


def PlotExp(magexp,rsexp,qexp,paexp,angle,xlim):


    xexp =  np.arange(xlim[0],xlim[1],0.1)


    alpha= angle - paexp


    num=1

    for idx, item in enumerate(magexp):

        radx,yexp=GalExp(magexp[idx],rsexp[idx],qexp[idx],paexp[idx],alpha[idx],xexp)

        strgas="exponential " + str(num)
        plt.plot(radx, yexp,'--',color='blue',markersize=0.7,label=strgas)
        num=num+1


def PlotMulExp(magexp,rsexp,qexp,paexp,angle,ax,row,xlim,flaglogx):


    xexp =  np.arange(xlim[0],xlim[1],0.1)

    alpha= angle - paexp


    num=1

    for idx, item in enumerate(magexp):

        radx,yexp=GalExp(magexp[idx],rsexp[idx],qexp[idx],paexp[idx],alpha[idx],xexp)
        strgas="exponential " + str(num)

        if flaglogx == False:
            ax[row,0].plot(radx, yexp,'--',color='blue',markersize=0.7,label=strgas)
        else:
            ax[row,0].semilogx(radx, yexp,'--',color='blue',markersize=0.7,label=strgas)

        num=num+1




def GalTotal(magauss,fwhmgauss,qgauss,pagauss,magser,reser,nser,qser,paser,magexp,rsexp,qexp,paexp,patot,rad):

    It=0

    ## rest angle
    anglegas= patot
    angleser= patot
    angleexp= patot
    # changing to rads!
    pagauss=pagauss*np.pi/180
    paser=paser*np.pi/180
    paexp=paexp*np.pi/180
    anglegas=anglegas*np.pi/180
    angleser=angleser*np.pi/180
    angleexp=angleexp*np.pi/180
    #

    ##  gauss
    sigma=fwhmgauss/2.354
    mags=(-magauss)/2.5
    ftot=10**(mags)
    I0=ftot/(2*np.pi*qgauss*sigma**2)
    ###

    ##  Sersic
    kser=GetKAprox(nser)
    mags=(-magser)/2.5
    ftot=10**(mags)
    Ie=ftot/(2*np.pi*(reser**2)*np.exp(kser)*nser*kser**(-2*nser)*(scipy.special.gamma(2*nser)*qser))
    ##



    ## exponential
    mags=(-magexp)/2.5
    ftot=10**(mags)
    Is=ftot/(2*np.pi*(rsexp**2)*qexp)

    ###############################################
    for idx, item in enumerate(magauss):
        alpha1=np.cos(anglegas)*np.cos(pagauss[idx])+np.sin(anglegas)*np.sin(pagauss[idx])
        alpha2=np.cos(anglegas)*np.sin(pagauss[idx])-np.sin(anglegas)*np.cos(pagauss[idx])

        agas2=(rad**2)*(alpha1**2+(1/(qgauss[idx]**2))*alpha2**2)
        Irgas=I0[idx]*np.exp(-(agas2)/(2*sigma[idx]**2))

        It=It+Irgas

    for idx, item in enumerate(magser):
        alpha1=np.cos(angleser)*np.cos(paser[idx])+np.sin(angleser)*np.sin(paser[idx])
        alpha2=np.cos(angleser)*np.sin(paser[idx])-np.sin(angleser)*np.cos(paser[idx])

        aser2=(rad**2)*(alpha1**2+(1/(qser[idx]**2))*alpha2**2)
        aser=np.sqrt(aser2)
        Irser=Ie[idx]*np.exp(-kser[idx]*((aser/reser[idx])**(1/nser[idx]) - 1 ))

        It=It+ Irser

    for idx, item in enumerate(magexp):

    #Exponential
        alpha1=np.cos(angleexp)*np.cos(paexp[idx])+np.sin(angleexp)*np.sin(paexp[idx])
        alpha2=np.cos(angleexp)*np.sin(paexp[idx])-np.sin(angleexp)*np.cos(paexp[idx])

        aexp2=(rad**2)*(alpha1**2+(1/(qexp[idx]**2))*alpha2**2)
        aexp=np.sqrt(aexp2)
        Irexp=Is*np.exp(-aexp/rsexp[idx])

        It = It + Irexp

    Magt = -2.5*np.log10(It)

    return rad,Magt

def GalSersic(magser,reser,nser,qser,paser,angle,radx):

    kser=GetKAprox(nser)


    mags=(-magser)/2.5
    ftot=10**(mags)

    angle=angle*np.pi/180
    paser=paser*np.pi/180

    Ie=ftot/(2*np.pi*(reser**2)*np.exp(kser)*nser*kser**(-2*nser)*(scipy.special.gamma(2*nser)*qser))

    xser=radx*np.cos(angle)*np.cos(paser) - radx*qser*np.sin(angle)*np.sin(paser)
    yser=radx*np.cos(angle)*np.sin(paser) + radx*qser*np.sin(angle)*np.cos(paser)

    radang= np.sqrt(xser**2 + yser**2)

    Ir=Ie*np.exp(-kser*((radx/reser)**(1/nser) - 1 ))

    yser = -2.5*np.log10(Ir)

    return radang,yser


# duplicated
#def GetKAprox(n):

 #   K=2*n-1/3+4/(405*n)+46/(25515*n**2)+131/(1148175*n**3)-2194697/(30690717750*n**4)


 #   return (K)




def GetK(n):
    "Solve the Sersic function to get the dependence of K over Sersic index"

    ## solve the Sersic equation
    # to get the dependence of K over
    # Sersic index

    count = 1

    #limits
    lima=0
    limb=100

    fxa = fx(n,lima)
    fxb = fx(n,limb)

    resk= (lima + limb)/2

    fxres=fx(n,resk)


    if(fxa * fxb < 0):

        while(np.abs(fxres) > 0.00000001):

            if(fxa * fxres > 0):
                lima=resk
            elif(fxa * fxres < 0):
                limb=resk
            elif(fxres==0):
                break
            resk= (lima + limb)/2
            fxres=fx(n,resk)

            count+=1

            if (count >= 10000):
                break

    else:
        print("no solution in the range: ({},{})\n".format(lima,limb))

    return (resk)


def fx(n,k):
    "function to solve to get the relation between Sersic index and K"


    func = np.exp(scipy.special.gammaln(2*n)) - 2 * np.exp(scipy.special.gammaln(2*n)) * scipy.special.gammainc(2*n,k)


    return(func)




def GalExp(magexp,rsexp,qexp,paexp,angle,radx):

    mags=(-magexp)/2.5
    ftot=10**(mags)

    Is=ftot/(2*np.pi*(rsexp**2)*qexp)


    angle=np.deg2rad(angle)
    paexp=np.deg2rad(paexp)


    xexp=radx*np.cos(angle)*np.cos(paexp) - radx*qexp*np.sin(angle)*np.sin(paexp)
    yexp=radx*np.cos(angle)*np.sin(paexp) + radx*qexp*np.sin(angle)*np.cos(paexp)

    radang= np.sqrt(xexp**2 + yexp**2)

    Ir=Is*np.exp(-radx/rsexp)

    yexp = -2.5*np.log10(Ir)

    return radang,yexp

def GalGauss(magauss,fwhmgauss,qgauss,pagauss,angle,radx):

    sigma=fwhmgauss/2.354


    mags=(-magauss)/2.5
    ftot=10**(mags)
    I0=ftot/(2*np.pi*qgauss*sigma**2)


    angle=np.deg2rad(angle)
    pagauss=np.deg2rad(pagauss)


    xgas=radx*np.cos(angle)*np.cos(pagauss) - radx*qgauss*np.sin(angle)*np.sin(pagauss)
    ygas=radx*np.cos(angle)*np.sin(pagauss) + radx*qgauss*np.sin(angle)*np.cos(pagauss)


    radang= np.sqrt(xgas**2 + ygas**2)

    Ir=I0*np.exp(-(radx**2)/(2*sigma**2))


    Magr = -2.5*np.log10(Ir)

    return radang,Magr

##############################################
##############################################
######## end of deprecated Functions  ########
##############################################
##############################################


if __name__ == '__main__':
    main()
