#! /usr/bin/env python3


import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path
import scipy
import matplotlib.pyplot as plt
import mimetypes

#import mgefit as mge
from mgefit.sectors_photometry import sectors_photometry


#from sys import argv
#option_handle_list = ['--script', '--lira_cbt', '--eur_hedge']
#options = {}
#for option_handle in option_handle_list:
#    options[option_handle[2:]] = argv[argv.index(option_handle) + 1] if option_handle in argv else None
#if options['eur_hedge'] == None:
    #Do x
#else:
    #Do y


def main():

    if (len(sys.argv[1:]) == 0):
        print ('Missing arguments')
        print ("Usage:\n %s [GALFITOutputFile] [--logx] [--q AxisRatio] [--pa PositionAngle] " % (sys.argv[0]))
        print ("Example:\n %s galfit.01 --logx" % (sys.argv[0]))
        print ("or Example:\n %s galfit.02 --q 0.35 --pa 60" % (sys.argv[0]))
        sys.exit()

#    galfile="galfit.01"
## reading arguments

    flaglogx=False
    flagq=False
    flagpa=False


    OptionHandleList = ['--logx', '--q', '--pa']
    options = {}
    for OptionHandle in OptionHandleList:
#        options[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1] if OptionHandle in sys.argv else None
        options[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)] if OptionHandle in sys.argv else None
    if options['logx'] != None:
        flaglogx=True
        print("X axis is logarithm")
    if options['q'] != None:
        flagq=True
#        qarg = np.float(options['q'])
    if options['pa'] != None:
        flagpa=True
#        parg=np.int(options['pa'])


    if flagpa == True:
        opt={}
        OptionHandle="--pa"
        opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
#        if opt['pa'] != None:
        parg=np.int(opt['pa'])

    if flagq == True:
        opt={}
        OptionHandle="--q"
        opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
#        if opt['q'] != None:
        qarg=np.float(opt['q'])









    galfile= sys.argv[1]



################################################
#################################################
##################################################

##############################
## default values  ignore this
############################

    # These parameters are given by find_galaxy for the mosaic image
    skylevel = 13.0
    sigmapsf = 0.4  # pixels
    eps = 0.28
    ang = 141.0  # major axis in the inner regions (gives a starting guess for the PA)
    xc = 974
    yc = 969
    ngauss = 11
    minlevel = 1.0
    scale = 0.100

#    numsectors=19
    #numsectors=36
    numsectors=15

    q=1-eps

    file = "ngc5831_f702w_mosaic.fits"

#########################################
#########################################
####### Read Galfit

    xc,yc,q,ang,skylevel,scale,file,mgzpt,exptime,mask=ReadGALFITout(galfile)

############
######################################
#    print(file)

    if flagq == True:
        q=qarg

    if flagpa == True:
#        q=qarg
        ang=parg

    str = "q = {} is used ".format(q)
    print(str)

    str = "pa = {} is used ".format(ang)
    print(str)


    (tmp)=file.split(".")

    namefile=tmp[0]

#    print(namefile)

# names for the different png
    namepng=namefile + ".png"
    namesec=namefile + "-gal.png"
    namemod=namefile + "-mod.png"
    namemul=namefile + "-mul.png"


    sbsky = mgzpt -2.5*np.log10(skylevel/exptime) + 2.5*np.log10(scale**2)


# hdu 1 => image   hdu 2 => model

    errmsg="file {} does not exist".format(file)

    assert os.path.isfile(file), errmsg

    hdu = fits.open(file)
    img = hdu[1].data
    model = hdu[2].data
    hdu.close()


    minlevel=-100  # minimun value for sky

    xradq,ysbq,ysberrq=elipsectors(img, mgzpt, exptime, scale, xc, yc, q, ang, skylevel=skylevel,
              n_sectors=numsectors, badpixels=mask, minlevel=minlevel, plot=1,nameplt=namesec)

    xradm,ysbm,ysberrm=elipsectors(model, mgzpt, exptime, scale, xc, yc, q, ang, skylevel=skylevel,
              n_sectors=numsectors, badpixels=mask, minlevel=minlevel, plot=1,nameplt=namemod)



    PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,sbsky,flaglogx)


    plt.savefig(namepng)


    Mulelipsectors(img, model, mgzpt, exptime, scale, xc, yc, q,
        ang, flaglogx,  skylevel=skylevel, n_sectors=numsectors, badpixels=mask, minlevel=minlevel, plot=1,nameplt="SectorGalaxy.png")

    plt.savefig(namemul)

    # deleting temp mask


    if mask != None:
        os.remove(mask)



##############
##############

def findrad(xarcq, ymgeq, numsectors):
# the xarcq array must be ordered
# use counts instead of mag
    xradq=[]
    ycntq=[]
    ycnterrq=[]
    xradq=np.array(xradq)
    ycntq=np.array(ycntq)
    ycnterrq=np.array(ycnterrq)

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
        ycntq=np.append(ycntq,np.mean(ymgeq[lima:limb]))
        ycnterrq=np.append(ycnterrq,np.std(ymgeq[lima:limb]))

        n=n-numsectors


    return xradq, ycntq, ycnterrq


def findsb(xarcq, ymgeq, numsectors):
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


def coordinates(q, pos_ang, xc, yc, s):

    ang = np.radians(90 - pos_ang)              # x-axis is major axis
    x, y = np.ogrid[:s[0], :s[1]] - np.array([xc, yc])
    x, y = x*np.cos(ang) - y*np.sin(ang), x*np.sin(ang) + y*np.cos(ang)
    x2, y2 = x**2, y**2
    rad = np.sqrt(x2 + y2)                      # Radius
    rell = np.sqrt(x2 + y2/q**2)                # Elliptical radius
    ecc = np.arctan2(np.abs(y/q), np.abs(x))    # Eccentric anomaly [0, pi/2]

    return rad, rell, ecc



def elipsectors(img, mgzpt, exptime, plate, xc, yc, q, ang, skylevel=0, badpixels=None,
              n_sectors=19, mask=None, minlevel=0, plot=False, nameplt=None):

    img = img - skylevel

    if badpixels is not None:

        errmsg="file {} does not exist".format(badpixels)
        assert os.path.isfile(badpixels), errmsg

        hdu = fits.open(badpixels)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()
    else:
        maskb=None

    eps=1-q

    if plot:
        plt.clf()
        print("")

############
# I have to switch x and y values because they are different axes for
# numpy

    xtemp=xc
    xc=yc
    yc=xtemp

    ang=90-ang
######################

    s = sectors_photometry(img, eps, ang, xc, yc,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


    if plot:
        plt.savefig(nameplt)
        plt.pause(1)  # Allow plot to appear on the screen


###################################################

    stidx = np.argsort(s.radius)

    mgerad=s.radius[stidx]
    mgecount=s.counts[stidx]
    mgeangle=s.angle[stidx]
    mgeanrad=mgeangle*np.pi/180

    ab=q

    aellab= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarc=aellab*plate

# formula according to cappellary mge manual
    mgesb= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1


    stidxq = np.argsort(aellarc)


    xarc = aellarc[stidxq]
    ymge = mgesb[stidxq]

#############
#############  Function


    xrad, ysb, ysberr    = findsb(xarc, ymge, n_sectors)

################
###############

    return xrad, ysb, ysberr


def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,sbsky,flag):

    """
    Produces final best-fitting plot

    """

        # Select an x and y plot range that is the same for all plots
        #

    plt.clf()
###  set limits

    minrad = np.min(xradq)
    maxrad = np.max(xradq)
    mincnt = np.min(ysbq)
    maxcnt = np.max(ysbq)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
    yran = mincnt * (maxcnt/mincnt)**np.array([-0.05, +1.05])


    plt.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='blue',markersize=0.7,label="galaxy")
#    plt.xlim(-1, 25)
#    plt.ylim(13, 23)


    plt.xlim(xran)
    plt.ylim(yran)

    plt.xlabel("arcsec")
    plt.ylabel("mag/''")

    plt.gca().invert_yaxis()

    plt.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='green',markersize=0.7,label="Model")


    if flag == True:
        plt.xscale("log")

    plt.legend(loc=1)



    plt.pause(1)

######################################################
#####################################################
######################################################

def Mulelipsectors(img, model, mgzpt, exptime, plate, xc, yc, q, ang, flag, skylevel=0, badpixels=None,
              n_sectors=19, mask=None, minlevel=0, plot=False, nameplt=None):


    img = img - skylevel
    model = model - skylevel

    fignum=1

    if badpixels is not None:

        errmsg="file {} does not exist".format(badpixels)
        assert os.path.isfile(badpixels), errmsg

        hdu = fits.open(badpixels)
        mask = hdu[0].data
        maskb=np.array(mask,dtype=bool)
        maskbt=maskb.T
        hdu.close()
    else:
        maskb=None


    eps=1-q

    if plot:
        plt.clf()
        print("")


############
# I have to switch x and y values because they are different axes for
# numpy
    xtemp=xc
    xc=yc
    yc=xtemp

    ang=90-ang
######################

    sg = sectors_photometry(img, eps, ang, xc, yc,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


    sm = sectors_photometry(model, eps, ang, xc, yc,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


###################################################

    stidx = np.argsort(sg.radius)

#   galaxy
    mgerad=sg.radius[stidx]*plate
    mgecount=sg.counts[stidx]
    mgeangle=sg.angle[stidx]
    mgeanrad=mgeangle*np.pi/180

# model
    stidx = np.argsort(sm.radius)

    mgemodrad=sm.radius[stidx]*plate
    mgemodcount=sm.counts[stidx]
    mgemodangle=sm.angle[stidx]
    mgemodanrad=mgemodangle*np.pi/180



# formula according to cappellary mge manual
    mgesb= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1
################# Model:
    mgemodsb= mgzpt - 2.5*np.log10(mgemodcount/exptime) + 2.5*np.log10(plate**2) + 0.1


    minrad = np.min(mgerad)
    maxrad = np.max(mgerad)
    minsb = np.min(mgesb)
    maxsb = np.max(mgesb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
#    yran = minsb * (maxsb/minsb)**np.array([-0.05, +1.05])
    yran = minsb * (maxsb/minsb)**np.array([+1.05,-0.05])


    sectors = np.unique(mgeangle)
    n = sectors.size
    dn = int(round(n/6.))
    nrows = (n-1)//dn + 1 # integer division

    plt.clf()

    fig, ax = plt.subplots(nrows, 2, sharex=True, sharey='col', num=fignum)
    fig.subplots_adjust(hspace=0.01)

    fig.text(0.04, 0.5, 'Surface brightness', va='center', rotation='vertical')
    fig.text(0.96, 0.5, 'error (%)', va='center', rotation='vertical')

    ax[-1, 0].set_xlabel("arcsec")
    ax[-1, 1].set_xlabel("arcsec")

    row = 0
    for j in range(0, n, dn):
        w = np.nonzero(mgeangle == sectors[j])[0]
        w = w[np.argsort(mgerad[w])]
        r = mgerad[w]
        r2 = mgemodrad[w]

        txt = "$%.f^\circ$" % sectors[j]



        ax[row, 0].set_xlim(xran)
        ax[row, 0].set_ylim(yran)

        if flag == False:
            ax[row, 0].plot(r, mgesb[w], 'C0o')
            ax[row, 0].plot(r2, mgemodsb[w], 'C1-', linewidth=2)
        else:

            ax[row, 0].semilogx(r, mgesb[w], 'C0o')
            ax[row, 0].semilogx(r2, mgemodsb[w], 'C1-', linewidth=2)



        ax[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=ax[row, 0].transAxes)

        sberr=1-mgemodsb[w]/mgesb[w]

        ax[row, 1].plot(r, sberr*100, 'C0o')

        ax[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        ax[row, 1].yaxis.tick_right()
        ax[row, 1].yaxis.set_label_position("right")
        ax[row, 1].set_ylim([-19.5, 20])

        row += 1




#    return xrad, ysb, ysberr


def ReadGALFITout(inputf):


#  It obtains the xc,yc,pa, q from the first components. It ignore the rest
    flagser = True
    flagexp = True
    flagauss = True

#    inputf = "fit.log"

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

#    function="sersic"

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

        # first sersic component
        if tmp[0] == "0)" and tmp[1] == "sersic" and flagser == True:     # plate scale

            flagser=False

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xc=float(tmp[1])
                    yc=float(tmp[2])

                if tmp[0] == "9)":    # axis ratio
                    q=float(tmp[1])

                if tmp[0] == "10)": # position angle
                    pa=float(tmp[1])

# second component exponential model
        if tmp[0] == "0)" and tmp[1] == "expdisk" and flagexp == True:     # plate scale

            flagexp=False

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcexp=float(tmp[1])
                    ycexp=float(tmp[2])

                    if flagser == False:
                        dist = np.sqrt((xc-xcexp)**2+(yc-ycexp)**2)
                    else:
                        dist=0
                        xc=xcexp
                        yc=ycexp


                if (tmp[0] == "9)") and (dist < 5):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 5): # position angle
                    pa=float(tmp[1])

# check if a third component exists

        if tmp[0] == "0)" and tmp[1] == "gaussian" and flagauss == True and flagexp == True:     # plate scale
        #  flagexp == True forces to take the value of exp instead of gauss

            flagauss=False
            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgauss=float(tmp[1])
                    ycgauss=float(tmp[2])

                    if flagser == False:
                        dist = np.sqrt((xc-xcgauss)**2+(yc-ycgauss)**2)
                    else:
                        dist=0
                        xc=xcgauss
                        yc=ycgauss


                if (tmp[0] == "9)") and (dist < 5):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 5): # position angle
                    pa=float(tmp[1])



        if tmp[0] == "0)" and tmp[1] == "sky" :     # plate scale

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":    # axis ratio
                    skylevel=float(tmp[1])

        index += 1

    GalfitFile.close()


    errmsg="file {} does not exist".format(inputimage)
    assert os.path.isfile(inputimage), errmsg

    exptime=GetExpTime(inputimage)


    errmsg="xc and yc are unknown "
    assert ("xc" in locals()) and ("yc" in locals())  , errmsg


    # correcting coordinates
    xc=xc-xmin+1
    yc=yc-ymin+1

##   mask image

### checking sane values


    if os.path.isfile(maskimage):
        mime=mimetypes.guess_type(maskimage)

        flagbm = not(mime[0] == "text/plain")

        errmsg="Sorry the mask file: {}  must be binary, not ASCII ".format(maskimage)
        assert flagbm is True, errmsg


        mask="tempmask.fits"
        GetFits(maskimage, mask, xmin, xmax, ymin, ymax)

    else:
        errmsg="Mask file does not exist"
        print(errmsg)
        mask=None

    return xc,yc,q,pa,skylevel,scale,outimage,mgzpt,exptime,mask



def GetExpTime(Image):
    # k Check
    "Get exposition time from the image"

    hdu = fits.open(Image)
    exptime = hdu[0].header["EXPTIME"]
    hdu.close()
    return exptime


def GetFits(Image, Imageout, xlo, xhi, ylo, yhi):
    "Get a piece from the image"
# k Check


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




if __name__ == '__main__':
    main()
