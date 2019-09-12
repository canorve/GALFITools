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
import mimetypes

from mgefit.sectors_photometry import sectors_photometry


def main():

    if (len(sys.argv[1:]) == 0):
        print ('Missing arguments')
        print ("Usage:\n %s [GALFITOutputFile] [--logx] [--q AxisRatio] [--pa PositionAngle] [--sub] " % (sys.argv[0]))
        print ("GALFITOutputFile: GALFIT output file ")
        print ("logx: activates X-axis as logarithm ")
        print ("q: introduce axis ratio ")
        print ("pa: introduce position angle (same as GALFIT) ")
        print ("sub: plots subcomponents ")
#        print ("init: plots initial parameter model from subcomponents ")
        print ("Example:\n %s galfit.01 --logx" % (sys.argv[0]))
        print ("or Example:\n %s galfit.02 --q 0.35 --pa 60 --sub" % (sys.argv[0]))
        sys.exit()


    flaglogx=False
    flagq=False
    flagpa=False
    flagsub=False

# init values
    qarg=1
    parg=0

## init sub values
    Comps=np.array([False])
    N=0

    OptionHandleList = ['--logx', '--q', '--pa','--sub']
    options = {}
    for OptionHandle in OptionHandleList:
        options[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)] if OptionHandle in sys.argv else None
    if options['logx'] != None:
        flaglogx=True
        print("X axis is logarithm")
    if options['q'] != None:
        flagq=True
    if options['pa'] != None:
        flagpa=True
    if options['sub'] != None:
        flagsub=True
        print("Plotting subcomponents ")
#        print("Warning: Subcomponents and Model will not match within the convolution box ")

################################
    if flagpa == True:
        opt={}
        OptionHandle="--pa"
        opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        parg=np.int(opt['pa'])

    if flagq == True:
        opt={}
        OptionHandle="--q"
        opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        qarg=np.float(opt['q'])


    galfile= sys.argv[1]

    print("angle in multi-plot is measured from the galaxy's major axis ")



################################################
################################################
################################################

######################################
####### Read Galfit File #############

    xc,yc,q,ang,skylevel,scale,file,mgzpt,exptime,mask=ReadGALFITout(galfile)

######################################
######################################
#    print(file)

    if flagq == True:
        q=qarg

    if flagpa == True:
        ang=parg

## correction
#    ang=90+ang

###

    str = "q = {} is used ".format(q)
    print(str)

    str = "pa = {} is used ".format(ang)
    print(str)


    (tmp)=file.split(".")

    namefile=tmp[0]


# names for the different png

    namepng=namefile + ".png"
    namesec=namefile + "-gal.png"
    namemod=namefile + "-mod.png"
    namemul=namefile + "-mul.png"
    namesub=namefile + "-sub.fits"

# hdu 1 => image   hdu 2 => model

    errmsg="file {} does not exist".format(file)

    assert os.path.isfile(file), errmsg

    hdu = fits.open(file)
    img = hdu[1].data
    model = hdu[2].data
    hdu.close()


#   numsectors=19
#   numsectors=36
    numsectors=15
    minlevel=-100  # minimun value for sky


    limx,limy=EllipSectors(img, model, mgzpt, exptime, scale, xc, yc, q, ang, galfile, namefile, flaglogx, flagsub, skylevel=skylevel,
              n_sectors=numsectors, badpixels=mask, minlevel=minlevel, plot=1)


##############################################
##############################################
#### Creating Subcomponents Galfit output file
    Comps=[]

    if flagsub:

#        if (not(os.path.isfile(namesub))):


        print("running galfit to create subcomponents...")

        rungal = "galfit -o3 {}".format(galfile)
        errgal = sp.run([rungal], shell=True, stdout=sp.PIPE,
        stderr=sp.PIPE, universal_newlines=True)

        runchg = "mv subcomps.fits {}".format(namesub)
        errchg = sp.run([runchg], shell=True, stdout=sp.PIPE,
        stderr=sp.PIPE, universal_newlines=True)


    # read number of components
        Comps,N=ReadNComp(galfile,xc,yc)


        SubComp(namesub,N,Comps,mgzpt,exptime,scale,xc,yc,q,ang,skylevel=skylevel,
              n_sectors=numsectors, badpixels=mask, minlevel=minlevel)


##############################################
##############################################
##############################################


    plt.pause(1)
    plt.savefig(namepng)


########################################################
################ Multiplots: ###########################
########################################################

    MulEllipSectors(img, model, mgzpt, exptime, scale, xc, yc, q,
        ang, galfile, limx, flaglogx, flagsub, namesub, N, Comps, skylevel=skylevel, n_sectors=numsectors, badpixels=mask, minlevel=minlevel, plot=1,nameplt="SectorGalaxy.png")



    plt.pause(1)
    plt.savefig(namemul)

    if mask != None:
        os.remove(mask)

    print("Fin.")

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


def EllipSectors(img, model, mgzpt, exptime, plate, xc, yc, q, ang, galfile, namefile, flaglogx, flagsub, skylevel=0, badpixels=None,
              n_sectors=19, mask=None, minlevel=0, plot=False):

    namesec=namefile + "-gal.png"
    namemod=namefile + "-mod.png"

# removing background:
    img = img - skylevel
    model = model - skylevel

    xradm = []
    ysbm = []
    ysberrm = []

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
# numpy:
    yctemp=xc
    xctemp=yc
# and angle is different as well:
    angsec=90-ang
#    angsec=ang

####################

###############################
#  galaxy:

    g = sectors_photometry(img, eps, angsec, xctemp, yctemp,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


    if plot:
        plt.savefig(namesec)
        plt.pause(1)  # Allow plot to appear on the screen


###################################################

    stidxg = np.argsort(g.radius)

    mgerad=g.radius[stidxg]
    mgecount=g.counts[stidxg]
    mgeangle=g.angle[stidxg]
#    mgeanrad=(mgeangle)*np.pi/180
    mgeanrad=np.deg2rad(mgeangle)

#    mgeanrad=(90-mgeangle)*np.pi/180 #converting back angle


    ab=q

    aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

    aellarcg=aellabg*plate

####
# formula according to cappellary mge manual:
    mgesbg= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1
####

    stidxq = np.argsort(aellarcg)


    xarcg = aellarcg[stidxq]
    ymgeg = mgesbg[stidxq]

#############
#############  Function to order SB

    xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)

################
###############
###############################

    #  model:
    m = sectors_photometry(model, eps, angsec, xctemp, yctemp,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


    if plot:
        plt.savefig(namemod)
        plt.pause(1)  # Allow plot to appear on the screen


    ###################################################

    stidxm = np.argsort(m.radius)

    mgerad=m.radius[stidxm]
    mgecount=m.counts[stidxm]
    mgeangle=m.angle[stidxm]
#        mgeanrad=mgeangle*np.pi/180
    mgeanrad=np.deg2rad(mgeangle)

#        mgeanrad=(90-mgeangle)*np.pi/180


    ab=q

    aellabm= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)


    aellarcm=aellabm*plate

    # formula according to cappellary mge manual

    mgesbm= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1

    stidxq = np.argsort(aellarcm)


    xarcm = aellarcm[stidxq]
    ymgem = mgesbm[stidxq]

    #############
    #############  Function

    xradm, ysbm, ysberrm    = FindSB(xarcm, ymgem, n_sectors)

    ################
    ###############
    ################

    limx,limy=PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,flaglogx)

    ####### Read Gaussians from GALFIT
#    (magas,fwhmgas,qgas,pagas)=ReadGauss(xc,yc,galfile)
#    (magser,reser,nser,qser,paser)=ReadSersic(xc,yc,galfile)
#    (magexp,rsexp,qexp,paexp)=ReadExp(xc,yc,galfile)


#    if flagsub == True:

########  Run galfit to create sub components

#        PlotGauss(magas,fwhmgas*plate,qgas,pagas,ang,limx)
#        PlotSersic(magser,reser*plate,nser,qser,paser,ang,limx)
#        PlotExp(magexp,rsexp*plate,qexp,paexp,ang,limx)

#        PlotGauss(magas,fwhmgas*plate,qgas,pagas,0,limx)
#        PlotSersic(magser,reser*plate,nser,qser,paser,0,limx)
#        PlotExp(magexp,rsexp*plate,qexp,0,0,limx)

#        PloTotal(magas,fwhmgas*plate,qgas,pagas,magser,reser*plate,nser,qser,paser,magexp,rsexp*plate,qexp,paexp,ang,limx)


    plt.legend(loc=1)


    return limx,limy




def PlotSB(xradq,ysbq,ysberrq,xradm,ysbm,ysberrm,flag):
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

    plt.errorbar(xradq, ysbq,yerr=ysberrq,fmt='o-',capsize=2,color='red',markersize=0.7,label="galaxy")


    plt.xlim(xran)
    plt.ylim(yran)

    plt.xlabel("arcsec")
    plt.ylabel("mag/''")

    plt.gca().invert_yaxis()

    plt.errorbar(xradm, ysbm,yerr=ysberrm,fmt='o-',capsize=2,color='blue',markersize=0.7,label="Model")

    if flag == True:
        plt.xscale("log")

    return xran,yran


##########################################################
##########################################################
##########################################################

def SubComp(namesub,N,Comps,mgzpt,exptime,plate,xc,yc,q,ang,skylevel=0,
    n_sectors=19, badpixels=None, minlevel=0):

    errmsg="file {} does not exist".format(namesub)
    assert os.path.isfile(namesub), errmsg

    hdu = fits.open(namesub)

    subimgs=[]

    cnt=0  # image =0 do not count
    while(cnt<len(Comps)):

        if Comps[cnt] == True:
            img = hdu[cnt+2].data
            subimgs.append(img)

        cnt=cnt+1

    hdu.close()

    xradm = []
    ysbm = []
    ysberrm = []

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


############
# I have to switch x and y values because they are different axes for
# numpy:
    yctemp=xc
    xctemp=yc
# and angle is different as well:
    angsec=90-ang
#    angsec=ang

####################

###############################
#  galaxy:
    ab=q

    n=0
    while(n<N):

        subim=subimgs[n]


        ### checar niveles de minlevel
        scmp = sectors_photometry(subim, eps, angsec, xctemp, yctemp,minlevel=minlevel,
                plot=0, badpixels=maskb, n_sectors=n_sectors)

###################################################

        stidxg = np.argsort(scmp.radius)

        mgerad=scmp.radius[stidxg]
        mgecount=scmp.counts[stidxg]
        mgeangle=scmp.angle[stidxg]
        mgeanrad=np.deg2rad(mgeangle)


        aellabg= mgerad * np.sqrt((np.sin(mgeanrad)**2)/ab**2 + np.cos(mgeanrad)**2)

        aellarcg=aellabg*plate

# formula according to cappellary mge manual

        mgesbg= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1


        stidxq = np.argsort(aellarcg)

        xarcg = aellarcg[stidxq]
        ymgeg = mgesbg[stidxq]

        xradq, ysbq, ysberrq    = FindSB(xarcg, ymgeg, n_sectors)



        PlotSub(xradq,ysbq,n)

        n=n+1


    return  xradq,ysbq,n




def PlotSub(xradq,ysbq,nsub):
    """
    Produces subcomponent plot

    """

    substr="component "+np.str(nsub+1)

    plt.plot(xradq, ysbq,'--',color='cyan',linewidth=4,markersize=0.7,label=substr)

    plt.legend(loc=1)


###########################
##########################
############################


def MulEllipSectors(img, model, mgzpt, exptime, plate, xc, yc, q, ang, galfile, limx, flag, flagsub, namesub, N, Comps, skylevel=0, badpixels=None,
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


        ####### Read Gaussians from GALFIT
#    (magas,fwhmgas,qgas,pagas)=ReadGauss(xc,yc,galfile)
#    (magser,reser,nser,qser,paser)=ReadSersic(xc,yc,galfile)
#    (magexp,rsexp,qexp,paexp)=ReadExp(xc,yc,galfile)


############
# I have to switch x and y values because they are different axes for
# numpy
    xtemp=xc
    xc=yc
    yc=xtemp

    angsec=90-ang
######################

    sg = sectors_photometry(img, eps, angsec, xc, yc,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)


    sm = sectors_photometry(model, eps, angsec, xc, yc,minlevel=minlevel,
            plot=1, badpixels=maskb, n_sectors=n_sectors)

###################################################

    stidx = np.argsort(sg.radius)

#   galaxy
    mgerad=sg.radius[stidx]*plate
    mgecount=sg.counts[stidx]
    mgeangle=sg.angle[stidx]
    mgeanrad=np.deg2rad(mgeangle)


# model

    stidx = np.argsort(sm.radius)

    mgemodrad=sm.radius[stidx]*plate
    mgemodcount=sm.counts[stidx]
    mgemodangle=sm.angle[stidx]
    mgemodanrad=np.deg2rad(mgemodangle)



# formula according to cappellary mge manual
    mgesb= mgzpt - 2.5*np.log10(mgecount/exptime) + 2.5*np.log10(plate**2) + 0.1
################# Model:

    mgemodsb= mgzpt - 2.5*np.log10(mgemodcount/exptime) + 2.5*np.log10(plate**2) + 0.1


    if flagsub:

        wtemp=[]
        mgesbsub=[]
        mgeradsub=[]
        mgeanglesub=[]
        rsub=[]
        hdusub = fits.open(namesub)

        subimgs=[]

        cnt=0
        while(cnt<len(Comps)):

            if Comps[cnt] == True:
                imgsub = hdusub[cnt+2].data
                subimgs.append(imgsub)
            cnt=cnt+1

        hdu.close()


    ###############################
    #  galaxy:
        ab=q
        ni=0
        while(ni<N):

            subim=subimgs[ni]

            subcmp = sectors_photometry(subim, eps, angsec, xc, yc, minlevel=minlevel,
                    plot=0, badpixels=maskb, n_sectors=n_sectors)

            subidx = np.argsort(subcmp.radius)

            temprad=subcmp.radius[subidx]*plate
            mgecountsub=subcmp.counts[subidx]

            tempangle=subcmp.angle[subidx]
            mgeanradsub=np.deg2rad(tempangle)


    # formula according to cappellary mge manual
            tempmge= mgzpt - 2.5*np.log10(mgecountsub/exptime) + 2.5*np.log10(plate**2) + 0.1

            mgesbsub.append(tempmge)
            mgeradsub.append(temprad)
            mgeanglesub.append(tempangle)

            ni+=1


    minrad = np.min(mgerad)
    maxrad = np.max(mgerad)
    minsb = np.min(mgesb)
    maxsb = np.max(mgesb)
    xran = minrad * (maxrad/minrad)**np.array([-0.02, +1.02])
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

#    row = 0
    row = 7
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

#        if flagsub:
#            ii=0
#            while(ii<N):
#                tsub = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
#                tsub = tsub[np.argsort(mgeradsub[ii][tsub])]

#                if (len(mgeradsub[ii]) < len(mgerad)):
#                    rtsub = mgeradsub[ii][tsub]
#                    rsub.append(rtsub)
#                else:
#                    tsub=w  # agregado
#                    rtsub = mgeradsub[ii][tsub]
#                    rsub.append(rtsub)
#                ii+=1
#                wsub.append(tsub)   # wsub corresponde a ii

        txtang= sectors[j]
        txt = "$%.f^\circ$" % txtang

        ax[row, 0].set_xlim(xran)
        ax[row, 0].set_ylim(yran)

        if flag == False:
            ax[row, 0].plot(r, mgesb[w], 'C3o')

            ax[row, 0].plot(r2, mgemodsb[wmod], 'C0-', linewidth=2)

        else:
            ax[row, 0].semilogx(r, mgesb[w], 'C3o')

            ax[row, 0].semilogx(r2, mgemodsb[wmod], 'C0-', linewidth=2)

#########################


        if flagsub == True:
#            PlotMulGauss(magas,fwhmgas*plate,qgas,pagas,90-sectors[j],ax,row,limx,flag)
#            PlotMulSersic(magser,reser*plate,nser,qser,paser,90-sectors[j],ax,row,limx,flag)
#            PlotMulExp(magexp,rsexp*plate,qexp,paexp,90-sectors[j],ax,row,limx,flag)

            ii=0
            while(ii<N):

                wtemp = np.nonzero(mgeanglesub[ii] == sectors[j])[0]
                wtemp = wtemp[np.argsort(mgeradsub[ii][wtemp])]

                rtemp = mgeradsub[ii][wtemp]


                if flag == False:

                    ax[row, 0].plot(rtemp, mgesbsub[ii][wtemp], 'C9--', linewidth=2)


                else:

                    ax[row, 0].semilogx(rtemp, mgesbsub[ii][wtemp], 'C9--', linewidth=2)


                ii+=1

#            SubMulComp(namesub,N,Comps,mgzpt,exptime,scale,xc,yc,q,ang,skylevel=skylevel,
#                n_sectors=numsectors, badpixels=mask, minlevel=minlevel)

        ax[row, 0].text(0.98, 0.95, txt, ha='right', va='top', transform=ax[row, 0].transAxes)

        if (len(mgemodrad) < len(mgerad)):
            sberr=1-mgemodsb[wmod]/mgesb[wmod]
            ax[row, 1].plot(r2, sberr*100, 'C0o')
        else:
            sberr=1-mgemodsb[w]/mgesb[w]
            ax[row, 1].plot(r, sberr*100, 'C0o')


        ax[row, 1].axhline(linestyle='--', color='C1', linewidth=2)
        ax[row, 1].yaxis.tick_right()
        ax[row, 1].yaxis.set_label_position("right")
        ax[row, 1].set_ylim([-19.5, 20])

#        row += 1
        row -= 1

#    return xrad, ysb, ysberr






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



def ReadGALFITout(inputf):

    flagfirst = True

    maskimage = ""

    skylevel=0

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


        # first component
        if tmp[0] == "0)" and flagfirst == True:     # plate scale

            flagfirst=False

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
                        dist = np.sqrt((xc-xcser)**2+(yc-ycser)**2)



                if tmp[0] == "9)" and (dist < 3):    # axis ratio
                    q=float(tmp[1])

                if tmp[0] == "10)" and (dist < 3): # position angle
                    pa=float(tmp[1])

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
                        dist = np.sqrt((xc-xcexp)**2+(yc-ycexp)**2)


                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
                    pa=float(tmp[1])


        if tmp[0] == "0)" and tmp[1] == "gaussian":

            while (tmp[0] != "Z)"):

                index += 1
                line = lines[index]
                (tmp) = line.split()

                if tmp[0] == "1)":   # center
                    xcgauss=float(tmp[1])
                    ycgauss=float(tmp[2])

                    if flagfirst == False:
                        dist = np.sqrt((xc-xcgauss)**2+(yc-ycgauss)**2)

                if (tmp[0] == "9)") and (dist < 3):    # axis ratio
                    q=float(tmp[1])

                if (tmp[0] == "10)") and (dist < 3): # position angle
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


    print("center is at xc, yc = ",xc,yc)


    # correcting coordinates
    xc=xc-xmin+1
    yc=yc-ymin+1


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

##################
##################


def ReadNComp(inputf,X,Y):
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

    Ncomp=0
    Comps=[]
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

#            print("xc, yc =",X,Y)


        if tmp[0] == "0)":

            if (tmp[1] != "sky"):
                while (tmp[0] != "Z)"):

                    index += 1
                    line = lines[index]
                    (tmp) = line.split()

                    if tmp[0] == "1)":   # center
                        xc=float(tmp[1])
                        yc=float(tmp[2])

                        dist = np.sqrt((xc-X)**2+(yc-Y)**2)
                        if (dist < 3):
                            Comps=np.append(Comps,True)
                            Ncomp=Ncomp + 1
                        else:
                            Comps=np.append(Comps,False)

                    if tmp[0] == "9)":    # axis ratio
                        q=float(tmp[1])

                    if tmp[0] == "10)": # position angle
                        pa=float(tmp[1])
            else:
                Comps=np.append(Comps,False)
#                Ncomp=Ncomp + 1


        index += 1

    GalfitFile.close()


    print("Number of components = ",Ncomp)



    return Comps,Ncomp

######
### This function is no longer used in this program
######
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


######################
######################




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



# function no longer used in this program
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


## function no longer used in this program
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
        #  flagexp == True forces to take the value of exp instead of gauss

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
#                    magauss = magauss + mgzpt

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




def GetKAprox(n):

    K=2*n-1/3+4/(405*n)+46/(25515*n**2)+131/(1148175*n**3)-2194697/(30690717750*n**4)


    return (K)




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



## function no longer used in this program
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

######




####################

## function no longer used in this  program

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
#                    magauss = magauss + mgzpt

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


def GetExpTime(Image):
    # k Check
    "Get exposition time from the image"

    hdu = fits.open(Image)
#    exptime = hdu[0].header["EXPTIME"]
    exptime = hdu[0].header.get("EXPTIME",1) # return 1 if not found

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
