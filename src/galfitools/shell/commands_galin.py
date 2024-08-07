

import argparse
import subprocess as sp
import sys
from collections import Counter

import numpy as np
from galfitools.galin.checkGalFile import checkFile
from galfitools.galin.galfit import GalComps, Galfit, GalHead, galfitLastFit
from galfitools.galin.getBoxSizeDs9 import getBoxSizeDs9
from galfitools.galin.getSersic import getSersic
from galfitools.galin.getStar import getStar
from galfitools.galin.imarith import imarith
from galfitools.galin.initgal import InitGal
from galfitools.galin.MakeMask import makeMask
from galfitools.galin.MakePSF import makePSF
from galfitools.galin.MaskDs9 import maskDs9
from galfitools.galin.MaskSky import skyRem
from galfitools.galin.xy2fits import xy2fits
from galfitools.mge.mge2galfit import (PrintHeader, PrintSersic, PrintSky,
                                       mge2gal)
from galfitools.shell.prt import printWelcome


def mainGetStar():

    printWelcome()

    parser = argparse.ArgumentParser(description="gets a image slice centered on the object peak")

    parser.add_argument("image", help="the image file to obtain the slice")
    parser.add_argument("Ds9regFile", help="the DS9 ellipse region file containing the star")
    parser.add_argument("size", type=int, help="the size of the new image in pixels")
    parser.add_argument("-c","--center", action="store_true", 
                        help="uses the center given in DS9 region file," + 
                        "otherwise it will find the x,y peaks within DS9 ellipse")
    parser.add_argument("-s","--sky", type=float, 
                        help="the sky background to be removed")
    parser.add_argument("-o","--out", type=str, 
                        help="the image output.",default="star.fits")

    parser.add_argument("-sig","--sigma", type=str, 
                        help="introduce the sigma image")

    parser.add_argument("-so","--sigout", type=str, 
                        help="the sigma image output.",default="sigma.fits")



    args = parser.parse_args()

    image = args.image
    regfile = args.Ds9regFile 
    imsize = args.size
    center = args.center

    sky = args.sky
    imout = args.out

    sigma = args.sigma
    sigout = args.sigout

    getStar(image, regfile, imsize, center, sky, imout, sigma, sigout)

    print("Done. Object fits file: {} created ".format(imout))

    if sigma:
        print("Done. sigma fits file: {} created ".format(sigout))


def mainInitGal(): 


    printWelcome()

   #parser argument section

    parser = argparse.ArgumentParser(description="Creates GALFIT's input files with different initial parameters")

    parser.add_argument("inFile", help="the galfit file galfit.XX ")
    parser.add_argument("-n","--number", type=int, help="the number of files generated. Default = 1",default=1)
        
    parser.add_argument("-p3","--param3", nargs=2, type=float, help="range of values to give to the 3) model's parameter in format [min max] ",)
    parser.add_argument("-p4", "--param4", nargs=2, type=float, help="range of values to give to the 4) model's parameter in format [min max] ")
    parser.add_argument("-p5", "--param5", nargs=2, type=float, help="range of values to give to the 5) model's parameter in format [min max] ")
    parser.add_argument("-p6", "--param6", nargs=2, type=float, help="range of values to give to the 6) model's parameter in format  [min max]")
    parser.add_argument("-p7", "--param7", nargs=2, type=float, help="range of values to give to the 7) model's parameter in format [min max]")
    parser.add_argument("-p8", "--param8", nargs=2, type=float, help="range of values to give to the 8) model's parameter in format [min max]")
    parser.add_argument("-p9", "--param9", nargs=2, type=float, help="range of values to give to the 9) model's parameter in format [min max]")
    parser.add_argument("-p10", "--param10", nargs=2, type=float, help="range of values to give to the 10) model's parameter in format [min max]")



    parser.add_argument("-nc","--numcomp", type=int, help="the component number which parameters will be changed")


    args = parser.parse_args()


    GalfitFile = args.inFile
    number = args.number
    param3 = args.param3
    param4 = args.param4
    param5 = args.param5
    param6 = args.param6
    param7 = args.param7

    param8 = args.param8
    param9 = args.param9
    param10 = args.param10


    numcomp = args.numcomp


    InitGal(GalfitFile, number, param3, param4, param5, param6, param7, param8, param9, param10, numcomp)


    print("rungalfit.sh has been created")


def mainMakeMask():


    printWelcome()

    parser = argparse.ArgumentParser(description="creates mask file from a catalog of Sextractor")

    parser.add_argument("Sexfile", help="Sextractor catalog file ")
    parser.add_argument("ImageFile", help="Image file")

    parser.add_argument("-o","--maskout", type=str, help="the output mask file name  ",default='masksex.fits')
    parser.add_argument("-sf","--satds9", type=str, help="ds9 saturation file",default='ds9sat.reg')
        
    parser.add_argument("-s","--scale", type=float, help="scale factor to increase the ellipses. Default=1",default = 1)
 


    args = parser.parse_args()


    sexfile = args.Sexfile 
    image = args.ImageFile  
    maskfile = args.maskout  
    scale = args.scale 
    satfileout = args.satds9





    makeMask(sexfile, image, maskfile, scale, satfileout)
    

    print("Done. Mask image created ")


def mainMaskDs9(): 

    printWelcome()

    parser = argparse.ArgumentParser(description="creates (or modify) a mask image for GALFIT using Ds9 regions such as Boxes, Ellipses and Polygons ")

    parser.add_argument("MaskFile", help="the Mask image file to modify or create")
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument("-f","--fill", type=int, help="the value in counts to fill into the Ds9 regions. Default = 0 (remove)",default=0)

    parser.add_argument("-i","--image", type=str, help="image to obtain the size  ")


    parser.add_argument("-b","--border", action="store_true", help="Mask the borders when their value is zero")
    parser.add_argument("-bv","--borValue",default=0,type=float, help="value of the border if it is different from zero")


    parser.add_argument("-sm","--skymean", type=float, help="sky mean to be used as a patch of sky")
    parser.add_argument("-sd","--skystd", type=float, help="sky standard deviation to be used with skymean as a patch of sky",default=1)


    args = parser.parse_args()

    MaskFile = args.MaskFile 
    RegFile = args.RegFile 
    fill = args.fill
    image = args.image

    bor_flag = args.border
    borValue = args.borValue

    skymean = args.skymean
    skystd = args.skystd

    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue, skymean=skymean, skystd=skystd) 

    print('Done. Mask {} created (or modified)'.format(MaskFile))



def mainMaskSky(): 


    printWelcome()

    parser = argparse.ArgumentParser(description="creates a mask image for GALFIT using original image and sky mean and sigma")


    parser.add_argument("ImageFile", help="original data image ")
    parser.add_argument("MaskFile", help="Name of the new Mask file")
    parser.add_argument("-sm","--skymean",default=0,type=float, help="mean of the sky background")
    parser.add_argument("-ss","--skysigma",default=0,type=float, help="sigma of the sky background")
    parser.add_argument("-ns","--numbersig",default=1,type=float, help="number of times that the sigma of the sky will be multiplied to remove the sky background")

    parser.add_argument("-r","--region", type=str, help="Ds9 region to remove from mask")

    parser.add_argument("-b","--border", action="store_true", help="Mask the borders when their value is zero")

    parser.add_argument("-bv","--borValue",default=0,type=float, help="value of the border if it is different from zero")
    args = parser.parse_args()

    image = args.ImageFile 
    mask = args.MaskFile
    sky_mean = args.skymean
    sky_sig = args.skysigma
    nsig = args.numbersig

    region = args.region

    bor_flag = args.border
    borValue = args.borValue


    skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)


    if region:
        maskDs9(mask, region, 0, None, False, 0)  #remove ds9 region if provided


    print('Done. Mask file {} created'.format(mask))


def mainxy2fits(): 

    printWelcome()

    parser = argparse.ArgumentParser(description="code to convert ASCII x,y positions to FTIS mask")

    parser.add_argument("ImageFile", help="The Image file ")
    parser.add_argument("AsciiMask", help="The ascii file with the x,y positions ")
    parser.add_argument("-c","--val", type=int, help="the value in counts for the masked pixels. Default = 1",default=1)
 
    args = parser.parse_args()




    ImageFile= args.ImageFile
    AsciiFile= args.AsciiMask
    Value = args.val 

    xy2fits().MakeFits(ImageFile, AsciiFile, Value)

    print("Ascii -> Fits done ")



#console scripts
def maincheckFile() -> None: 

    printWelcome()
    #reading arguments parsing

    parser = argparse.ArgumentParser(description = "checkFile: check that the parameters and file names inside the GALFIT input file are correct")


    # required arguments
    parser.add_argument("GalfitFile", help = "GALFIT input File ")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance in pixels among components. Default = 10", default=10)


    ## parsing variables

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

   
    headinfo, galax, mag, freepar = checkFile(galfitFile, dis)

    testflag = False

    if not(headinfo.inputimageflag):
        line="File {} not found ".format(headinfo.inputimage)
        print(line)
        testflag = True

    #if not(headinfo.outimageflag):
    #    line="File {} not found ".format(headinfo.outimage)
    #    print(line)
    #    testflag = True

    if not(headinfo.sigimageflag):
        line="File {} not found ".format(headinfo.sigimage)
        print(line)
        #testflag = True
        line="GALFIT will create a sigma file "
        print(line)
 
    if not(headinfo.psfimageflag):
        line="File {} not found ".format(headinfo.psfimage)
        print(line)
        testflag = True

    if (headinfo.psfimageflag):
        if not(headinfo.convxflag):
            line="Warning: x-size of convolution is smaller than psf image x-axis"
            print(line)
            testflag = True

        if not(headinfo.convyflag):
            line="Warning: y-size of convolution is smaller than psf image y-axis"
            print(line)
            testflag = True




    if not(headinfo.maskimageflag):
        line="File {} not found ".format(headinfo.maskimage)
        print(line)
        testflag = True

    if not(headinfo.constraintsflag):
        line="File {} not found ".format(headinfo.constraints)
        print(line)
        testflag = True

    if not(headinfo.xsizeflag):
        line="format of size in x-axis not valid. xmax < xmin"
        print(line)
        testflag = True

    if not(headinfo.ysizeflag):
        line="Format of size in y-axis not valid. ymax < ymin"
        print(line)
        testflag = True



    if not(testflag):
        print("No issues with the header information")
    else:
        print("There are some issues with the header information. Check messages above.")

    print('Total number of model components: ', len(galax))
    print('Total number of galaxies: ', len(np.unique(galax)))


    print('Components per galaxy: ')

    cnt = Counter(galax)


    Flux = 10**((25 - mag)/2.5) #25 is just a constant to avoid small numbers. I substract it later

    for idx,item in enumerate(np.unique(galax)): 
        totcomp = cnt[item]  
        maskgal  = galax == item 
         

        totFlux = Flux[maskgal].sum()

        totmag = -2.5*np.log10(totFlux) + 25 

        line="galaxy {} has {} components and a total mag of: {:.2f} ".format(int(item), totcomp, totmag)
        print(line)
 

    line="Total number of free parameters: {} ".format(freepar)
    print(line)
 




def mainGetBoxSizeDs9(): 


    printWelcome()

   #parser argument section

    parser = argparse.ArgumentParser(description="Computes the Box size from a Ds9 region file for galfit header")

    parser.add_argument("RegFile", help="Ds9 region file containing the box region ")
    args = parser.parse_args()


    RegFile = args.RegFile

    xmin, xmax, ymin, ymax = getBoxSizeDs9(RegFile)


    line="xmin, xmax, ymin, ymax: {} {} {} {} ".format(xmin, xmax, ymin, ymax)
    print(line)
 

def maingetSersic(): 

    printWelcome()

    parser = argparse.ArgumentParser(description="prints the Sersic function from DS9 ellipse region")


    parser.add_argument("Image", help = "image fits file")
    parser.add_argument("RegFile", help = "DS9 ellipse region file")


    parser.add_argument("-zp","--zeropoint", type=float, help="The value of the zero point. Default = 25 ",default=25)


    parser.add_argument("-sk","--sky", default=0, type=float, help="Sky background value to be removed from image before photometry. Default = 0")


    parser.add_argument("-bt","--bulgetot", type=float, help="Bulge to total ratio. If set it will print two sersics: one for the bulge and the other for the disk")

    parser.add_argument("-c","--center", action="store_true", help="takes center of ds9 region file")

    parser.add_argument("-n","--noprint", action="store_true", help="avoids to print Sersic functionts to stdout")
    parser.add_argument("-m","--mask", type=str, help="the mask file")

    parser.add_argument("-b","--bards9", type=str, help="DS9 ellipse region file that containts the bar region. bulgetot flag must be activated")

    args = parser.parse_args()

    image = args.Image
    regfile = args.RegFile
    center = args.center
    maskfile = args.mask
    zeropoint = args.zeropoint
    sky = args.sky
    bulgetot = args.bulgetot
    noprint = args.noprint
    bards9 = args.bards9


    galcomps = getSersic(image, regfile, center, maskfile, zeropoint, sky, noprint, bulgetot, bards9)



def main_imarith() -> None: 

    printWelcome()
    #reading arguments parsing

    parser = argparse.ArgumentParser(description = "imarith: makes arithmetics operations on images")

    # required arguments
    parser.add_argument("ImageFile", help = "The input image")


    parser.add_argument("-o", "--output", type=str, default="output.fits", help = "The output image")

    parser.add_argument("-i", "--image2", type=str, help = "second input image to make arithmetic operations with ImageFile. Image2 must be of the same size of ImageFile. If this second image is provided it will make operations indicated by arithmetic flag ignoring its constant input")

    parser.add_argument("-a","--add", type=float, help="add constant to image pixels" )
    parser.add_argument("-d","--div", type=float, help="divide all pixels by constant" )
    parser.add_argument("-m","--mul", type=float, help="multiply all pixels by constant" )
    parser.add_argument("-s","--sub", type=float, help="substract constant to all pixels")



    ## parsing variables

    args = parser.parse_args()

    ImageFile = args.ImageFile
    output = args.output
    image2 = args.image2

    add = args.add
    mul = args.mul
    div = args.div
    sub = args.sub

    imarith(ImageFile, output, image2, add, mul, div, sub)


    print("done")

def mainMakePSF():

    printWelcome()

    parser = argparse.ArgumentParser(description="Makes a PSF model of a star using MGE")

    parser.add_argument("image", help="the image file where it contains the star to be modelled ")
    parser.add_argument("GalfitFile", help="GALFIT file to obtain the header options")

    parser.add_argument("Ds9regFile", help="the DS9 ellipse region file containing the star to model")

    parser.add_argument("-c","--center", action="store_true", 
                        help="uses the center given in DS9 region file, " + 
                        "otherwise it will find the (x,y) peak within DS9 ellipse")

    parser.add_argument("-o","--out", type=str, 
                        help="the PSF model image",default="psf.fits")

    parser.add_argument("-sig","--sigma", type=str, 
                        help="introduce the sigma image")

    parser.add_argument("-t","--twist", action="store_true", help="uses twist option for mge ")

    parser.add_argument("-ng","--numgauss", type=int, help="number of gaussians that will be used for galfit.")

    #######################
    args = parser.parse_args()
    ##

    image = args.image
    regfile = args.Ds9regFile 
    center = args.center

    psfout = args.out

    sigma = args.sigma


    galfitFile = args.GalfitFile
    twist = args.twist

    numgauss = args.numgauss 

    #######################


    makePSF(galfitFile, image, regfile, center, psfout, sigma, twist, numgauss) 



