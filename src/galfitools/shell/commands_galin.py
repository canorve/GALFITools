

import argparse

from galfitools.galin.getStar import getStar
from galfitools.galin.initgal import InitGal
from galfitools.galin.MaskDs9 import maskDs9

from galfitools.galin.MakeMask import makeMask
from galfitools.galin.MaskSky import skyRem
from galfitools.galin.xy2fits import xy2fits

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
                        help="the sky background to be removed. Default = 0",default=0)
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




    args = parser.parse_args()

    MaskFile = args.MaskFile 
    RegFile = args.RegFile 
    fill = args.fill
    image = args.image

    bor_flag = args.border
    borValue = args.borValue



    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue) 

    print('Done. Mask {} created (or modified)'.format(MaskFile))



def mainMaskSky(): 


    printWelcome()

    parser = argparse.ArgumentParser(description="creates a mask image for GALFIT using original image and sky mean and sigma")


    parser.add_argument("ImageFile", help="original data image ")
    parser.add_argument("MaskFile", help="Name of the new Mask file")
    parser.add_argument("-sm","--skymean",default=0,type=float, help="mean of the sky background")
    parser.add_argument("-ss","--skysigma",default=0,type=float, help="sigma of the sky background")
    parser.add_argument("-ns","--numbersig",default=1,type=float, help="number of times that the sigma of the sky will be multiplied to remove the sky background")

    parser.add_argument("-b","--border", action="store_true", help="Mask the borders when their value is zero")

    parser.add_argument("-bv","--borValue",default=0,type=float, help="value of the border if it is different from zero")
    args = parser.parse_args()

    image = args.ImageFile 
    mask = args.MaskFile
    sky_mean = args.skymean
    sky_sig = args.skysigma
    nsig = args.numbersig

    bor_flag = args.border
    borValue = args.borValue


    skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)

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



