
import argparse



from galfitools.sky.GalfitSky import galfitSky

from galfitools.sky.Sky import sky



def mainGalfitSky():

    parser = argparse.ArgumentParser(description="computes the sky using GALFIT")

    parser.add_argument("image", help="the image file")
    parser.add_argument("mask", help="the GALFIT mask file")


    parser.add_argument("-s","--scale", type=float, help="the plate scale. default = 1", default=1)

    parser.add_argument("-zp","--mgzpt", type=float, help="the magnitud zero point. default=25",default = 25)

    parser.add_argument("-x","--xpos", type=float, help="the x position. default=1",default = 1)
    parser.add_argument("-y","--ypos", type=float, help="the y position. default=1",default = 1)

    parser.add_argument("-is","--initsky", type=float, help="the initial sky value default=0",default = 0)


    args = parser.parse_args()

    imgname =  args.image
    maskfile = args.mask 

    mgzpt = args.mgzpt 
    scale = args.scale

    X = args.xpos
    Y = args.ypos

    initsky = args.initsky


    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)



def mainSky():


    parser = argparse.ArgumentParser(description="computes sky from a ds9 region box file")

    parser.add_argument("image", help="the Mask image file to modify or create")
    parser.add_argument("maskfile", help="the Mask image file to modify or create")
    parser.add_argument("Ds9regFile", help="the DS9 ellipse region file containing the galaxy")
 

   
    args = parser.parse_args()

    imgname = args.image 
    maskimage = args.maskfile 
    filereg = args.Ds9regFile


    mean, sig = sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:") 

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))


