
import argparse



from galfitools.sky.GalfitSky import galfitSky

from galfitools.sky.Sky import Sky



def mainGalfitSky():

    parser = argparse.ArgumentParser(description="computes the sky using GALFIT")

    parser.add_argument("image", help="the image file")
    parser.add_argument("mask", help="the GALFIT mask file")


    parser.add_argument("-s","--scale", type=float, help="the plate scale. default = 1", default=1)

    parser.add_argument("-zp","--mgzpt", type=float, help="the magnitud zero point. default=25",default = 25)

    parser.add_argument("-x","--xpos", type=float, help="the x position. default=1",default = 1)
    parser.add_argument("-y","--ypos", type=float, help="the y position. default=1",default = 1)



    args = parser.parse_args()

    imgname =  args.image
    maskfile = args.mask 

    mgzpt = args.mgzpt 
    scale = args.scale

    X = args.xpos
    Y = args.ypos




    galfitSky(imgname, maskfile, mgzpt, scale, X, Y)



def mainSky():


    parser = argparse.ArgumentParser(description="computes sky from a ds9 region box file")

    parser.add_argument("image", help="the Mask image file to modify or create")
    parser.add_argument("maskfile", help="the Mask image file to modify or create")
    parser.add_argument("Ds9regFile", help="the DS9 ellipse region file containing the galaxy")
 

    if len(sys.argv[1:]) != 2 and len(sys.argv[1:]) != 6 and len(sys.argv[1:]) != 3:
        print ('Missing arguments')
        print ("Usage:\n %s [ImageFile] [MaskImage] [optional Ds9RegFile] [optional: Xmin Xmax Ymin Ymax]" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits ds9.reg" % sys.argv[0])
        print ("Example:\n %s image.fits mask.fits 330 450 200 700" % sys.argv[0])

        sys.exit()


    flagseg=False
    flagreg=False
    xmin=xmax=ymin=ymax=0
    if len(sys.argv[1:]) == 6:
        flagseg=True
        xmin=int(sys.argv[3])
        xmax=int(sys.argv[4])
        ymin=int(sys.argv[5])
        ymax=int(sys.argv[6])

    if len(sys.argv[1:]) == 3:
        flagreg=True
        filereg=sys.argv[3]

   
    args = parser.parse_args()

    imgname = args.image 
    maskimage = args.maskfile 
    filereg = args.Ds9regFile


    mean, sig = Sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:") 

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))


