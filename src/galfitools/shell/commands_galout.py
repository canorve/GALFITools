

import argparse

from galfitools.galout.getRads import getBreak
from galfitools.galout.getRads import getFWHM
from galfitools.galout.getRads import getKappa
from galfitools.galout.getRads import getReComp
from galfitools.galout.getRads import getSlope
from galfitools.galout.getRads import getSlope
from galfitools.galout.getN import getN
from galfitools.galout.showcube import displayCube


#console scripts
def mainGetBreak() -> None: 

    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getBreak: gets the break radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-er","--effrad", type=float, 
                        help="percentage of light to compute for radius. default=.5 for effective radius ", default=.5)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-ni","--numinitial", type=int, help="Number of component where it'll obtain the initial parameter to search break radius or to generated random initial radius. ", default=2)


    parser.add_argument("-q","--quick", action="store_true", help='evaluate in position only (given by -ni parameter') 


    parser.add_argument("-r","--random", type=int, help="Number of random radius as initial parameters to search for the minimum. It will generated random radius from 0 to effective radius of the component indicated by parameter -ni ")




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    inicomp = args.numinitial

    quick = args.quick
    random = args.random

    num_comp =  args.numcomp

    angle = args.angle


    rbreak, N, theta = getBreak(galfitfile, dis, eff, inicomp, quick, random, num_comp, angle)
        
    print('number of model components: ',N)

    line = 'Using a theta value of: {:.2f} degrees\n'.format(theta)
    print(line)


    line = 'The break radius is {:.2f} pixels \n'.format(rbreak)
    print(line)



def mainFWHM() -> None: 

    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getFWHM: gets the FWHM from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")





    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    angle = args.angle


   
    num_comp =  args.numcomp

    fwhm, N, theta = getFWHM(galfitFile, dis, angle, num_comp)


    print('number of model components: ',N)

    line = 'Using a theta value of : {:.2f} degrees\n'.format(theta)
    print(line)

    line = 'The FWHM is {:.2f} pixels \n'.format(fwhm)
    print(line)


#console scripts
def mainKappa() -> None: 

    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getKappa: gets the Kappa radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-er","--effrad", type=float, 
                        help="percentage of light to compute for radius. default=.5 for effective radius ", default=.5)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-ni","--numinitial", type=int, help="Number of component where it'll obtain the initial parameter to search break radius or to generated random initial radius. ", default=2)


    parser.add_argument("-q","--quick", action="store_true", help='evaluate in position only (given by -ni parameter') 


    parser.add_argument("-r","--random", type=int, help="Number of random radius as initial parameters to search for the minimum. It will generated random radius from 0 to effective radius of the component indicated by parameter -ni ")


    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    inicomp = args.numinitial

    quick = args.quick
    random = args.random
    angle = args.angle


    num_comp =  args.numcomp


    rkappa, N, theta = getKappa(galfitFile, dis, eff, inicomp, quick, random, angle, num_comp) 

    print('number of model components: ',N)

    line = 'The Kappa radius  is {:.2f} pixels \n'.format(rkappa)
    print(line)


def mainGetReComp() -> None: 

    #reading arguments parsing
    parser = argparse.ArgumentParser(description = "getReComp: gets the effective radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-er","--effrad", type=float, 
                        help="percentage of light to compute for radius. default=.5 for effective radius ", default=.5)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-pa","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components. Angle measured from Y-Axis as same as GALFIT. ")



    args = parser.parse_args()


    galfitFile = args.GalfitFile
    dis = args.dis
    eff = args.effrad
    num_comp =  args.numcomp
    angle = args.angle

    Effrad, totmag, meanme, me, N, theta = getReComp(galfitFile, dis, eff, angle, num_comp)

    print('number of model components: ', N)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta)
    print(line)

    line = 'Total Magnitude of the galaxy: {:.2f} \n'.format(totmag)
    print(line)

    line = 'Mean Surface Brightness at effective radius: {:.2f} mag/\" \n'.format(meanme)
    print(line)

    line = 'Surface brightness at effective radius {:.2f} mag/\" \n'.format(me)
    print(line)

    line = 'The radius at {:.0f}% of light is {:.2f} pixels \n'.format(eff*100,EffRad)
    print(line)



def maingetSlope() -> None: 

    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getSlope: gets the slope radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-er","--effrad", type=float, 
                        help="percentage of light to compute for radius. default=.5 for effective radius ", default=.5)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-s","--slope", type=float, 
                        help="value of slope to find. default=.5 ", default=.5)




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    slope = args.slope

    num_comp =  args.numcomp

    angle = args.angle
   
    rgam, N, theta = getSlope(galfitFile, dis, eff, slope, angle, num_comp)


    print('number of model components: ',N)

    line = 'Using a theta value of : {:.2f} degrees\n'.format(theta)
    print(line)

    line = 'The radius with slope {:.2f} is {:.2f} pixels \n'.format(slope,rgam)
    print(line)



#check modify
#console scripts
def maingetN() -> None: 

    #reading arguments parsing

    parser = argparse.ArgumentParser(description = "getN: computes the Sersic index from surface brightness at effective radius")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-pa","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components. Angle measured from Y-Axis as same as GALFIT. ")


    parser.add_argument("-rf","--radfrac", type=float, help="fraction of light radius. Default = .2 ", default=.2)

    ## parsing variables

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis


   
    num_comp =  args.numcomp

    frac = args.radfrac


    angle = args.angle


    sersic, meanser, stdser, totmag, N, theta = getN(galfitFile, dis, frac, angle, num_comp)


    print('number of model components: ', N)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta)
    print(line)

    line = 'Total Magnitude of the galaxy: {:.2f} \n'.format(totmag)
    print(line)

    #line = 'Mean Surface Brightness at effective radius: {:.2f} mag/\" \n'.format(meanme)
    #print(line)

    #line = 'Surface brightness at effective radius {:.2f} mag/\" \n'.format(me)
    #print(line)

    #line = 'The radius at {:.0f}% of light is {:.2f} pixels \n'.format(eff*100,EffRad)
    #print(line)

    line = 'Sersic index with the method of Mean Surface Brightness at effective radius: {:.2f}  \n'.format(sersic)
    print(line)


    line = 'Sersic index with the method of fraction of light (at different radius)  \n'
    print(line)



    line = '\nSersic index mean: {:.2f}  Standard deviation: {:.2f}  '.format(meanser, stdser))
    print(line)






def mainShowCube():


    parser = argparse.ArgumentParser(description="show the cube fits of the galfit output")

    parser.add_argument("cubeimage", help="the cube GALFIT image")
    parser.add_argument("-o","--outimage", type=str, help="the output png file",default="cube.png")


    #### 
    parser.add_argument("-br","--brightness", type=float, 
                        help="brightness of the image. Only for galaxy and model. Default = 0. Preferible range goes from -1 to 1", default=0)
    parser.add_argument("-co","--contrast", type=float, 
                        help="contrast of the image. Only for galaxy and model. Default = 1. Preferible range goes from 0 to 1",default=1)

    parser.add_argument("-cm","--cmap", type=str, help="cmap to be used for the cube image ",default="viridis")
    parser.add_argument("-dpi","--dotsinch", type=int, help="dots per inch used for images files ",default=100)
    parser.add_argument("-s","--scale", type=float, help="plate scale of the image. Default = 1",default=1)
    parser.add_argument("-np","--noplot", action="store_true", help="it doesn\'t show plotting window")

    args = parser.parse_args()

    cubeimage = args.cubeimage
    namecube = args.outimage 
    dpival = args.dotsinch 
    brightness = args.brightness
    contrast = args.contrast 
    cmap = args.cmap 
    scale = args.scale  
    noplot = args.noplot


    displayCube(cubeimage, namecube, dpival, brightness, contrast, cmap, scale, noplot)



