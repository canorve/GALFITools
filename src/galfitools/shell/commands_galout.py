

import argparse

from galfitools.galout.getRads import getBreak
from galfitools.galout.getRads import getBreak2
from galfitools.galout.getRads import getFWHM
from galfitools.galout.getRads import getKappa
from galfitools.galout.getRads import getReComp
from galfitools.galout.getRads import getSlope
from galfitools.galout.getRads import getBulgeRad
from galfitools.galout.getMissingLight import getMissLight
from galfitools.galout.getN import getN
from galfitools.galout.showcube import displayCube

from galfitools.galout.PhotDs9 import photDs9 


from galfitools.shell.prt import printWelcome


def mainPhotDs9(): 

    printWelcome()

    parser = argparse.ArgumentParser(description="computes photometry from a Ds9 region file: Box, Ellipses and Polygons ")

    parser.add_argument("ImageFile", help="the image file where the photometry will be computed")
    parser.add_argument("RegFile", help="the DS9 region file")

    parser.add_argument("-zp","--zeropoint", type=float, help="The value of the zero point. Default = 25 ",default=25)


    parser.add_argument("-sk","--sky", default=0, type=float, help="the value of the sky background to be removed")




    args = parser.parse_args()

    ImageFile = args.ImageFile 
    RegFile = args.RegFile 
    zeropoint = args.zeropoint
    sky = args.sky


    mag = photDs9(ImageFile, RegFile, zeropoint, sky)


    line = "the magnitude from the ds9 region is: {:.2f} \n".format(mag)
    print(line)



#console scripts
def mainGetBreak() -> None: 

    printWelcome()
    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getBreak: gets the break radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-ni","--numinitial", type=int, help="Number of component where it'll obtain the initial parameter to search break radius or to generated random initial radius. ", default=2)


    parser.add_argument("-q","--quick", action="store_true", help='evaluate in position only (given by -ni parameter') 



    parser.add_argument("-r","--random", type=int, help="Number of random radius as initial parameters to search for the minimum. It will generated random radius from 0 to effective radius of the component indicated by parameter -ni ")


    parser.add_argument("-p","--plot", action="store_true", help='makes plot of double derivative vs. radius ') 
    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="provide a range for the plot x-axis: xmin - xmax ")




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    inicomp = args.numinitial

    quick = args.quick
    random = args.random

    num_comp =  args.numcomp

    angle = args.angle


    ranx = args.ranx
    plot = args.plot


    rbreak, N, theta = getBreak(galfitFile, dis, inicomp, quick, random, num_comp, angle, plot, ranx)
        
    print('number of model components: ',N)

    line = 'Using a theta value of: {:.2f} degrees\n'.format(theta)
    print(line)


    line = 'The break radius is {:.2f} pixels \n'.format(rbreak)
    print(line)



def mainGetBreak2():
    '''gets the break radius from a set of Sersics'''


    printWelcome()

    parser = argparse.ArgumentParser(description = "getBreak2: gets the break radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)
    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-p","--plot", action="store_true", help='makes plot of double derivative vs. radius ') 

    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="x-axis range to search for the Break radius: xmin - xmax ")


    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis
    angle = args.angle
    num_comp =  args.numcomp
    plot = args.plot
    ranx = args.ranx

    rbreak, N, theta =  getBreak2(galfitFile, dis, angle, num_comp, plot, ranx)


    print('number of model components: ',N)

    line = 'Using a theta value of: {:.2f} degrees\n'.format(theta)
    print(line)


    line = 'The break radius is {:.2f} pixels \n'.format(rbreak)
    print(line)




def mainFWHM() -> None: 


    printWelcome()

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

    printWelcome()
    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getKappa: gets the Kappa radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-ni","--numinitial", type=int, help="Number of component where it'll obtain the initial parameter to search break radius or to generated random initial radius. ", default=2)


    parser.add_argument("-q","--quick", action="store_true", help='evaluate in position only (given by -ni parameter') 


    parser.add_argument("-r","--random", type=int, help="Number of random radius as initial parameters to search for the minimum. It will generated random radius from 0 to effective radius of the component indicated by parameter -ni ")


    parser.add_argument("-p","--plot", action="store_true", help='makes plot of maximum curvature vs. radius ') 
    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="provide a range for x-axis: xmin - xmax ")




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    inicomp = args.numinitial

    quick = args.quick
    random = args.random
    angle = args.angle

    ranx = args.ranx
    plot = args.plot


    num_comp =  args.numcomp


    rkappa, N, theta = getKappa(galfitFile, dis, inicomp, quick, random, angle, num_comp, plot, ranx) 

    print('number of model components: ',N)

    line = 'The Kappa radius  is {:.2f} pixels \n'.format(rkappa)
    print(line)


def mainGetReComp() -> None: 

    printWelcome()
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

    EffRad, totmag, meanme, me, N, theta = getReComp(galfitFile, dis, eff, angle, num_comp)

    print('number of model components: ', N)

    line = 'Using a theta value of : {:.2f} degrees \n'.format(theta)
    print(line)

    line = 'Total Magnitude of the galaxy: {:.2f} \n'.format(totmag)
    print(line)

    line = 'Surface brightness at effective radius (\u03BCe): {:.2f} mag/\" \n'.format(me)
    print(line)


    line = 'Mean Surface Brightness at effective radius (<\u03BC>e): {:.2f} mag/\" \n'.format(meanme)
    print(line)

    line = 'The radius at {:.0f}% of light is {:.2f} pixels \n'.format(eff*100,EffRad)
    print(line)



def maingetSlope() -> None: 

    printWelcome()
    #reading argument parsing

    parser = argparse.ArgumentParser(description = "getSlope: gets the slope radius from a set of Sersics ")


    # required arguments
    parser.add_argument("GalfitFile", help = "Galfit File containing the Sersics or gaussians components")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-a","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components")


    parser.add_argument("-s","--slope", type=float, 
                        help="value of slope to find. default=.5 ", default=.5)



    parser.add_argument("-p","--plot", action="store_true", help='makes plot of first derivative vs. radius ') 
    parser.add_argument("-rx","--ranx",nargs=2, type=float, help="provide a range for x-axis: xmin - xmax ")




    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    slope = args.slope

    num_comp =  args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot




    rgam, N, theta = getSlope(galfitFile, dis, slope, angle, num_comp, plot, ranx)


    print('number of model components: ',N)

    line = 'Using a theta value of : {:.2f} degrees\n'.format(theta)
    print(line)

    line = 'The radius with slope {:.2f} is {:.2f} pixels \n'.format(slope,rgam)
    print(line)



#check modify
#console scripts
def maingetN() -> None: 

    printWelcome()
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



    parser.add_argument("-p","--plot", action="store_true", help='makes plot of Sersic index vs. fraction of light ') 




    ## parsing variables

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis


   
    num_comp =  args.numcomp

    frac = args.radfrac

    angle = args.angle

    plot = args.plot


    sersic, meanser, stdser, totmag, N, theta = getN(galfitFile, dis, frac, angle, num_comp, plot)


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



    line = 'Sersic index mean: {:.2f}  Standard deviation: {:.2f}  '.format(meanser, stdser)
    print(line)




def mainGetBulgeRad() -> None: 
    '''gets the bulge radius or the radius where two models of surface brightness models are equal'''

    printWelcome()
    #reading arguments parsing

    parser = argparse.ArgumentParser(description = "getBulgeRad: gets the bulge radius or the radius where two models of surface brightness models are equal")


    # required arguments
    parser.add_argument("GalfitFile1", help = "Galfit File containing the coreless surface brightness model")
    parser.add_argument("GalfitFile2", help = "Galfit File containing the core surface brightness model")

    #parser.add_argument("gammarad", type=float, help="gamma radius")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)

    #parser.add_argument("-ser","--sersic", action="store_true", help="uses sersic function for galfit file")

    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)

    parser.add_argument("-pa","--angle", type=float, 
                        help="Angle of the major axis of the galaxy. Default= it will take the angle of the last components. Angle measured from Y-Axis as same as GALFIT. ")


    parser.add_argument("-p","--plot", action="store_true", help='makes plot of GalfitFile1 - GalfitFile2 vs. radius ') 
    parser.add_argument("-rx","--ranx", nargs=2, type=float, help="provide a range for x-axis: xmin - xmax ")


    ## parsing variables

    args = parser.parse_args()

    galfitFile1 = args.GalfitFile1
    galfitFile2 = args.GalfitFile2

    dis = args.dis

    num_comp =  args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot


    rbulge, N1, N2, theta = getBulgeRad(galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx)


    print('number of model components for the bulge: ',N1)
    print('number of model components for the rest of the galaxy: ',N2)


    line = 'Using a theta value of: {:.2f} degrees\n'.format(theta)
    print(line)


    line = 'The bulge radius is {:.2f} pixels \n'.format(rbulge)
    print(line)




def mainMissingLight() -> None: 
    '''gets the missing light from a two models'''

    printWelcome()
    #reading arguments parsing

    parser = argparse.ArgumentParser(description = "getMissLight: computes the missing light from two surface brightness models")


    # required arguments
    parser.add_argument("GalfitFile1", help = "Galfit File containing the coreless surface brightness model")
    parser.add_argument("GalfitFile2", help = "Galfit File containing the core surface brightness model")

    parser.add_argument("rad", type=float, help="upper limit of radius to integrate the missing light in pixels")

    parser.add_argument("-d","--dis", type=int, help="Maximum distance among components", default=10)


    parser.add_argument("-n","--numcomp", type=int, help="Number of component where it'll obtain center of all components, default = 1 ", default=1)


    ## parsing variables

    args = parser.parse_args()

    galfitFile1 = args.GalfitFile1
    galfitFile2 = args.GalfitFile2

    dis = args.dis

    num_comp =  args.numcomp

    rad = args.rad



    magmiss, N1, N2 = getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad)


    print('number of model components coreless model: ',N1)
    print('number of model components core model: ',N2)


    line = 'the missing light is {:.2f} mag \n'.format(magmiss)
    print(line)








def mainShowCube():


    printWelcome()
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



