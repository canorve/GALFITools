

.. contents::
   :depth: 3
..

-------------------

**API**
============================

To see how to call the routines
for your scripts, check
the python scrips in src/galfitools/shell/commands_*.py


Here we explain everyone how to call
those routines.

Again, here it is divided in 5 sections
as it was done in 
`HowTo <docs/howto.rst>`__

**GALFIT INPUT**
------------------
Routines that aid the GALFIT's user to
prepare the necessary files for GALFIT input 



**getStar** gets a image slice centered on the object peak

::


    from galfitools.galin.getStar import getStar

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





**initGal** Creates GALFIT's input files with different initial parameters


::

    from galfitools.galin.initgal import InitGal

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


      


**gtmakeMask**  creates mask file from a SExtractor's catalog 

::


    from galfitools.galin.MakeMask import makeMask


    args = parser.parse_args()


    sexfile = args.Sexfile 
    image = args.ImageFile  
    maskfile = args.maskout  
    scale = args.scale 
    satfileout = args.satds9





    makeMask(sexfile, image, maskfile, scale, satfileout)
 

**maskDs9**  creates (or modify) a mask image for GALFIT using Ds9 regions such as Boxes, Ellipses and Polygons

::


    from galfitools.galin.MaskDs9 import maskDs9



    args = parser.parse_args()

    MaskFile = args.MaskFile 
    RegFile = args.RegFile 
    fill = args.fill
    image = args.image

    bor_flag = args.border
    borValue = args.borValue



    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue) 




**maskSky** creates a mask image for GALFIT using original image and sky mean and sigma

::

    from galfitools.galin.MaskSky import skyRem

    args = parser.parse_args()

    image = args.ImageFile 
    mask = args.MaskFile
    sky_mean = args.skymean
    sky_sig = args.skysigma
    nsig = args.numbersig

    bor_flag = args.border
    borValue = args.borValue


    skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)



**xy2fits** code to convert ASCII x,y positions to FTIS mask

::


    from galfitools.galin.xy2fits import xy2fits

    args = parser.parse_args()


    ImageFile= args.ImageFile
    AsciiFile= args.AsciiMask
    Value = args.val 

    xy2fits().MakeFits(ImageFile, AsciiFile, Value)





**GALFIT OUTPUT**
-------------------
Routines that computes photometric variables from 
the surface brightness models fitted by GALFIT 


**getBreak** gets the break radius from a set of Sersics

::

      from galfitools.galout.getRads import getBreak
      args = parser.parse_args()

      galfitFile = args.GalfitFile
      dis = args.dis

      eff = args.effrad
      inicomp = args.numinitial

      quick = args.quick
      random = args.random

      num_comp =  args.numcomp

      angle = args.angle


      ranx = args.ranx
      plot = args.plot


      rbreak, N, theta = getBreak(galfitFile, dis, eff, inicomp, quick, random, num_comp, angle, plot, ranx)
          
      print('number of model components: ',N)

      line = 'Using a theta value of: {:.2f} degrees\n'.format(theta)
      print(line)


      line = 'The break radius is {:.2f} pixels \n'.format(rbreak)
      print(line)



**getBreak2** gets the break radius from a set of Sersics using an 
alternative method to getBreak

::

    from galfitools.galout.getRads import getBreak2

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




**getFWHM** gets the FWHM from a set of Sersics
::

    from galfitools.galout.getRads import getFWHM

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





**getKappa** gets the Kappa radius from a set of Sersics

::


    from galfitools.galout.getRads import getKappa

    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    inicomp = args.numinitial

    quick = args.quick
    random = args.random
    angle = args.angle

    ranx = args.ranx
    plot = args.plot


    num_comp =  args.numcomp


    rkappa, N, theta = getKappa(galfitFile, dis, eff, inicomp, quick, random, angle, num_comp, plot, ranx) 

    print('number of model components: ',N)

    line = 'The Kappa radius  is {:.2f} pixels \n'.format(rkappa)
    print(line)





**getReComp** gets the effective radius from a set of Sersics
::


    from galfitools.galout.getRads import getReComp


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





**getSlope** gets the slope radius from a set of Sersics
::

    from galfitools.galout.getRads import getSlope


    args = parser.parse_args()

    galfitFile = args.GalfitFile
    dis = args.dis

    eff = args.effrad
    slope = args.slope

    num_comp =  args.numcomp

    angle = args.angle

    ranx = args.ranx
    plot = args.plot



    rgam, N, theta = getSlope(galfitFile, dis, eff, slope, angle, num_comp, plot, ranx)


    print('number of model components: ',N)

    line = 'Using a theta value of : {:.2f} degrees\n'.format(theta)
    print(line)

    line = 'The radius with slope {:.2f} is {:.2f} pixels \n'.format(slope,rgam)
    print(line)




**getN** computes the Sersic index from surface brightness at effective radius
::


    from galfitools.galout.getN import getN


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

    line = 'Sersic index with the method of Mean Surface Brightness at effective radius: {:.2f}  \n'.format(sersic)
    print(line)


    line = 'Sersic index with the method of fraction of light (at different radius)  \n'
    print(line)


    line = 'Sersic index mean: {:.2f}  Standard deviation: {:.2f}  '.format(meanser, stdser)
    print(line)




**getMissLight** computes the missing light from two surface brightness models
::

    from galfitools.galout.getMissingLight import getMissLight

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





**getBulgeRad** gets the bulge radius or the radius where two models of surface brightness models are
equal
::

    from galfitools.galout.getRads import getBulgeRad

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




**showCube** takes the GALFIT output and creates an image that shows galaxy, model and residual 
::

    from galfitools.galout.showcube import displayCube


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





**photDs9** computes photometry from a Ds9 region file: Box, Ellipses and Polygons
::


    from galfitools.galout.PhotDs9 import photDs9 

    args = parser.parse_args()

    ImageFile = args.ImageFile 
    RegFile = args.RegFile 
    zeropoint = args.zeropoint
    sky = args.sky


    mag = photDs9(ImageFile, RegFile, zeropoint, sky)


    line = "the magnitude from the ds9 region is: {:.2f} \n".format(mag)
    print(line)






**MGE**
---------------

Routines that use the Multi-Gaussian Expansion

**mge2galfit** fits multi-gaussian expansion of Cappellari (2002) and formats to GALFIT
::

    from galfitools.mge.mge2galfit import mge2gal
    args = parser.parse_args()

    mge2gal(args) 


**SbProf** creates a surface brightness profile from a ellipse ds9 region
::


    from galfitools.mge.SbProf import sbProf

    args = parser.parse_args()

    sbProf(args)



 
**SIM**
---------------
Routines that make a simulated galaxy image.

**makeSim** simulates a observed galaxy from a GALFIT model

::


    from galfitools.sim.MakeSim import makeSim

    args = parser.parse_args()

    image = args.image 
    GAIN = args.gain 

    skymean = args.sky
    skystd = args.std 

    newimage = args.newimage 


    makeSim(image, GAIN, skymean, skystd, newimage)


**SKY**
-------------

Routines that compute the sky background

**galSky** computes the sky using GALFIT
::

    from galfitools.sky.GalfitSky import galfitSky

    args = parser.parse_args()

    imgname =  args.image
    maskfile = args.mask 

    mgzpt = args.mgzpt 
    scale = args.scale

    X = args.xpos
    Y = args.ypos

    initsky = args.initsky


    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)



**getSky** computes sky from a ds9 region box file
::


    from galfitools.sky.Sky import sky
   
    args = parser.parse_args()

    imgname = args.image 
    maskimage = args.maskfile 
    filereg = args.Ds9regFile


    mean, sig = sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:") 

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))



--------------

**Questions?**
--------------


Something is not clear for you or do you have further questions?
write to me at canorve [at] gmail [dot] com 

