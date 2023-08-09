

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

Some of the routines are optional arguments
for argparse library and can be ignored. Those 
optional arguments can be empty 
variables (e.g. None or False values but check -h option when calling
the shell command version) 

**GALFIT INPUT**
------------------
Routines that aid the GALFIT's user to
prepare the necessary files for GALFIT input 



**getStar** gets a image slice centered on the object peak

::


    from galfitools.galin.getStar import getStar


    # image:             the image file to obtain the slice
    #regfile:            the DS9 ellipse region file containing the star
    #imsize:   the size of the new image in pixels

    
    #optional for argparse

    #center: boolean flag that indicates to use the center given in 
    #        DS9 region file, otherwise it will find the x,y peaks within DS9 ellipse
 
    #sky:   the value of  the sky background to be removed
    #imout:   the image output name.

    #sigma: the sigma image to obtain the slice if any 
    #sigout: the sigma image output name.


    getStar(image, regfile, imsize, center, sky, imout, sigma, sigout)





**initGal** Creates GALFIT's input files with different initial parameters


::

    from galfitools.galin.initgal import InitGal


    #GalfitFile: the galfit file galfit.XX

    #optional for argparse

    #number:  the number of files generated.
    #param3:  range of values to give to the 3) model's parameter in format [min max]
    #param4: range of values to give to the 4) model's parameter in format [min max]
    #param5: range of values to give to the 5) model's parameter in format [min max]
    #param6: range of values to give to the 6) model's parameter in format [min max]
    #param7: range of values to give to the 7) model's parameter in format [min max]
    #param8: range of values to give to the 8) model's parameter in format [min max]
    #param9: range of values to give to the 9) model's parameter in format [min max]
    #param10: range of values to give to the 10) model's parameter in format [min max] 
    #numcomp: the component number which parameters will be changed


    InitGal(GalfitFile, number, param3, param4, param5, param6, param7, param8, param9, param10, numcomp)


      


**gtmakeMask**  creates mask file from a SExtractor's catalog 

::


    from galfitools.galin.MakeMask import makeMask



    #sexfile:    Sextractor catalog file
    #image: Image file

    #optional for argparse

    #maskfile:  the output mask file name
    #scale: ds9 saturation file
    #satfileout: scale factor to increase the ellipses. Default=1



    makeMask(sexfile, image, maskfile, scale, satfileout)


*Note* The Sextractor catalog must have the following
columns: 

::

    #   1 NUMBER                 Running object number
    #   2 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
    #   3 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
    #   4 X_IMAGE                Object position along x                                    [pixel]
    #   5 Y_IMAGE                Object position along y                                    [pixel]
    #   6 MAG_APER               Fixed aperture magnitude vector                            [mag]
    #   7 KRON_RADIUS            Kron apertures in units of A or B
    #   8 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
    #   9 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
    #  10 A_IMAGE                Profile RMS along major axis                               [pixel]
    #  11 ELLIPTICITY            1 - B_IMAGE/A_IMAGE
    #  12 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
    #  13 BACKGROUND             Background at centroid position                            [count]
    #  14 CLASS_STAR             S/G classifier output
    #  15 FLAGS                  Extraction flags




**maskDs9**  creates (or modify) a mask image for GALFIT using Ds9 regions such as Boxes, Ellipses and Polygons

::


    from galfitools.galin.MaskDs9 import maskDs9

    
    #MaskFile:              the Mask image file to modify or create
    #RegFile:               the DS9 region file

    #optional arguments for argparse
    
    # fill: the value in counts to fill into the Ds9 regions
    #image: image to obtain the size

    #bor_flag:    Mask the borders when their value of this regions is zero
    #borValue:    value of the border if this region has values different from zero 

               

    maskDs9(MaskFile, RegFile, fill, image, bor_flag, borValue) 




**maskSky** creates a mask image for GALFIT using original image and sky mean and sigma

::

    from galfitools.galin.MaskSky import skyRem


    #image:        original data image
    #mask:    Name of the new Mask file

    #optional arguments from argparse

    #sky_mean: mean of the sky background
    #sky_sig:  background
    #nsig:  number of times that the sigma of the sky will be multiplied to remove the
    #        sky background


    #bor_flag:  Mask the borders when their value is zero
    #borValue: value of the border if it is different from zero
                  

    skyRem(image, mask, sky_mean, sky_sig, nsig, borValue, bor_flag)



**xy2fits** code to convert ASCII x,y positions to FTIS mask

::


    from galfitools.galin.xy2fits import xy2fits



    #ImageFile: The Image file
    #AsciiFile: The ascii file with the x,y positions

 

    #optional argument from argparse

    #Value: the value in counts for the masked pixels

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

