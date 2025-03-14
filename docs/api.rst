.. _api:

**API Reference**
====================



.. contents::
   :depth: 3
..


.. image:: https://img.shields.io/pypi/v/GALFITools.svg
    :alt: PyPI-Server
    :target: https://pypi.org/project/GALFITools/

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8248042.svg
   :target: https://doi.org/10.5281/zenodo.8248042





--------------------------------

**API Integration Guidelines**
-------------------------------

The integration of GALFITools routines 
within your custom scripts is facilitated through 
the utilization of Application Programming Interfaces (APIs). 
These APIs offer a systematic and efficient means of 
invoking the routines for specific tasks. 

To comprehensively comprehend the invocational 
mechanics, it is advisable to reference the Python 
scripts located in the directory "src/galfitools/shell/commands_*.py". 
These scripts serve as examples, illustrating 
the precise syntax and methodologies for calling 
the respective routines. Analyzing these scripts 
aids in the assimilation of the necessary syntax 
and parameters for API invocation.


The API integration process is explained across 
five distinct sections, mirroring the structural 
delineation observed in the  :doc:`usage` section. 

Some of the routines are optional arguments
for argparse library and can be ignored. Those 
optional arguments can be empty 
variables (e.g. None or False values but check -h option when calling
the shell command version).


While these optional arguments 
can be omitted, it is recommended 
to consult the help documentation ("-h" option) 
associated with the respective shell command version. 




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

**checkFile** check that the parameters and file names inside the GALFIT input file are correct 

::



    from galfitools.galin.checkGalFile import checkFile 

    #galfitFile is the galfit input parameter file
    #dis is the maximum distance among components of the same galaxy

    headinfo, galax, mag = checkFile(galfitFile, dis)

    #output:

    # galax: is an array with a size of the number of components. It indicates the
    # galaxy number which belongs to the galaxy. This has the same order as 
    #the galfit input file

    #mag: is the magnitud of every component

    #headinfo is a class that contains the name of the files 
    #which comes in the galfit header file. It contains a flag that indicates
    # if the file exists (True) or not (False). Check galhead class below

**Galfit**, **galfit.ReadHead**, **galfit.ReadComps** and  **galfit.ReadSky**. Class 
functions to read the galfit output file galfit.XX. The class functions return a data
class with the parameters read from sky. 


::

    from galfitools.galin.galfit import Galfit



    galfit = Galfit(galfitFile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()
 


**galhead** is a data class that stores the variables of the header of the galfit file:

::


    from galfitools.galin.galfit import GalHead


    class GalHead():
        '''store the header of galfit file'''

        inputimage = "none.fits"     # Input data image (FITS file)
        outimage = "none-out.fits"   # Output data image block
        sigimage = "none"            # Sigma image name (made from data if blank or "none") 
        psfimage = "none"            # Input PSF image and (optional) diffusion kernel
        psfsamp = 1                  # PSF fine sampling factor relative to data 
        maskimage = "none"           # Bad pixel mask (FITS image or ASCII coord list)
        constraints = "none"         # File with parameter constraints (ASCII file) 
        xmin = 0                     # Image region to fit (xmin)
        xmax = 1                     # Image region to fit (xmax)
        ymin = 0                     # Image region to fit (ymin)
        ymax = 1                     # Image region to fit (ymax)
        convx = 1                    # Size of the convolution box (x)
        convy = 1                    # Size of the convolution box (y)
        mgzpt = 25                   # Magnitude photometric zeropoint
        scale = 1                    # Plate scale (dx)   [arcsec per pixel]
        scaley = 1                   # Plate scale (dy)   [arcsec per pixel]
        display = "regular"          # Display type (regular, curses, both)
        P = 0                        # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

        # internal variables of the data class

        imgidx = "sci"
        flagidx = False
        num = 1
        flagnum = False
        exptime = 1
        tempmask = "tempmask.fits"


**galcomps** is a data class that stores the variables of every component of the galfit file:

::


    from galfitools.galin.galfit import GalComps

    #lastmod
    class GalComps:
        '''stores the components of galfit file'''

        #all the variables are stored as an array.

        N = np.array([])               #   number of the component
        NameComp = np.array([])        #0) Name of the component
        PosX = np.array([])            #1) X - position in pixels   
        PosY = np.array([])            #2) Y - position in pixels 
        Mag = np.array([])             #3) magnitud of the component
        Rad = np.array([])             #4) Radius. If Sersic this is Re, and so on for every model
        Exp = np.array([])             #5) Exponent. If Sersic this is for 
        Exp2 = np.array([])            #6) exponent for moffat
        Exp3 = np.array([])            #7) exponnent for moffat
                                       #8) There is No 8 in any galfit model
        AxRat = np.array([])           #9) Axis ratio of the component
        PosAng = np.array([])          #10) angular position of the component
        skip = np.array([])            #z) skip model from output

        Active = np.array([])          # For simultaneous fitting, this paramters tells 
                                       # which components belong to the galaxy of interest
                                       # Activate = True/False 

        # The params below correspond to the variables above and
        # tells to GALFIT  whether this parameter must keep fixed during the fitting 
        PosXFree = np.array([])            #1)   
        PosYFree = np.array([])            #2)   
        MagFree = np.array([])             #3)
        RadFree = np.array([])             #4)
        ExpFree = np.array([])             #5)
        Exp2Free = np.array([])            #6)  for moffat
        Exp3Free = np.array([])            #7)  for moffat
                                           #8)  There is No 8 in any galfit model
        AxRatFree = np.array([])           #9)  AxisRatio
        PosAngFree = np.array([])          #10) position angle

        # the parameters below are not from galfit file, but computed with
        # the routines of this library
        Rad50 = np.array([])            # Radius that keeps the 50% of light
        SerInd = np.array([])           # Computed Sersic index 
        Rad50kpc = np.array([])         # Radius that keeps the 50% of light in kpc
        Rad50sec = np.array([])         # Radius that keeps the 50% of light in arc sec
        Rad90 = np.array([])            # Radius that keeps the 90% of light
        AbsMag = np.array([])           # absolute magnitude
        Lum = np.array([])              # Luminosity
        Flux = np.array([])             # Flux
        PerLight = np.array([])         # Percentage of light that have this component with respect to galaxy
        me = np.array([])               # surface brightness at Re
        mme = np.array([])              # mean surface brightness at Re
        kser = np.array([])             # K parameter related to n to keep Ie at Re

        KronRad = np.array([])          # kron radius
        PetRad = np.array([])           # petrosian radius




**galsky** is a data class that stores the variables of the sky component of the galfit file:

::


    from galfitools.galin.galfit import GalSky

    class GalSky:
        '''stores the value of the GALFIT file'''

        sky = 0  #sky background
        dskyx = 0 # sky gradient in x
        dskyy = 0 #sky gradient in y
        skip = 0 #skip component from model output 

        skyfree = 1  #keep varying this parameter for sky background
        dskyxfree = 0  #keep varying this parameter for sky gradient in x
        dskyyfree = 0 #keep varying this parameter for sky  gradient in y

     





**conver2Sersic**
::

    from galfitools.galin.galfit import conver2Sersic


**SelectGal**
::

    from galfitools.galin.galfit import  SelectGal



**numComps**
::

    from galfitools.galin.galfit import numComps


**GetRadAnd**
::

    from galfitools.galin.galfit import GetRadAng


**getBoxSizeDs9**

::

    from galfitools.galin.getBoxSizeDs9 import getBoxSizeDs9 


    xmin, xmax, ymin, ymax = getBoxSizeDs9(RegFile)

    #RegFile: ds9 box region file
    #xmin, xmax, ymin, ymax box size for fitting region for galfit file option H) 

    
**imarith** makes arithmetic operations on images 

::


    from galfitools.galin.imarith import imarith 


    #Imagefile  input image file
    #output  output image file
    #image2 second optional input image to perform arithmetic operations

    #add add constant to all pixels on the image
    #mul multiply constant to all pixels on the image
    #div divide constant to all pixels on the image 
    #sub substract constant to all pixels on the image 

    #please perform one arithmetic operation per call to imarith

    imarith(ImageFile, output, image2, add, mul, div, sub)


**getSersic** Its estimates and prints initial parameters for Sersic components. It
              addtion if proved options for single Sersic, bulge/disk and bulge/bar/disk

::

    from galfitools.galin.getSersic import getSersic


    #image: fits  image of the galaxy
    #regfile: Ds9 ellipse region containing the galaxy
    #center: if activated, it then uses the center of the Ds9 ellipse region
    #maskfile: fits mask  file (the same file provided for GALFIT)
    #zeropoint: magnitude zero point
    #sky: value of the sky background in pixel units 
    #bulgetot: if provided it divides the magnitud in two components: bulge and disk according
    # to the value of bulgetot (value expected to be between 0 and 1)
    #noprint: avoids to print to stdout and just returns the galcomps data class
    #bards9: if provided it used the info of the ds9 ellipse region to estimate the
    # initial parameter of the bar. Note: This is a different file of the one provided in regfile

    #galcomps: data class containing the initial parameters of every component

    galcomps = getSersic(image, regfile, center, maskfile, zeropoint, sky, noprint, bulgetot, bards9)


**makePSF** Makes a PSF fits model from a star using Multi Gaussian Expansion (MGE) 
::

    # galfitFile: GALFIT file. Used to read header and sky component value
    # image: Fits image containing the galaxy 
    # regfile: Ds9 ellipse region file containing the star
    # center: If activated, it used the center of the Ds9 region instead of 
    # the peak (default mode)
    # psfout: name of the output PSF fits model
    # sigma: Sigma image used for galfit (if any)
    # twist: use twist mode for MGE 

    numgauss: maximum number of gaussians used for MGE



    makePSF(galfitFile, image, regfile, center, psfout, sigma, twist, numgauss) 




**GALFIT OUTPUT**
-------------------
Routines that computes photometric variables from 
the surface brightness models fitted by GALFIT 


**getBreak** gets the break radius from a set of Sersics

::

      from galfitools.galout.getRads import getBreak

      #galfitFile: Galfit File containing the Sersics or gaussians components

      #optional from argparse:
      #dis: Maximum distance among components

      #inicomp: Number of component where it'll obtain the initial parameter to search break
      #                  radius or to generated random initial radius.

      #quick: evaluate in the position given by inicomp parameter

      #random: Number of random radius as initial parameters to search for the minimum. It
      #        will generated random radius from 0 to effective radius of the component
      #        indicated by parameter -ni
      
      #num_comp: Number of component where it'll obtain center of all components, default = 1
      #angle:  Angle of the major axis of the galaxy measured from the image Y-axis
      #ranx: list that indicates the range for the plot x-axis: xmin - xmax
      #plot: boolean flag that indicates to  make a plot of double derivative vs. radius



      rbreak, N, theta = getBreak(galfitFile, dis, eff, inicomp, quick, random, num_comp, angle, plot, ranx)

      # output variables:

      #rbreak: the break radius in pixels  
      #N: number of surface brightness model components of the galaxy
      #theta: the angle used to determine the break radius. Break radius
      #  is computed in that angle direction.


**getBreak2** gets the break radius from a set of Sersics using an 
alternative method to getBreak.

::

    from galfitools.galout.getRads import getBreak2


    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:
    #dis: Maximum distance among components
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis
    #num_comp: Number of component where it'll obtain center of all components, default = 1
    #plot: boolean flag that indicates to  make a plot of double derivative vs. radius
    #ranx: list that indicates the range for the plot x-axis: xmin - xmax

    rbreak, N, theta =  getBreak2(galfitFile, dis, angle, num_comp, plot, ranx)

    # output variables:

    #rbreak: the break radius in pixels  
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the break radius. Break radius
    #  is computed in that angle orientation



**getFWHM** gets the FWHM from a set of Sersics
::

    from galfitools.galout.getRads import getFWHM


    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:

    #dis: Maximum distance among components
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis 
    #num_comp: Number of component where it'll obtain center of all components, default = 1



    fwhm, N, theta = getFWHM(galfitFile, dis, angle, num_comp)

    # output variables:

    #fwhm: the fwhm 
    #N: number of surface brightness model components of the star
    #theta: the angle used to determine the FWHM. it 
    #  is computed in that angle orientation





**getKappa** gets the Kappa radius from a set of Sersics

::


    from galfitools.galout.getRads import getKappa

    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:
    #dis: Maximum distance among components

    #inicomp: Number of component where it'll obtain the initial parameter to search break
    #                  radius or to generated random initial radius.

    #quick: evaluate in the position given by inicomp parameter

    #random: Number of random radius as initial parameters to search for the minimum. It
    #        will generated random radius from 0 to effective radius of the component
    #        indicated by parameter -ni
      
    #num_comp: Number of component where it'll obtain center of all components, default = 1
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis 
    #ranx: list that indicates the range for the plot x-axis: xmin - xmax
    #plot: boolean flag that indicates to  make a plot of maximum curvature vs. radius




    rkappa, N, theta = getKappa(galfitFile, dis, eff, inicomp, quick, random, angle, num_comp, plot, ranx) 


    # output variables:

    #rkappa: the kappa radius in pixels  
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the kappa radius. It 
    #  is computed in that angle orientation


**getKappa2** gets the kappa radius from a set of Sersics using an 
alternative method to getKappa.

::

    from galfitools.galout.getRads import getBreak2


    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:
    #dis: Maximum distance among components
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis
    #num_comp: Number of component where it'll obtain center of all components, default = 1
    #plot: boolean flag that indicates to  make a plot of double derivative vs. radius
    #ranx: list that indicates the range for the plot x-axis: xmin - xmax

    rkappa, N, theta =  getKappa2(galfitFile, dis, angle, num_comp, plot, ranx)

    # output variables:

    #rkappa: the kappa radius in pixels  
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the break radius. Break radius
    #  is computed in that angle orientation






**getReComp** gets the effective radius from a set of Sersics
::


    from galfitools.galout.getRads import getReComp


    #galfitFile = Galfit File containing the Sersics or gaussians components
    #dis: Maximum distance among components
    #eff: percentage of light of the radius to be computed. Effective radius = 0.5  
    #num_comp:Number of component where it'll obtain center of all components, default = 1
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis 

    EffRad, totmag, meanme, me, N, theta = getReComp(galfitFile, dis, eff, angle, num_comp)


    # output variables:

    #EffRad: the computed fraction of light radius  in pixels  
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the kappa radius. It 
    #  is computed in that angle orientation
    #totmag: total magnitude of the galaxy.
    #me: Surface brightness at effective radius
    #mme: Mean surface brightness at effective radius



**getSlope** gets the slope radius from a set of Sersics
::

    from galfitools.galout.getRads import getSlope


    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:
    #dis: Maximum distance among components

    #slope = value of slope at which the radius is to be found. 

    #num_comp: Number of component where it'll obtain center of all components, default = 1

    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis 

    #ranx: list that indicates the range for the plot x-axis: xmin - xmax
    #plot: boolean flag that indicates to make a plot of first derivative vs. radius



    rgam, N, theta = getSlope(galfitFile, dis, eff, slope, angle, num_comp, plot, ranx)


    # output 

    #rgam: the pixel radius at which the model have the specified slope value 
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the break radius. Break radius
    #  is computed in that angle direction.




**getN** computes the Sersic index from surface brightness at effective radius
::


    from galfitools.galout.getN import getN



    #galfitFile: Galfit File containing the Sersics or gaussians components

    #optional from argparse:
    #dis: Maximum distance among components
    #num_comp: Number of component where it'll obtain center of all components, default = 1

   
    #frac: fraction of light radius 
    #angle:  Angle of the major axis of the galaxy measured from the image Y-axis 
    #plot: boolean flag that indicates  to make plot of Sersic index vs. fraction of light
    #const: constant to be substracted from plot


    sersic, meanser, stdser, totmag, N, theta = getN(galfitFile, dis, frac, angle, num_comp, plot, const = 0)


    # output

    #sersic: sersic index obtained with the method of the surface brightness at effective radius
    #meanser: mean of the sersic index obtained with the method of effective radius  
    #stdser: standard deviation of the sersic index obtained with the method of effective radius  
    #totmag: total magnitud of the galaxy
    #N: number of surface brightness model components of the galaxy
    #theta: the angle used to determine the break radius. Break radius
    #  is computed in that angle direction.




**getMissLight** computes the missing light from two surface brightness models
::

    from galfitools.galout.getMissingLight import getMissLight




    #GalfitFile1           Galfit File containing the coreless surface brightness model
    #GalfitFile2           Galfit File containing the core surface brightness model
    #rad                   upper limit of radius to integrate the missing light in pixels 


    #optional from argparse:

    #dis: Maximum distance among components
    #num_comp: Number of component where it'll obtain center of all components, default = 1


    magmiss, N1, N2 = getMissLight(galfitFile1, galfitFile2, dis, num_comp, rad)


    # output

    #N1: number of surface brightness model components of the coreless model
    #N2: number of surface brightness model components of the core model
    #magmiss: magnitude of the missing light



**getBulgeRad** gets the bulge radius or the radius where two models of surface brightness models are
equal
::

    from galfitools.galout.getRads import getBulgeRad

    #GalfitFile1           Galfit File containing the coreless surface brightness model
    #GalfitFile2           Galfit File containing the core surface brightness model

    #optional from argparse

    #dis: Maximum distance among components
    #num_comp: Number of component where it'll obtain center of all components, default = 1

    #angle:  Angle of the major axis of the galaxy. Default= it will take the angle of the
    #plot: boolean flag that indicates  to make  a plot of GalfitFile1 - GalfitFile2 vs. radius 
    #plot: boolean flag that indicates to make a plot of first derivative vs. radius
    #ranx: list that indicates the range for the plot x-axis: xmin - xmax



    rbulge, N1, N2, theta = getBulgeRad(galfitFile1, galfitFile2, dis, num_comp, angle, plot, ranx)


    # output

    #N1: number of surface brightness model components of the coreless model
    #N2: number of surface brightness model components of the core model
    #rbulge: bulge radius  




**showCube** takes the GALFIT output and creates an image that shows galaxy, model and residual 
::

    from galfitools.galout.showcube import displayCube


    #cubeimage: the cube GALFIT image 

    #optional arguments from argparse

    #namecube: name of the output image 
    #dpival: value of dpi (dots per inch)
    #brightness: brightness of the image. Only for galaxy and model. Default = 0. Preferible
    #                    range goes from -1 to 1
    #contrast: contrast of the image. Only for galaxy and model. Default = 1. Preferible
    #           range goes from 0 to 1
    #cmap: colormap to be used for the cube image (based on the matplotlib)
    #scale: plate scale of the image
    #noplot:  avoids to show the plotting window


    displayCube(cubeimage, namecube, dpival, brightness, contrast, cmap, scale, noplot)





**photDs9** computes photometry from a Ds9 region file: Box, Ellipses and Polygons
::


    from galfitools.galout.PhotDs9 import photDs9 

    args = parser.parse_args()

    #ImageFile =  the image file where the photometry will be computed

    #RegFile = the DS9 region file



    #optional for argparse

    #zeropoint: magnitude zeropoint 
    #sky: sky background value to be removed from computation 


    mag, exptime = photDs9(ImageFile, RegFile, zeropoint, sky)

    #output

    #mag: magnitud of the Ds9 regions 
    #exptime: exposition time of the image 



**fitlog2csv** converts fit.log file into a CSV file 
::


    from galfitools.galout.fitlog2csv  import log2csv 

    args = parser.parse_args()


    #optional for argparse

    #zeropoint: number of the fit to be extracted. Default: last one 
    #fileout: name of the output file 


    log2csv(num, fileout)



**getPeak**  Obtains the center, axis ratio and angular position from DS9 region
::


    from galfitools.galout.getPeak import getPeak 

    args = parser.parse_args()

    # image: image fits file
    #regfile: DS9 ellipse region file 

    #optional for argparse
    # center: optional flag to indicate that center of ds9 file will used instead 
    # maskfile: mask fits file 


    X, Y, AxRat, PA = getPeak(image, regfile, center, maskfile)


    #output

    # X, Y: position of the center (peak: coordinate with the highest value) 
    # AxRat: axis ratio of ellipse
    # PA: angular position  measured from Y-axis


**getBarSize**  Computes the bar size assuming the model is a composed of three 
Sersic function: Bulge, bar and disk. 

::


    from galfitools.galout.getBarSize import getBarSize

    #galfitFile: galfit file: galfit.XX  
    #dis: maximum distance among components to be considered as part of the same galaxy
    #out: name of the Ds9 ellipse region output file
    #num_comp: number of component to obtain center and find galaxy (if simultaneaous fitting) 
    #plot: if true, it creates plot file of radio break and radio kappa
    #ranx: search range for bar size

    rbar, N, theta =  getBarSize(galfitFile, dis, num_comp, plot, ranx, out)

    #rbar: size of bar in pixels
    #N: number of components of the galaxy
    #theta: angle of bar measured from Y-axis (same as GALFIT)


**MGE**
---------------

Routines that use the Multi-Gaussian Expansion

**mge2galfit** fits multi-gaussian expansion of Cappellari (2002) and formats to GALFIT
::

    from galfitools.mge.mge2galfit import mge2gal

    #args is an class object from argparse

    mge2gal(args) 


    #to check the args options check the -h option (shown below):

    #  positional arguments:
    #    image                 the Mask image file to modify or create
    #    Ds9regFile            the DS9 ellipse region file containing the galaxy
    #    magzpt                the magnitude zero point
    #
    #  options:
    #    -h, --help            show this help message and exit
    #    -t, --twist           uses twist option for mge
    #    -r, --regu            regularized mode for mge_fit_sectors
    #    -c, --center          uses the center given in DS9 region file,otherwise it will found the x,y
    #                          peaks within DS9 ellipse
    #    -p PSF, --psf PSF     the value of PSF sigma
    #    -s SKY, --sky SKY     the sky background value
    #    -m MASK, --mask MASK  the mask file
    #    -ps PLATE, --plate PLATE
    #                          plate scale of the image
    #    -gas, --gauss         uses gauss function for galfit file
    #    -fser, --freeser      leaves the sersic index as a free parameter to fit
    #    -fsk, --freesky       leaves the sky as a free parameter to fit
    #    -pf PSFILE, --psfile PSFILE
    #                          name of the psf file for GALFIT. default = psf.fits
    #    -sf SIGFILE, --sigfile SIGFILE
    #                          name of the sigma image for GALFIT. default = sigma.fits
    #    -ng NUMGAUSS, --numgauss NUMGAUSS
    #                          number of gaussians that will be used for galfit.Starting from the first one
    #



**SbProf** creates a surface brightness profile from a ellipse ds9 region
::


    from galfitools.mge.SbProf import sbProf

    #args is an class object from argparse

    sbProf(args)


    #to check the args options check the -h option (shown below):

    #positional arguments:
    #  Image                 image fits file
    #  Ds9Region             Ds9 ellipse region file

    #options:
    #  -h, --help            show this help message and exit
    #  -q AXRAT, --axrat AXRAT
    #                        axis ratio
    #  -pa ANGLE, --angle ANGLE
    #                        angular position (same as GALFIT)
    #  -mz MGZPT, --mgzpt MGZPT
    #                        Magnitud zero point
    #  -m MASK, --mask MASK  mask fits file
    #  -s SKY, --sky SKY     sky value. Default = 0
    #  -p PLATE, --plate PLATE
    #                        plate scale
    #  -o OUTPUT, --output OUTPUT
    #                        output file
    #  -c, --center          uses the center given in DS9 region file,otherwise it will found the x,y
    #                        peaks within DS9 ellipse
    #  -rx RANX RANX, --ranx RANX RANX
    #                        provide a range for x-axis: xmin - xmax
    #  -ry RANY RANY, --rany RANY RANY
    #                        provide a range for y-axis: ymin - ymax
    #  -lx, --logx           turn the X-axis to logarithm
    #  -px, --pix            turn the top x-axis in pixels
    #  -g, --grid            display a grid in the plot
    #  -r RAD, --rad RAD     value for a vertical line to add into the plot
    #  -r2 RAD2, --rad2 RAD2
    #                        value for a second vertical line to add into the plot

 
**SIM**
---------------
Routines that make a simulated galaxy image.

**makeSim** simulates a observed galaxy from a GALFIT model

::


    from galfitools.sim.MakeSim import makeSim

    #args = parser.parse_args()

    #image:  the GALFIT galaxy model
    #newimage:  the name of the new galaxy image

    #optional arguments from argparse

    #GAIN: the gain value of the image.

    #skymean: the sky background value.
    #skystd: the sky background value


    makeSim(image, GAIN, skymean, skystd, newimage)


**SKY**
-------------

Routines that compute the sky background

**galSky** computes the sky using GALFIT
::

    from galfitools.sky.GalfitSky import galfitSky


    # imgname: the image file
    # maskfile: the galfit mask file

    # mgzpt: magnitude zero point
    # scale: the plate scale

    # X:  the X position
    # Y: the Y position

    # initsky: the initial sky value


    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)



**getSky** computes sky from a ds9 region box file
::


    from galfitools.sky.Sky import sky
   

    # imgname: the image file
    # maskimage: The galfit mask file
    # filereg: Ds9 box region file containing the area to compute


     mean, sig = sky(imgname, maskimage, filereg)

    # mean: the mean value of the sky background
    # sig: the standard deviation the sky backround




**skyDs9** computes sky using ds9 region file
::



    from galfitools.sky.SkyDs9 import SkyDs9 

    # imgname: the image file
    # filereg: Ds9 box region file containing the area to compute



    mean, sig = SkyDs9(imgname, filereg) 





**skyRing** computes sky computing the gradient over concentric rings around the galaxy.

::


    from galfitools.sky.SkyRing import SkyRing

    # image: the image file
    # mask: The galfit mask file
    # ds9regfile: Ds9 box region file containing the area to compute


    #width: width of the rings
    #center: if True, it uses  the center indicated by the ellipse in 'ds9regfile' 



    mean, std, median, rad = SkyRing(image, mask, ds9regfile, width, center)





