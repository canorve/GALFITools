

.. contents::
   :depth: 3
..

-------------------

**GALFITools: Usage Guide**
============================


GALFITools constitutes a collection of specialized 
routines aimed at facilitating the utilization of 
GALFIT. The toolkit encompasses a spectrum of routines, 
spanning from mask creation, Point Spread Function (PSF) 
creation, Multi-gaussian expansion, sky computation, 
initial parameter estimation, 
to the comprehensive analysis of GALFIT's output.

Upon successful installation of GALFITools, an 
ensemble of shell routines becomes seamlessly 
accessible through your bash (or zsh) environment. 
The subsequent discourse elucidates the function and 
operation of each of these routines.

For an enhanced comprehension of the intricacies 
of each routine, a explanation of their 
respective parameters can be summoned through 
the deployment of the "-h" option.


The instructional material presented herein is 
structured into five distinct sections, each 
devoted to a specific facet of GALFIT's application. 
These sections are delineated as follows: GALFIT Input, 
GALFIT Output, Multi-Gaussian Expansion (MGE), 
Simulation (Sim), and Sky Computation.

**Note**: It is paramount to underscore that 
routines necessitating the utilization of a Ds9 
region input mandate its conservation in image 
(or physical) coordinates. It is imperative to 
refrain from saving this input in the World
Coordinate System (WCS) coordinates, in order to ensure 
functionality.





**GALFIT INPUT**
------------------
Routines that aid the GALFIT's user to
prepare the necessary files for GALFIT input 


**getStar** gets a image slice centered on the object peak

::

  positional arguments:
    image                 the image file to obtain the slice
    Ds9regFile            the DS9 ellipse region file containing the star 
    size                  the size of the new image in pixels

  options:
    -h, --help            show this help message and exit
    -c, --center          uses the center given in DS9 region file,otherwise it will find the x,y
                          peaks within DS9 ellipse
    -s SKY, --sky SKY     the sky background to be removed. Default = 0
    -o OUT, --out OUT     the image output.
    -sig SIGMA, --sigma SIGMA
                          introduce the sigma image
    -so SIGOUT, --sigout SIGOUT
                          the sigma image output

**initGal** Creates GALFIT's input files with different initial parameters


::

  positional arguments:
    inFile                the galfit file galfit.XX

  options:
    -h, --help            show this help message and exit
    -n NUMBER, --number NUMBER
                          the number of files generated. Default = 1
    -p3 PARAM3 PARAM3, --param3 PARAM3 PARAM3
                          range of values to give to the 3) model's parameter in format [min max]
    -p4 PARAM4 PARAM4, --param4 PARAM4 PARAM4
                          range of values to give to the 4) model's parameter in format [min max]
    -p5 PARAM5 PARAM5, --param5 PARAM5 PARAM5
                          range of values to give to the 5) model's parameter in format [min max]
    -p6 PARAM6 PARAM6, --param6 PARAM6 PARAM6
                          range of values to give to the 6) model's parameter in format [min max]
    -p7 PARAM7 PARAM7, --param7 PARAM7 PARAM7
                          range of values to give to the 7) model's parameter in format [min max]
    -p8 PARAM8 PARAM8, --param8 PARAM8 PARAM8
                          range of values to give to the 8) model's parameter in format [min max]
    -p9 PARAM9 PARAM9, --param9 PARAM9 PARAM9
                          range of values to give to the 9) model's parameter in format [min max]
    -p10 PARAM10 PARAM10, --param10 PARAM10 PARAM10
                          range of values to give to the 10) model's parameter in format [min max]
    -nc NUMCOMP, --numcomp NUMCOMP
                          the component number which parameters will be changed
       


**gtmakeMask**  creates mask file from a SExtractor's catalog 

::

    positional arguments:
      Sexfile               Sextractor catalog file
      ImageFile             Image file

    options:
      -h, --help            show this help message and exit
      -o MASKOUT, --maskout MASKOUT
                            the output mask file name
      -sf SATDS9, --satds9 SATDS9
                            ds9 saturation file
      -s SCALE, --scale SCALE
                            scale factor to increase the ellipses. Default=1


                            *Note* The Sextractor catalog must have the following columns: 



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

  positional arguments:
    MaskFile              the Mask image file to modify or create
    RegFile               the DS9 region file

  options:
    -h, --help            show this help message and exit
    -f FILL, --fill FILL  the value in counts to fill into the Ds9 regions. Default = 0 (remove)
    -i IMAGE, --image IMAGE
                          image to obtain the size
    -b, --border          Mask the borders when their value is zero
    -bv BORVALUE, --borValue BORVALUE
                          value of the border if it is different from zero


**maskSky** creates a mask image for GALFIT using original image and sky mean and sigma

::

  positional arguments:
    ImageFile             original data image
    MaskFile              Name of the new Mask file

  options:
    -h, --help            show this help message and exit
    -sm SKYMEAN, --skymean SKYMEAN
                          mean of the sky background
    -ss SKYSIGMA, --skysigma SKYSIGMA
                          sigma of the sky background
    -ns NUMBERSIG, --numbersig NUMBERSIG
                          number of times that the sigma of the sky will be multiplied to remove the
                          sky background
    -b, --border          Mask the borders when their value is zero
    -bv BORVALUE, --borValue BORVALUE
                          value of the border if it is different from zero

**xy2fits** code to convert ASCII x,y positions to FTIS mask

::

  positional arguments:
    ImageFile          The Image file
    AsciiMask          The ascii file with the x,y positions

  options:
    -h, --help         show this help message and exit
    -c VAL, --val VAL  the value in counts for the masked pixels


**GALFIT OUTPUT**
-------------------
Routines that computes photometric variables from 
the surface brightness models fitted by GALFIT 


**getBreak** gets the break radius from a set of Sersics

::


  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -a ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components
    -ni NUMINITIAL, --numinitial NUMINITIAL
                          Number of component where it'll obtain the initial parameter to search break
                          radius or to generated random initial radius.
    -q, --quick           evaluate in position only (given by -ni parameter
    -r RANDOM, --random RANDOM
                          Number of random radius as initial parameters to search for the minimum. It
                          will generated random radius from 0 to effective radius of the component
                          indicated by parameter -ni
    -p, --plot            makes plot of double derivative vs. radius
    -rx RANX RANX, --ranx RANX RANX
                          provide a range for the plot x-axis: xmin - xmax


**getBreak2** gets the break radius from a set of Sersics using an 
alternative method to getBreak

::

  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -a ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components
    -p, --plot            makes plot of double derivative vs. radius
    -rx RANX RANX, --ranx RANX RANX
                          x-axis range to search for the Break radius: xmin - xmax



**getFWHM** gets the FWHM from a set of Sersics
::


  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -a ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components


**getKappa** gets the Kappa radius from a set of Sersics

::

  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -a ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components
    -ni NUMINITIAL, --numinitial NUMINITIAL
                          Number of component where it'll obtain the initial parameter to search break
                          radius or to generated random initial radius.
    -q, --quick           evaluate in position only (given by -ni parameter
    -r RANDOM, --random RANDOM
                          Number of random radius as initial parameters to search for the minimum. It
                          will generated random radius from 0 to effective radius of the component
                          indicated by parameter -ni
    -p, --plot            makes plot of double derivative vs. radius
    -rx RANX RANX, --ranx RANX RANX
                          provide a range for x-axis: xmin - xmax




**getReComp** gets the effective radius from a set of Sersics
::

  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -er EFFRAD, --effrad EFFRAD
                          percentage of light to compute for radius. default=.5 for effective radius
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -pa ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components. Angle measured from Y-Axis as same as GALFIT.



**getSlope** gets the slope radius from a set of Sersics
::


  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -a ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components
    -s SLOPE, --slope SLOPE
                          value of slope to find. default=.5
    -p, --plot            makes plot of double derivative vs. radius
    -rx RANX RANX, --ranx RANX RANX
                          provide a range for x-axis: xmin - xmax




**getN** computes the Sersic index from surface brightness at effective radius
::

  positional arguments:
    GalfitFile            Galfit File containing the Sersics or gaussians components

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -pa ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components. Angle measured from Y-Axis as same as GALFIT.
    -rf RADFRAC, --radfrac RADFRAC
                          fraction of light radius. Default = .2
    -p, --plot            makes plot of double derivative vs. radius



**getMissLight** computes the missing light from two surface brightness models
::

  positional arguments:
    GalfitFile1           Galfit File containing the coreless surface brightness model
    GalfitFile2           Galfit File containing the core surface brightness model
    rad                   upper limit of radius to integrate the missing light in pixels 

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1



**getBulgeRad** gets the bulge radius or the radius where two models of surface brightness models are
equal
::

  positional arguments:
    GalfitFile1           Galfit File containing the coreless surface brightness model
    GalfitFile2           Galfit File containing the core surface brightness model

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -pa ANGLE, --angle ANGLE
                          Angle of the major axis of the galaxy. Default= it will take the angle of the
                          last components. Angle measured from Y-Axis as same as GALFIT.
    -p, --plot            makes plot of double derivative vs. radius
    -rx RANX RANX, --ranx RANX RANX
                          provide a range for x-axis: xmin - xmax


**showCube** takes the GALFIT output and creates an image that shows galaxy, model and residual 
::

  positional arguments:
    cubeimage             the cube GALFIT image

  options:
    -h, --help            show this help message and exit
    -o OUTIMAGE, --outimage OUTIMAGE
                          the output png file
    -br BRIGHTNESS, --brightness BRIGHTNESS
                          brightness of the image. Only for galaxy and model. Default = 0. Preferible
                          range goes from -1 to 1
    -co CONTRAST, --contrast CONTRAST
                          contrast of the image. Only for galaxy and model. Default = 1. Preferible
                          range goes from 0 to 1
    -cm CMAP, --cmap CMAP
                          cmap to be used for the cube image
    -dpi DOTSINCH, --dotsinch DOTSINCH
                          dots per inch used for images files
    -s SCALE, --scale SCALE
                          plate scale of the image. Default = 1
    -np, --noplot         it doesn't show plotting window


**photDs9** computes photometry from a Ds9 region file: Box, Ellipses and Polygons
::


  positional arguments:
    ImageFile             the image file where the photometry will be computed
    RegFile               the DS9 region file

  options:
    -h, --help            show this help message and exit
    -zp ZEROPOINT, --zeropoint ZEROPOINT
                          The value of the zero point. Default = 25
    -sk SKY, --sky SKY    the value of the sky background to be removed


**MGE**
---------------

Routines that use the Multi-Gaussian Expansion.

**mge2galfit** fits multi-gaussian expansion of Cappellari (2002) and formats to GALFIT
::

  positional arguments:
    image                 the Mask image file to modify or create
    Ds9regFile            the DS9 ellipse region file containing the galaxy
    magzpt                the magnitude zero point

  options:
    -h, --help            show this help message and exit
    -t, --twist           uses twist option for mge
    -r, --regu            regularized mode for mge_fit_sectors
    -c, --center          uses the center given in DS9 region file,otherwise it will found the x,y
                          peaks within DS9 ellipse
    -p PSF, --psf PSF     the value of PSF sigma
    -s SKY, --sky SKY     the sky background value
    -m MASK, --mask MASK  the mask file
    -ps PLATE, --plate PLATE
                          plate scale of the image
    -gas, --gauss         uses gauss function for galfit file
    -fser, --freeser      leaves the sersic index as a free parameter to fit
    -fsk, --freesky       leaves the sky as a free parameter to fit
    -pf PSFILE, --psfile PSFILE
                          name of the psf file for GALFIT. default = psf.fits
    -sf SIGFILE, --sigfile SIGFILE
                          name of the sigma image for GALFIT. default = sigma.fits
    -ng NUMGAUSS, --numgauss NUMGAUSS
                          number of gaussians that will be used for galfit.Starting from the first one

**SbProf** creates a surface brightness profile from a ellipse ds9 region
::

  positional arguments:
    Image                 image fits file
    Ds9Region             Ds9 ellipse region file

  options:
    -h, --help            show this help message and exit
    -q AXRAT, --axrat AXRAT
                          axis ratio
    -pa ANGLE, --angle ANGLE
                          angular position (same as GALFIT)
    -mz MGZPT, --mgzpt MGZPT
                          Magnitud zero point
    -m MASK, --mask MASK  mask fits file
    -s SKY, --sky SKY     sky value. Default = 0
    -p PLATE, --plate PLATE
                          plate scale
    -o OUTPUT, --output OUTPUT
                          output file
    -c, --center          uses the center given in DS9 region file,otherwise it will found the x,y
                          peaks within DS9 ellipse
    -rx RANX RANX, --ranx RANX RANX
                          provide a range for x-axis: xmin - xmax
    -ry RANY RANY, --rany RANY RANY
                          provide a range for y-axis: ymin - ymax
    -lx, --logx           turn the X-axis to logarithm
    -px, --pix            turn the top x-axis in pixels
    -g, --grid            display a grid in the plot
    -r RAD, --rad RAD     value for a vertical line to add into the plot
    -r2 RAD2, --rad2 RAD2
                          value for a second vertical line to add into the plot

**SIM**
---------------
Routines that make a simulated galaxy image using GALFIT.

**makeSim** simulates a observed galaxy from a GALFIT model. It 
adds Poisson and sky noise to the image.
::

  positional arguments:
    image                 the GALFIT galaxy model
    newimage              the name of the new galaxy image

  options:
    -h, --help            show this help message and exit
    -s SKY, --sky SKY     the sky background value. default = 0
    -std STD, --std STD   the sky standard deviation. default = 1
    -g GAIN, --gain GAIN  the gain value of the image. default = 1

**SKY**
-------------

Routines that compute the sky background.

**galSky** computes the sky using GALFIT
::

  positional arguments:
    image                 the image file
    mask                  the GALFIT mask file

  options:
    -h, --help            show this help message and exit
    -s SCALE, --scale SCALE
                          the plate scale. default = 1
    -zp MGZPT, --mgzpt MGZPT
                          the magnitud zero point. default=25
    -x XPOS, --xpos XPOS  the x position. default=1
    -y YPOS, --ypos YPOS  the y position. default=1
    -is INITSKY, --initsky INITSKY
                          the initial sky value default=0

**getSky** computes sky from a ds9 region box file
::

  positional arguments:
    image       the image file 
    maskfile    the Mask image file 
    Ds9regFile  the DS9 box region file containing the galaxy

  options:
    -h, --help  show this help message and exit


**skyDs9** computes sky using ds9 region file
::

  positional arguments:
    image       the image file 
    Ds9regFile  the DS9 box region file containing the galaxy

  options:
    -h, --help  show this help message and exit



**skyRing** computes sky computing the gradient over concentric rings
around the galaxy.

::

  positional arguments:
    image       the image file 
    maskfile    the Mask image file 
    Ds9regFile  the DS9 box region file containing the galaxy


  options:
    -h, --help  show this help message and exit
    -c, --center  use the center of the ellipse. Otherwise it will use the (x,y) position with the highest value of the ellipse


--------------

**Questions?**
--------------


Something is not clear for you or do you have further questions?
write to me at canorve [at] gmail [dot] com 

