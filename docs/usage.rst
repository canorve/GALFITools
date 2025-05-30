
.. _usage:

**GALFITools: Usage Guide**
============================


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



---------------------

**Introduction**
---------------------

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
The text below elucidates the function and 
operation of each of these routines.

For a better understanding of the intricacies 
of each routine, an explanation of their respective 
parameters can be accessed by using the '-h' option

This usage guide explains how to use the GALFITools routines. 
All of these routines are *shell commands*. If you wish to use 
the Python functions underlying these routines within 
your own scripts, please check the :doc:`api` or :doc:`api/modules` module
reference.

Many of the routines (particularly the **GALFIT Output** routines) 
require the GALFIT output file, typically named "galfit.XX", where 
"XX" refers to a numerical identifier. These files are generated 
after GALFIT completes a fitting process. They serve as input for 
parsing data within these routines.

Components within the GALFIT file refer to individual surface brightness 
models, each following a specific order within 
the GALFIT file. A galaxy can be composed of a 
single component or multiple components. For instance,
a barred galaxy can be model by 3 components: bulge, bar
and disk. 

To differentiate the components of different galaxies, 
components belonging to the same galaxy must share the 
same center. This is determined by a tolerance parameter 
specified by the '-d' flag, which indicates the maximum 
allowable distance among components (default 10 pixels).


.. important:: Soome routines require the use of a `DS9 <http://tdc-www.harvard.edu/saoimage/>`_ 
               region file as a input. The region file needs to be saved in 
               **image** (or **physical**) coordinates. DS9 regions 
               saved in **World Coordinate System (WCS)** coordinates 
               will not work.



The instructional material presented herein is 
structured into five distinct sections, each 
devoted to a specific facet of GALFIT's application. 
These sections are delineated as follows: GALFIT Input, 
GALFIT Output, Multi-Gaussian Expansion (MGE), 
Simulation of galaxy images, and Sky Computation.


Before presenting the five sections of the instructional 
material, check the  three video tutorials below which show 
the full capabilities of GALFITools in action.

---------------------

**Video Tutorials**
---------------------

For the videos, click on the images and turn on English subtitles.

Here is a tutorial video for galaxy NGC720 using GALFITools:


.. image:: https://img.youtube.com/vi/2npeGmC1mCg/maxresdefault.jpg
    :alt: IMAGE ALT TEXT HERE
    :target: https://www.youtube.com/watch?v=2npeGmC1mCg



Here is another for galaxy NGC1198 using GALFITools:

.. image:: https://img.youtube.com/vi/VmJJkKVd37U/maxresdefault.jpg
    :alt: IMAGE ALT TEXT HERE
    :target: https://www.youtube.com/watch?v=VmJJkKVd37U


Here we model the barred spiral galaxy PGC 34232
with GALFITools v1.11.0:

.. image:: https://img.youtube.com/vi/wUA-sigVSts/maxresdefault.jpg
    :alt: IMAGE ALT TEXT HERE
    :target: https://www.youtube.com/watch?v=wUA-sigVSts



---------------------

**GALFIT Input**
--------------------

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
                          peak within DS9 ellipse
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




**maskDs9**  creates (or modify) a mask image for GALFIT using DS9 regions 
such as Boxes, Ellipses and Polygons

::

  usage: maskDs9 [-h] [-f FILL] [-i IMAGE] [-b] [-bv BORVALUE] MaskFile RegFile


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

  usage: maskSky [-h] [-sm SKYMEAN] [-ss SKYSIGMA] [-ns NUMBERSIG] [-b] [-bv BORVALUE]
               ImageFile MaskFile


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

  usage: xy2fits [-h] [-c VAL] ImageFile AsciiMask

  positional arguments:
    ImageFile          The Image file
    AsciiMask          The ascii file with the x,y positions

  options:
    -h, --help         show this help message and exit
    -c VAL, --val VAL  the value in counts for the masked pixels




**checkFile** check that the parameters and file names inside the GALFIT input file are correct 

::

  usage: checkFile [-h] [-d DIS] GalfitFile

  positional arguments:
    GalfitFile         GALFIT input File

  options:
    -h, --help         show this help message and exit
    -d DIS, --dis DIS  Maximum distance in pixels among components. Default = 10


**boxSize** computes the box size from a ds9 box region for galfit header option H) 

::

    
  usage: boxSize [-h] RegFile

  Computes the Box size from a Ds9 region file for galfit header

  positional arguments:
    RegFile     Ds9 region file containing the box region

  options:
    -h, --help  show this help message and exit



**getPeak**  Obtains the center, axis ratio and angular position from DS9 region

::

    positional arguments:
      Image                 image fits file
      RegFile               DS9 ellipse region file

    options:
      -h, --help            show this help message and exit
      -c, --center          takes center of ds9 region file
      -m MASK, --mask MASK  the mask file


**imarith** makes arithmetic operations on image 

::

    
  usage: imarith [-h] ImageFile 

  makes arithmetic operations on image 

  positional arguments:
    ImageFile   The input image

  options:

   -h, --help  show this help message and exit


   -o --output    The output image

   -i  --image2   second input image to make arithmetic operations with ImageFile. Image2 must be of the same size of ImageFile. If this second image is provided it will make operations indicated by arithmetic flag ignoring its constant input

   -a   --add   add constant to image pixels
   -d   --div   divide all pixels by constant
   -m   --mul   multiply all pixels by constant
   -s   --sub   substract constant to all pixels


 

**getSersic** Its estimates and prints initial parameters for Sersic components. It
              addtion if proved options for single Sersic, bulge/disk and bulge/bar/disk

::



  usage: getSersic [-h] [-zp ZEROPOINT] [-sk SKY] [-bt BULGETOT] [-c] [-n] [-m MASK] [-b BARDS9]
                 Image RegFile

  prints the Sersic function from DS9 ellipse region

  positional arguments:
    Image                 image fits file
    RegFile               DS9 ellipse region file

  options:
    -h, --help            show this help message and exit
    -zp ZEROPOINT, --zeropoint ZEROPOINT
                          The value of the zero point. Default = 25
    -sk SKY, --sky SKY    Sky background value to be removed from image before photometry. Default = 0
    -bt BULGETOT, --bulgetot BULGETOT
                          Bulge to total ratio. If set it will print two sersics: one for the bulge and
                          the other for the disk
    -c, --center          takes center of ds9 region file
    -n, --noprint         avoids to print Sersic functionts to stdout
    -m MASK, --mask MASK  the mask file
    -b BARDS9, --bards9 BARDS9
                          DS9 ellipse region file that containts the bar region. bulgetot flag must be
                          activated


**MakePSF** Makes a PSF model of a star using Multi Gaussian Expansion

::

  usage: makePSF [-h] [-c] [-o OUT] [-sig SIGMA] [-t] [-ng NUMGAUSS] image GalfitFile Ds9regFile


  positional arguments:
    image                 the image file where it contains the star to be modelled
    GalfitFile            GALFIT file to obtain the header options
    Ds9regFile            the DS9 ellipse region file containing the star to model

  options:
    -h, --help            show this help message and exit
    -c, --center          uses the center given in DS9 region file, otherwise it will find the (x,y) peak
                          within DS9 ellipse
    -o OUT, --out OUT     the PSF model image
    -sig SIGMA, --sigma SIGMA
                          introduce the sigma image
    -t, --twist           uses twist option for mge
    -ng NUMGAUSS, --numgauss NUMGAUSS
                          number of gaussians that will be used for galfit.


---------------------

**GALFIT Output**
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



**getKappa2** gets the kappa radius from a set of Sersics using an 
alternative method to getKappa

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


  usage: photDs9 [-h] [-zp ZEROPOINT] [-m MASK] [-sk SKY] ImageFile RegFile

  positional arguments:
    ImageFile             the image file where the photometry will be computed
    RegFile               the DS9 region file

  options:
    -h, --help            show this help message and exit
    -zp ZEROPOINT, --zeropoint ZEROPOINT
                          The value of the zero point. Default = 25
    -sk SKY, --sky SKY    the value of the sky background to be removed



**fitlog2csv** converts fit.log file into a comma separated values file 
::

  usage: fitlog2csv [-h] [-o FILEOUT] [-n NUM]

    -h, --help            show this help message and exit
    -n NUM, --NUM NUM     the number of the fit to be extracted 
    -o OUTFILE, --fileout OUTFILE 
                          the name of the output file 



**getBarSize** computes the barsize from a composed Sersic model: Bulge, bar and disk.
. It assumes that bar is a gaussian (Sersic index = 0.5) and it is positioned as a 
second component in galfit file galfit.XX

::


  usage: getBarSize [-h] [-d DIS] [-n NUMCOMP] [-o OUT] [-p] [-rx RANX RANX] GalfitFile

  getBarSize: gets the bar size from Sersic models: Bulge, Bar and disk. It assumes that bar is galfit
  component 2

  positional arguments:
    GalfitFile            Galfit File containing the Sersics components bulge, bar, disk

  options:
    -h, --help            show this help message and exit
    -d DIS, --dis DIS     Maximum distance among components
    -n NUMCOMP, --numcomp NUMCOMP
                          Number of component where it'll obtain center of all components, default = 1
    -o OUT, --out OUT     output DS9 ellipse region file
    -p, --plot            makes plot of double derivatives and Kappa radius
    -rx RANX RANX, --ranx RANX RANX
                          x-axis range to search for the Break and Kappa radius: xmin - xmax


-----------------

**MGE**
---------------

Routines that use the Multi-Gaussian Expansion.

**mge2galfit** fits multi-gaussian expansion of Cappellari (2002) and formats to GALFIT
::

  positional arguments:
    GalfitFile            GALFIT file to obtain the header options
    Ds9regFile            the DS9 ellipse region file containing the galaxy

  options:
    -h, --help            show this help message and exit
    -t, --twist           uses twist option for mge
    -c, --center          uses the center given in DS9 region file,otherwise it will found the x,y peak within DS9
                          ellipse
    -p PSF, --psf PSF     the value of PSF sigma
    -gas, --gauss         uses gauss function for galfit file
    -fser, --freeser      leaves the sersic index as a free parameter to fit
    -fsk, --freesky       leaves the sky as a free parameter to fit
    -ng NUMGAUSS, --numgauss NUMGAUSS
                         




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
                          peak within DS9 ellipse
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



-----------------

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



-------------

**Sky**
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



