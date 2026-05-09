=========
Changelog
=========


Version 0.15.2
===============


GALFITools serves as a comprehensive 
library designed to enhance the functionality 
of GALFIT, a powerful tool for galaxy 
modeling. With GALFITools, you can achieve the following:

- Generate mask images, calculate sky backgrounds, and extract stars from your images.
- Perform essential computations such as deriving effective radius, Sersic index, slope, FWHM, and more.
- Conduct photometry analysis utilizing Ds9 regions.
- Create various visual representations including galaxy images, models, and residual images.
- Construct detailed surface brightness profiles.
- Estimate Multi-Gaussian Expansion (MGE) for optimal GALFIT utilization.
- Generate simulated galaxy data images to refine your modeling techniques.



Version 1.0.0
===============


GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the comprehensive analysis of GALFIT's output



Version 1.7.0
===============

GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output

In this release:

- fitlog2csv function was added to convert the fit.log file to a CSV file
- Computation of sky with SkyDs9 and SkyRing includes an option to remove the top 80% and bottom 20% of pixel values
- checkFile command added to check GALFIT file before run
- mge2galfit reads directly from GALFIT file to reduce arg parsing options
- creates curve of growth plot
- minor changes in names and plots
- bugs fixed


Version 1.7.7
===============


GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output

In this release:

- getPeak function was added to indicate the coordinates of the galaxy's pixel with the highest count value
- skyring function now writes a fits file to indicate the ring where the sky was computed
- bugs fixed




Version 1.15.0
================

GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output

In this release:

- The functions have been documented.
- Console scripts return results in units of pixels and arc sec
- New function GetSersic to estimate initial parameters for GALFIT input file
- New function imarith to make arithmetic operations on images
- New Function makePSF to create PSF images using multi gaussian expansion
- New function getKappa2 as an alternative to getKappa to estimate the curvature of a function
- New function getBarSize to estimate the bar size
- Bugs fixed


Version 1.17.0
================

GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output

In this release:

- New function (getBT) to compute the bulge to total luminosity ratio
- New function (getMeRad) to compute the surface brightess (and 
- mean surface brightness at a given radio

Version 1.18.0
================

GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output

In this release:

- New function (getBT) to compute the bulge to total luminosity ratio
- New function (getMeRad) to compute the surface brightess and 
  mean surface brightness at a any radio
- skyRing and skyDs9 functions computes the surface brightness of the sky 
- better documentation 
- bugs fixed


Version 1.30.0
================

GALFITools: A library for GALFIT

GALFITools constitutes a collection of specialized routines to facilitate
the utilization of GALFIT. The toolkit encompasses a spectrum of routines,
spanning from mask creation, Point Spread Function (PSF) creation,
Multi-gaussian expansion, sky computation, initial parameter estimation,
to the analysis of GALFIT's output


This release includes many changes:

- function getSersic improved
- function to convert ferrer from sersic added 
- print output format changed
- Many galout functions now have batch versions for processing multiple files.
- The makePSF function has been improved and now provides two PSF models.
- Added a batch GALFIT function that allows GALFIT to run in parallel on multiple input files.
- exp2disk function added to convert exponential function to edge-on disk function
- An alternative for mge2galfit has been added: sersic2mge. MGE from a Sersic function 
- toSersic function added to convert de Vaucoulers, exponential gaussian functions to Sersic in a GALFIT file
- New functions to download DESI data and convert innvar image to sigma image for GALFIT
- photDs9 function improved
- The maskDs9 function has been improved and now provides options to mask pixels with a specified value and invert the mask 
  outside of the DS9 region.
- new function to correct galfit files for magnitude extinction and K correction 
- new function to correct for magnitude exctincion and K correction for fitlog table created with fitlog2csv
- surface brightness can be now corrected for expansion of Universe for getrecomp and getmerad
- new functions to estimate errors on parameters from multiple GALFIT files
- getChinu function added to compute the Chinu square inside a fraction of light radius
- getChinu also computes the Akaike and Bayesian information criterions
- DESI: new function to detect central ellipse in DESI maskbits file 
- getsegvalue function added to obtain the value on a determined pixel of Segmentation Sextractor file 
- New funciton to rotate and scale a ds9 region ellipse. megamask incorporated
- megamask function for DESI added which combines three masks to create a new mask
- Added the supermask function, which creates a mask by combining the output of the maskSky function with a SExtractor segmentation image.
- Added function to convert decam 2 sdss filters
- Added new function to compute magnitude of edge disk in a GALFIT file
- pyopensci and joss badges included
- bugs fixed






