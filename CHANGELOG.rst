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
- New function (getMeRad) to compute the surface brightess (and 
- mean surface brightness at a any radio
- skyRing and skyDs9 functions computes the surface brightness of the sky 
- better documentation 
- bugs fixed





