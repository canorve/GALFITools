
Concepts
========

This page provides a short overview of common astronomy and image–analysis
terms used throughout the GALFITools documentation. The goal is to make
the package more accessible for users who may be less familiar with
astronomical jargon.



-------------------------
**Computational Tools**
-------------------------

.. _concept-galfit:

GALFIT
------
GALFIT is a widely used program for two–dimensional fitting of galaxy
images (Peng et al. 2002). It models galaxies with analytic functions
(e.g., Sérsic, exponential disks, Nuker profiles) convolved with the
point–spread function (PSF). GALFIT is the external engine that
GALFITools interfaces with.

**Related GALFITools API**

- :py:func:`galfitools.galin.MakePSF.makePSF`
- :py:func:`galfitools.mge.mge2galfit.mge2gal`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`makePSF <routine-makePSF>`
   - :ref:`mge2galfit <routine-mge2galfit>`



.. _concept-galfit-file:

GALFIT file
-----------
A GALFIT file (often called a *feedme file*) is a text file that
contains all of the information required for GALFIT to run a model fit.
It specifies the image properties (such as image filename, plate scale,
and photometric zero point), the fitting and convolution regions, and
the initial parameters of each component (e.g., magnitude, effective
radius, Sérsic index). After a run, the file can also record the final
fitted parameters. GALFITools provides utilities to create and manage
these files programmatically.

**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`


Example of a GALFIT file
------------------------

Below is a simplified snippet of a GALFIT input file (*feedme file*).  
Each line specifies a property of the image or a parameter of a model
component. The symbols in brackets (A, B, etc.) indicate whether a
parameter is free or fixed during fitting.

.. code-block:: text

   =============================================================================
   # IMAGE and GALFIT CONTROL PARAMETERS
   A) galaxy.fits            # Input data image (FITS file)
   B) output.fits            # Output data image block
   C) none                   # Sigma image name (or "none")
   D) psf.fits               # Input PSF image
   E) 1                      # PSF oversampling factor
   F) mask.fits              # Bad pixel mask (or "none")
   G) none                   # Parameter constraint file (or "none")
   H) 1    512   1    512    # Image fitting region (xmin xmax ymin ymax)
   I) 100    100             # Size of convolution box (x y)
   J) 25.0000             # Magnitude photometric zeropoint 
   K) 0.2500  0.2500      # Plate scale (dx dy)   [arcsec per pixel]
   O) regular             # Display type (regular, curses, both)
   P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps


   # Component 1: Sérsic profile
   0) sersic                 # Object type
   1) 256.0   256.0   1 1    # Position x, y [pixel]
   3) 18.5     1             # Total magnitude
   4) 20.0     1             # Effective radius [pixels]
   5) 2.5      1             # Sérsic index
   9) 0.9      1             # Axis ratio (b/a)
  10) 45.0     1             # Position angle (degrees)
   Z) 0                      # Skip this model in output image? (no=0)



**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`

.. _concept-ds9:

SaoImage DS9
------------
SaoImage DS9 is an astronomical imaging and visualization application.
It is commonly used to inspect FITS images, define regions of interest,
and interactively examine astronomical data.

**Related GALFITools API**

- :py:func:`galfitools.galin.MaskDs9.maskDs9`
- :py:func:`galfitools.galout.PhotDs9.photDs9`
- :py:func:`galfitools.sky.SkyDs9.SkyDs9`

.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`maskDs9 <routine-maskDs9>`
   - :ref:`photDs9 <routine-photDs9>`
   - :ref:`skyDs9 <routine-skyDs9>`



.. _concept-ds9-regions:

DS9 regions
-----------
In DS9, a *region* is a user–defined geometric shape (circle, box,
polygon, etc.) drawn on an image. Regions can mark sources, masks,
or fitting boundaries, and can be saved to files that GALFITools 
can read. In GALFITools, DS9 regions must be saved in physical or image
coordinates.


**Related GALFITools API**

- :py:func:`galfitools.galin.MaskDs9.maskDs9`
- :py:func:`galfitools.galout.PhotDs9.photDs9`
- :py:func:`galfitools.sky.SkyDs9.SkyDs9`

.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`maskDs9 <routine-maskDs9>`
   - :ref:`photDs9 <routine-photDs9>`
   - :ref:`skyDs9 <routine-skyDs9>`




.. _concept-mask:

Mask
----
A mask is an image that flags which pixels should be ignored during
fitting (for example, bright foreground stars, cosmic rays, or image
defects). GALFITools includes utilities to create masks automatically or
manually.

**Related GALFITools API**

- :py:func:`galfitools.galin.MaskDs9.maskDs9`
- :py:func:`galfitools.galin.MaskSky.skyRem`
- :py:func:`galfitools.galin.MakeMask.makeMask`



.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`maskDs9 <routine-maskDs9>`
   - :ref:`maskSky <routine-maskSky>`
   - :ref:`gtmakeMask <routine-gtmakeMask>`


.. _concept-initial-params:

Initial parameters
------------------
Initial parameters are starting guesses for the fitting procedure, such
as galaxy magnitude, effective radius, axis ratio, and Sérsic index.
Providing reasonable initial values helps GALFIT converge to the correct
solution.


**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`



.. _concept-sb-model:

Surface brightness model
------------------------
A surface brightness model is a mathematical description of how the
light distribution of a galaxy is represented. Models are constructed
by combining analytic functions such as Sérsic profiles, de Vaucouleurs
laws, exponential disks, Gaussians, or Nuker profiles. Each function
contributes to the total brightness distribution, and the complete
model can include one or many components.

**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`



.. _concept-component:

Component
---------
A component refers to a single analytic function used within a surface
brightness model, for example one Sérsic profile, one exponential disk,
or one Gaussian. A model may consist of a single component (e.g. one
Sérsic function) or several components combined (e.g. bulge + disk,
or bulge + disk + bar). Components allow complex galaxy structures to be
described in a modular way.

**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`



.. _concept-sbp:

Surface brightness profile
--------------------------
A surface brightness profile is a one–dimensional curve showing how the
surface brightness of a galaxy changes with radius. Profiles are often
used to characterize galaxy structure and to fit analytic models.


**Related GALFITools API**

- :py:func:`galfitools.mge.SbProf.sbProf`
- :py:func:`galfitools.galout.showcube.displayCube`
- :py:func:`galfitools.galout.getCOW.getCOW`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`SbProf <routine-SbProf>`
   - :ref:`ShowCube <routine-showCube>`
   - :ref:`getCOW <routine-getCOW>`


.. _concept-simultaneous-fitting:

Simultaneous fitting
--------------------
In GALFIT, *simultaneous fitting* refers to fitting **two or more galaxies at
the same time within the same image**. Each galaxy may have its own surface
brightness model (bulge, disk, bar, etc.), and all models are optimized
together in a single run.

**Why use it.**  
- To correctly account for light contamination when galaxies overlap or lie
  close to one another in projection.  
- To avoid bias in parameter estimates that would occur if neighboring galaxies
  were ignored or subtracted beforehand.  
- To ensure that the global sky, PSF, and noise properties are treated
  consistently across the field.

**How GALFITools interprets simultaneous fitting.**  
- Components belonging to the *same* galaxy are grouped by a common center
  (within a distance tolerance given by ``-d/--dis``).  
- Components belonging to *different* galaxies have distinct centers.  
- Many GALFITools routines allow you to specify which galaxy (set of components)
  to analyze via ``--numcomp`` or similar options.  
- The software can therefore compute quantities for one galaxy (e.g., bulge to
  total ratio) even if several galaxies were fit simultaneously.

.. note::
   This is different from fitting multiple components of a single galaxy.
   Here, **entire galaxies** (each possibly with multiple components) are fit
   together.

**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getReComp`
- :py:func:`galfitools.galout.getMeRad.getMeRad`

.. admonition:: Related GALFITools CLI commands
   :class: dropdown

   - :ref:`routine-getReComp`
   - :ref:`routine-getMeRad`


------------------------------
**Photometric Parameters**
------------------------------

.. _concept-star-image:

Star (in an image)
------------------
In the context of astronomical images, a *star* usually appears as a
point–like source broadened by the PSF. Stars are often used to
construct PSFs or to calibrate the image.

**Related GALFITools API**

- :py:func:`galfitools.galin.getStar.getStar`
- :py:func:`galfitools.galin.MakePSF.makePSF`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getStar <routine-getStar>`
   - :ref:`makePSF <routine-makePSF>`


.. _concept-psf:

Point–spread function (PSF)
---------------------------
The point–spread function represents the response of the telescope and
detector to a point source of light (such as a star). It describes the
blurring of the image due to the instrument and atmosphere (for
ground–based observations). GALFIT requires a PSF to deconvolve models.


**Related GALFITools API**

- :py:func:`galfitools.galin.getStar.getStar`
- :py:func:`galfitools.galin.MakePSF.makePSF`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getStar <routine-getStar>`
   - :ref:`makePSF <routine-makePSF>`




.. _concept-fwhm:

Full Width at Half Maximum (FWHM)
---------------------------------
FWHM is a measure of the width of a profile. It is the distance between
the two points on the profile where the value falls to half of its
maximum. In astronomy, the FWHM of a stellar image provides an estimate
of the PSF size and image resolution.

**Related GALFITools API**

- :py:func:`galfitools.galin.getStar.getStar`
- :py:func:`galfitools.galin.MakePSF.makePSF`
- :py:func:`galfitools.galout.getRads.getFWHM`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getStar <routine-getStar>`
   - :ref:`makePSF <routine-makePSF>`
   - :ref:`getFWHM <routine-getFWHM>`


.. _concept-surface-brightness:

Surface brightness (SB)
-----------------------
Surface brightness is the flux received per unit area on the sky,
commonly reported in magnitudes per square arcsecond (mag arcsec⁻²).
It describes how light is distributed across a galaxy image.

**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getReComp`
- :py:func:`galfitools.galout.getMeRad.getMeRad`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getReComp <routine-getReComp>`
   - :ref:`getMeRad <routine-getMeRad>`



.. _concept-photometry:

Photometry
----------
Photometry is the measurement of fluxes or magnitudes of astronomical
objects. It can be performed with apertures, PSF fitting, or model
fitting methods such as GALFIT.

**Related GALFITools API**

- :py:func:`galfitools.galout.PhotDs9.photDs9`

.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`photDs9 <routine-photDs9>`





.. _concept-simulated-galaxy:

Simulated galaxy (photometric)
------------------------------
A simulated galaxy is a synthetic image constructed using analytic
profiles (e.g., Sérsic, exponential disk) and observational effects such
as PSF convolution and noise. Simulated galaxies are used for testing,
teaching, and validating analysis pipelines.

**Related GALFITools API**

- :py:func:`galfitools.sim.MakeSim.makeSim`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`makeSim <routine-makeSim>`




.. _concept-sky:



Sky background
----------------
The sky background is the level of light in an image not associated with
the target object. Accurate background estimation is important because
it strongly affects measured magnitudes and profiles.

**Related GALFITools API**

- :py:func:`galfitools.sky.SkyDs9.SkyDs9`
- :py:func:`galfitools.sky.GalfitSky.galfitsky`
- :py:func:`galfitools.sky.SkyRing.SkyRing`



.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`skyDs9 <routine-skyDs9>`
   - :ref:`galSky <routine-galSky>`
   - :ref:`skyRing <routine-skyRing>`


.. _concept-magnitude:

Magnitude
-----------
Magnitude is a logarithmic measure of the brightness of an astronomical
object. A decrease of 1 magnitude corresponds to an increase in
brightness by a factor of about 2.512. Fainter objects have larger
magnitude values, while brighter objects have smaller values.


**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getReComp`
- :py:func:`galfitools.galout.PhotDs9.photDs9`
- :py:func:`galfitools.galout.getMeRad.getMeRad`
- :py:func:`galfitools.galout.getN.getN`
- :py:func:`galfitools.galout.get.getCOW`

.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getReComp <routine-getReComp>`
   - :ref:`photDs9 <routine-photDs9>`
   - :ref:`getMeRad <routine-getMeRad>`
   - :ref:`getN <routine-getN>`
   - :ref:`getBT <routine-getBT>`
   - :ref:`getCOW <routine-getCOW>`



.. _concept-convolution:

Convolution
------------
In image analysis, convolution is the process of combining two
functions, such as a model galaxy image and the point–spread function
(PSF), to simulate how the model would appear through a telescope and
detector. GALFIT uses convolution to compare model components with the
observed data.

**Related GALFITools API**

- :py:func:`galfitools.galin.getStar.getStar`
- :py:func:`galfitools.galin.MakePSF.makePSF`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getStar <routine-getStar>`
   - :ref:`makePSF <routine-makePSF>`


.. _concept-sigma-image:

Sigma image
-----------
A sigma image is an auxiliary image where each pixel value represents
the estimated standard deviation (uncertainty) of the corresponding
pixel in the science image. GALFIT can use a sigma image to weight the
fit, giving less importance to noisy pixels.



.. _concept-mag-zero:

Magnitude zero point
---------------------
The magnitude zero point is a calibration constant that converts between
instrumental fluxes (in counts or electrons) and standard magnitudes. It
depends on the instrument, filter, and exposure time. A correct zero
point ensures that fitted magnitudes can be compared with standard
photometric systems.



.. _concept-plate-scale:

Plate scale
-----------
The plate scale is the conversion factor between pixel units in the
image and angular units on the sky, usually expressed in arcseconds per
pixel. It depends on the telescope optics and detector.

**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getReComp`
- :py:func:`galfitools.galout.getMeRad.getMeRad`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getReComp <routine-getReComp>`
   - :ref:`getMeRad <routine-getMeRad>`





-------------------------
**Galactic Components**
-------------------------


.. _concept-bulge:

Bulge
-------

The bulge is the central, spheroidal component of a galaxy. It is
generally more concentrated and has higher surface brightness than the
surrounding disk.

**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`



.. _concept-disk:

Disk
-------

The disk is the flattened, rotating component of a galaxy, typically
hosting spiral arms and ongoing star formation. Its brightness profile
is often well described by an exponential law.


**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`



.. _concept-bar:

Galactic bar
----------------

A bar is an elongated structure of stars crossing the central region of
a disk galaxy. Bars redistribute angular momentum and can drive gas
inflows toward the galaxy center.


**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`
- :py:func:`galfitools.galin.checkGalFile.checkFile`
- :py:func:`galfitools.galout.getBarSize.getBarSize`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`
   - :ref:`checkFile <routine-checkFile>`
   - :ref:`getBarSize <routine-getBarSize>`





--------------------------------
**Surface Brightness Modeling**
--------------------------------

.. _concept-sersic:

Sérsic profile and index
------------------------

A Sérsic profile is a mathematical function that describes how the
brightness of a galaxy varies with radius. The **Sérsic index** (*n*)
controls the shape of the profile: low *n* values represent disk–like,
exponential profiles, while high *n* values represent more concentrated,
bulge–like profiles.

**Related GALFITools API**

- :py:func:`galfitools.galin.getSersic.getSersic`
- :py:func:`galfitools.galin.initgal.InitGal`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSersic <routine-getSersic>`
   - :ref:`initgal <routine-initGal>`




.. _concept-mge:

Multi–Gaussian Expansion (MGE)
------------------------------
The Multi–Gaussian Expansion (MGE) method represents a complex
two–dimensional light distribution as a sum of multiple two–dimensional
Gaussian functions. This approach provides a flexible but compact way to
model galaxy surface brightness profiles and is often used as input for
dynamical modeling.

**Related GALFITools API**

- :py:func:`galfitools.galin.MakePSF.makePSF`
- :py:func:`galfitools.mge.mge2galfit.mge2gal`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`makePSF <routine-makePSF>`
   - :ref:`mge2galfit <routine-mge2galfit>`




.. _concept-nuker:

Nuker profile
---------------
The Nuker profile is a broken power–law function used to describe the
inner surface brightness distribution of galaxies, especially elliptical
galaxies and bulges. It is defined by an inner slope, an outer slope,
a break radius where the transition occurs, and a smoothness parameter
that controls how sharp the transition is. The Nuker profile was
introduced by Lauer et al. (1995) to fit the central light profiles of
early–type galaxies observed with the Hubble Space Telescope.

**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getBreak`
- :py:func:`galfitools.galout.getRads.getBreak2`
- :py:func:`galfitools.galout.getRads.getSlope`

.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getBreak <routine-getBreak>`
   - :ref:`getBreak2 <routine-getBreak2>`
   - :ref:`getSlope <routine-getSlope>`


---------------------------------
**Structural Parameters**
---------------------------------

.. _concept-bt:

Bulge–to–total luminosity ratio (B/T)
---------------------------------------

The bulge–to–total luminosity ratio is the fraction of a galaxy’s total
light that comes from the bulge compared to the sum of bulge and disk.
It is commonly used to quantify galaxy morphology.

**Related GALFITools API**

- :py:func:`galfitools.galout.getBT.getBT`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getBT <routine-getBT>`


.. _concept-effective-radius:

Effective radius (Re)
-----------------------

The effective radius is the radius of a circular aperture that contains
half of the total light of a galaxy or model component. It is a standard
measure of galaxy size.


**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getReComp`
- :py:func:`galfitools.galout.getMeRad.getMeRad`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getReComp <routine-getReComp>`
   - :ref:`getMeRad <routine-getMeRad>`



.. _concept-axis-ratio:

Axis ratio
----------
The axis ratio (*b/a*) is the ratio of the minor axis length (*b*) to
the major axis length (*a*) of an ellipse that describes the projected
shape of a galaxy component. Values near 1 correspond to nearly circular
objects, while smaller values indicate more elongated shapes.





.. _concept-break-radius:

Break radius (Nuker function)
-----------------------------
In the Nuker profile, the break radius is the scale at which the slope
of the surface brightness profile changes from the inner power–law to
the outer power–law regime.


**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getBreak`
- :py:func:`galfitools.galout.getRads.getBreak2`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getBreak <routine-getBreak>`
   - :ref:`getBreak2 <routine-getBreak2>`



.. _concept-gamma-radius:

Gamma radius (Nuker function)
-----------------------------
The gamma radius is defined as the radius where the negative logarithmic
slope of the Nuker profile equals 0.5. It is used as a scale indicator
for the transition between the core and outer regions.

**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getSlope`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSlope <routine-getSlope>`


.. _concept-slope:

Slope of a function
-------------------
The slope of a function is the rate at which the function changes with
respect to its variable. In logarithmic surface brightness profiles,
slopes describe how steeply brightness declines with radius.


**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getSlope`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getSlope <routine-getSlope>`




.. _concept-kappa:

Kappa (κ)
---------
In mathematics, κ (kappa) is commonly used to denote curvature. For a
curve, curvature measures how quickly the direction of the tangent
changes with position. In galaxy dynamics, κ often appears as the
*epicyclic frequency*, describing the radial oscillations of stars in a
disk around their guiding center orbit. See
`Wikipedia: Curvature <https://en.wikipedia.org/wiki/Curvature>`_ for
the mathematical definition.




**Related GALFITools API**

- :py:func:`galfitools.galout.getRads.getKappa`
- :py:func:`galfitools.galout.getRads.getKappa2`


.. admonition:: Related GALFITools CLI commands
   :class: seealso

   - :ref:`getKappa <routine-getKappa>`
   - :ref:`getKappa2 <routine-getKappa2>`

