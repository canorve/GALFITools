
Concepts
========

This page provides a short overview of common astronomy and image–analysis
terms used throughout the GALFITools documentation. The goal is to make
the package more accessible for users who may be less familiar with
astronomical jargon.

Surface brightness (SB)
-----------------------
Surface brightness is the flux received per unit area on the sky,
commonly reported in magnitudes per square arcsecond (mag arcsec⁻²).
It describes how light is distributed across a galaxy image.

Point–spread function (PSF)
---------------------------
The point–spread function represents the response of the telescope and
detector to a point source of light (such as a star). It describes the
blurring of the image due to the instrument and atmosphere (for
ground–based observations). GALFIT requires a PSF to deconvolve models.

Sérsic profile and index
------------------------
A Sérsic profile is a mathematical function that describes how the
brightness of a galaxy varies with radius. The **Sérsic index** (*n*)
controls the shape of the profile: low *n* values represent disk–like,
exponential profiles, while high *n* values represent more concentrated,
bulge–like profiles.

Mask
----
A mask is an image that flags which pixels should be ignored during
fitting (for example, bright foreground stars, cosmic rays, or image
defects). GALFITools includes utilities to create masks automatically or
manually.

Initial parameters
------------------
Initial parameters are starting guesses for the fitting procedure, such
as galaxy magnitude, effective radius, axis ratio, and Sérsic index.
Providing reasonable initial values helps GALFIT converge to the correct
solution.

Sky background
--------------
The sky background is the level of light in an image not associated with
the target object. Accurate background estimation is important because
it strongly affects measured magnitudes and profiles.

References
----------
- Peng, C. Y., et al. (2002). *Detailed structural decomposition of
  galaxy images*. AJ, 124, 266.
- Graham, A., & Driver, S. (2005). *A review of galaxy structural
  parameters*. PASA, 22, 118.
