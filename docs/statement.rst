
.. _statement:



====================
Statement of need
====================


The analysis of galaxy morphology through image fitting is a fundamental task 
in extragalactic astronomy. GALFIT (Peng et al., 2002; Peng et al. 2010) is a 
well-established tool that performs parametric two-dimensional modeling of galaxy 
surface brightness profiles. However, GALFIT does not provide utilities for 
related tasks such as generating mask images, selecting stars for PSF modeling, 
estimating initial fit parameters, or interpreting its output formats in a 
structured way. These tasks are typically handled manually or with ad hoc scripts, 
which reduces reproducibility, becomes time-consuming and error-prone, 
raises the barrier for new users, and makes the application to large surveys difficult.


**GALFITools** (Añorve 2024) addresses this gap by offering a cohesive suite of 
Python-based tools that extend the functionality of GALFIT. It facilitates the full 
modeling pipeline, from input preparation to result interpretation. This includes 
routines to construct PSFs, galaxy modeling via Multi-Gaussian Expansion (MGE, 
Cappellari, 2002), estimate sky backgrounds,  simulate galaxies, estimate initial 
parameters, generate diagnostic plots, model visualization, and compute photometric 
parameters from multiple Sérsic components such as effective radius, bar size, Sérsic 
index,  (among others). These functionalities are available as command-line tools, in 
line with GALFIT's command-line interface, but can also be imported as Python modules.


GALFITools lowers the technical barriers for new users and increases efficiency 
for experienced researchers. The package is designed for astronomers who analyze 
galaxy images and require flexible, scriptable tools that complement GALFIT analysis.




