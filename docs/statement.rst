
.. _statement:



---------------------

**Statement of need**
---------------------




The analysis of galaxy morphology through image fitting is a fundamental task 
in extragalactic astronomy. GALFIT [@peng02; @peng10] is a well-established tool that performs parametric 
two-dimensional modeling of galaxy surface brightness profiles. However, GALFIT does 
not provide utilities for related tasks such as generating mask images, selecting 
stars for PSF modeling, estimating initial fit parameters, or interpreting its output 
formats in a structured way. These tasks are typically handled manually or with ad hoc
scripts, which can reduce reproducibility and increase the barrier for new users.

**GALFITools** [@anorve24] addresses this gap by offering a cohesive suite of Python-based tools 
that extend the functionality of GALFIT. It facilitates the full modeling pipeline, 
from input preparation to result interpretation. This includes routines to construct 
PSFs, galaxy modeling via Multi-Gaussian Expansion [MGE, @cappellari02], estimate sky backgrounds,  simulate 
galaxies, estimate initial parameters, generate diagnostic plots, model visualization, 
and compute photometric parameters from multiple Sérsic
components such as effective radius, bar size, Sérsic index,  (among others). 
These functionalities are available as command-line tools, in line with GALFIT's command-line 
interface, but can also be imported as Python modules.

While other packages such as `IMFIT` [@erwin15] and `ProFit` [@robotham16] provide similar modeling capabilities, 
GALFIT remains one of the most widely used tools in the astronomical community for galaxy 
image modeling, supported by an extensive user base and numerous legacy workflows. GALFITools 
is not a replacement for GALFIT, but a complementary toolset specifically tailored to its 
ecosystem. By reducing the technical overhead and providing automation, GALFITools supports 
researchers in conducting large-scale studies of galaxy structure and photometry with 
improved efficiency and reproducibility.


