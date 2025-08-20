
.. _statement:



============
Statement of need
============


GALFITools addresses a key gap in astronomical image analysis: while
GALFIT is a powerful 2D fitting engine, its input/output management,
mask creation, PSF star selection, and sky background estimation
require repetitive manual work. GALFITools provides Python routines
to streamline these tasks, making the fitting workflow reproducible
and efficient for both research and teaching.


GALFIT is a widely used tool for two–dimensional surface–brightness
modeling of galaxies. However, preparing inputs and handling outputs
from GALFIT often involves repetitive manual steps, such as mask
construction, PSF star selection, generation of initial parameters,
and organization of fitting results. These tasks can become
time–consuming and error–prone, particularly in large surveys.

GALFITools provides a set of Python routines that automate and
streamline these tasks. By integrating common pre– and post–processing
steps into a reproducible workflow, GALFITools lowers the technical
barriers for new users and increases efficiency for experienced
researchers. The package is designed for astronomers who analyze
galaxy images and require flexible, scriptable tools that complement
GALFIT in both research and teaching contexts.


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


