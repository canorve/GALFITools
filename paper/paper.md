---
title: 'GALFITools: A Python library to enhance GALFIT usage in galaxy image modeling'
tags:
  - Python
  - astronomy
  - galaxy modeling
  - GALFIT
  - photometry
authors:
  - name: Christopher Añorve
    orcid:  0000-0002-3721-8869
    affiliation: "1"
    corresponding: true
affiliations:
  - name: Facultad de Ciencias de la Tierra y el Espacio, Universidad Autónoma de Sinaloa, México
    index: 1
date: 1 May 2025
bibliography: paper.bib
---

# Summary

Understanding how galaxies form and evolve requires measuring their light distributions 
in images taken by telescopes. This process often involves fitting mathematical models 
to galaxy images to extract properties such as size, brightness, components, and shape. 
GALFIT is a widely used tool for this purpose, but it requires careful 
preparation of input files and interpretation of results, which can be a barrier to efficient use.

GALFITools is a Python library that streamlines this workflow by automating many 
of the tasks surrounding the use of GALFIT. These include generating image masks, 
estimating sky background, modeling the telescope’s point spread function (PSF), and 
extracting physical parameters from GALFIT outputs. The software is designed for 
researchers and students who work with galaxy image modeling and aims to make the 
process more reproducible, accessible, and scalable.


# Statement of need

The analysis of galaxy morphology through image fitting is a fundamental task 
in extragalactic astronomy. GALFIT [@peng02; @peng10] is a well-established tool that performs parametric 
two-dimensional modeling of galaxy surface brightness profiles. However, GALFIT does 
not provide utilities for related tasks such as generating mask images, selecting 
stars for PSF modeling, estimating initial fit parameters, or interpreting its output 
formats in a structured way. These tasks are typically handled manually or with ad-hoc 
scripts, which reduces reproducibility, becomes time-consuming and error-prone, 
raises the barrier for new users, and makes the application to large surveys difficult.

GALFITools [@anorve24] addresses this gap by offering a cohesive suite of Python-based tools 
that extend the functionality of GALFIT. It facilitates the full modeling pipeline, 
from input preparation to result interpretation. This includes routines to construct 
PSFs, model galaxies via Multi-Gaussian Expansion [MGE, @cappellari02], estimate sky 
backgrounds, simulate galaxies, estimate initial parameters, generate diagnostic plots, 
visualise models, and compute photometric parameters from multiple Sérsic
components such as effective radius, bar size and Sérsic index (among others). 
These functionalities are available as command-line tools, in line with GALFIT's command-line 
interface, but can also be imported as Python modules.

GALFITools lowers the technical barriers for new users and increases efficiency 
for experienced researchers. The package is designed for astronomers who analyze 
galaxy images and require flexible, scriptable tools that complement GALFIT analysis.
While other packages such as `Imfit` [@erwin15] and `ProFit` [@robotham16] provide 
similar modeling capabilities, GALFIT remains one of the most widely used tools 
in the astronomical community for galaxy image modeling, supported by an extensive 
user base and numerous legacy workflows. GALFITools 
is not a replacement for GALFIT, but a complementary toolset specifically tailored to its 
ecosystem. By reducing the technical overhead and providing automation, GALFITools supports 
researchers in conducting large-scale studies of galaxy structure and photometry with 
improved efficiency and reproducibility.


# Acknowledgements

The author acknowledges support from Universidad Autónoma de Sinaloa through 
project PROFAPI 2022, with the project key A1009. GALFITools was developed using 
PyScaffold. GALFIT itself is maintained by Chien Peng. 

# References

