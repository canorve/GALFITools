.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/GALFITools.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/GALFITools
    .. image:: https://readthedocs.org/projects/GALFITools/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://GALFITools.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/GALFITools/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/GALFITools
    .. image:: https://img.shields.io/conda/vn/conda-forge/GALFITools.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/GALFITools
    .. image:: https://pepy.tech/badge/GALFITools/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/GALFITools
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/GALFITools

.. image:: https://img.shields.io/pypi/v/GALFITools.svg
    :alt: PyPI-Server
    :target: https://pypi.org/project/GALFITools/

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8216472.svg 
  :target: https://doi.org/10.5281/zenodo.8216472 

.. image:: https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg 
   :target: CODE_OF_CONDUCT.md 


============
GALFITOOLS
============

`GALFIT`_, a well-established two-dimensional image fitting algorithm, 
as outlined in the work by Peng et al. (2002, AJ, 124, 266), 
serves as the basis for precise astronomical image surface brightness 
analysis. In pursuit of optimizing the utilization of GALFIT, GALFITools emerges 
as an collection of Python routines. These routines 
enhances the input and output parsing associated with GALFIT.


GALFITools extends its utility through an array of functionalities, 
including the facilitation of mask creation, star selection for PSFs, generation 
of multiple initial parameters, simulate galaxy images, multigaussian 
expansion (MGE) fitting, as well as computation of sky background 
and other pertinent photometric variables.


.. _GALFIT: https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html



-----------------------------------

**INSTALLATION INSTRUCTIONS**
-----------------------------------

First of all, install `GALFIT`_ if you haven't done so. Check
instructions `here <https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html`__.
Make sure that GALFIT can run in any path in your *linux/macOS* terminal.

*Optionally*, consider setting up a virtual environment.

Install GALFITools via pip:


::

   pip install GALFITools 


Here is video tutorial for the installation of GALFITools:


.. image:: https://img.youtube.com/vi/rqZLxR1yRCs/maxresdefault.jpg
    :alt: IMAGE ALT TEXT HERE
    :target: https://www.youtube.com/watch?v=rqZLxR1yRCs



--------------------------------

**GET STARTED WITH GALFITOOLS**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once installed, you can test GALFITools using any of 
the commands in the console. 

In the example below, command *getReComp* is used with the output GALFIT file 
`galfit.01 <img/galfit.01>`__  to calculate the effective radius for a composite 
model consisting of three Sersic components:

::
   
  ❯ getReComp galfit.01
  GALFITools: a library for GALFIT
  Version: 1.14.1
  webpage: https://github.com/canorve/GALFITools

  number of model components:  3
  Using a theta value of : 14.13 degrees 

  Total Magnitude of the galaxy: 9.79 

  Surface brightness at radius of 50% of light (μr): 21.15 mag/'' 

  Mean Surface Brightness at effective radius (<μ>e): 20.27 mag/'' 

  The radius at 50% of light is 199.64 pixels or 49.91 " 


-----------------

**DOCUMENTATION**
~~~~~~~~~~~~~~~~~

full documentation along usage guide can be found in `https://galfitools.readthedocs.io/en/latest/ <https://galfitools.readthedocs.io/en/latest/>`__.




------------------------------

**CONTRIBUTING AND SUPPORT**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us know if you have issues, questions or innovative suggestions arise, 
we encourage you to reach out via email to  canorve [at] gmail [dot] com  
or initiate a discussion by opening an  `issue <https://github.com/canorve/GALFITools/issues>`__.  
Your input is invaluable in fostering the continual refinement of 
GALFITools, for the betterment of the *GALFIT* community and beyond.


- Issue tracker: `https://github.com/canorve/GALFITools/issues <https://github.com/canorve/GALFITools/issues>`__. 

- Source code: `https://github.com/canorve/GALFITools <https://github.com/canorve/GALFITools>`__. 

- `CODE OF CONDUCT <CODE_OF_CONDUCT.md>`__


--------------

**LICENSE**
--------------

The codebase of GALFITools is governed by the terms of the **MIT** license.


---------------------

**API REFERENCE**
~~~~~~~~~~~~~~~~~~~~


For the customization of these functions 
to align with your specific scripting requirements, 
the API documentation serves as an indispensable resource. 
The detailed instructions for utilizing these 
functions within your own scripts can be found here: 

`API <docs/api.rst>`__


----------------------

**CITING GALFITOOLS**
-----------------------

If you find this code useful in your research, 
we kindly request that you cite it as follows:

Añorve, C. (2024). canorve/GALFITools: GALFITools v1.15.0 (v1.15.0). 
Zenodo. https://doi.org/10.5281/zenodo.13994492


cite all versions using the DOI: https://doi.org/10.5281/zenodo.8216472




--------------------------

**Additional Resources**
--------------------------

Check EllipSect to create surface brightness profiles
from GALFIT output and estimate other photometric parameters:

For further capabilities and valuable extensions 
pertaining to GALFIT output, such as the generation of 
surface brightness profiles and estimation of other 
photometric parameters, we invite you to explore the 
EllipSect tool: 

`EllipSect <https://github.com/canorve/EllipSect>`__




.. _pyscaffold-notes:


====
Note
====

This project has been set up using PyScaffold 4.2.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.



