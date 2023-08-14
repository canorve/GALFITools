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

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8216473.svg
   :target: https://doi.org/10.5281/zenodo.8216473



==========
GALFITools
==========


A library for  `GALFIT`_ 

`GALFIT`_, a well-established two-dimensional image fitting algorithm, 
as outlined in the work by Peng et al. (2002, AJ, 124, 266), 
serves as the basis for precise astronomical image surface brightness 
analysis. In pursuit of optimizing the utilization of GALFIT, GALFITools emerges 
as an collection of Python routines. These routines 
enhances the input and output parsing associated with GALFIT.



.. _GALFIT: https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html


GALFITools extends its utility through an array of functionalities, 
including the facilitation of mask creation, star selection for PSFs, generation 
of multiple initial parameters, simulate galaxy images, multigaussian 
expansion (MGE) fitting, as well as computation of sky background 
and other pertinent photometric variables.




--------------------------------

**Installation Instructions**
-------------------------------



The python libraries required are:

-  numpy
-  astropy
-  scipy
-  matplotlib
-  mgefit



Install `GALFIT`_ if you haven't done so.


*Optionally*, the establishment of a virtual environment can be considered.


The latest release of GALFITools can be downloaded 
and subsequently installed via one of the following methods:


::

   cd GALFITools 
   pip install . 

or

::

   cd GALFITools 
   python setup.py install


Alternatively, you can install it via pip:


::

   pip install GALFITools 




In conjunction with the installation, a compilation of pertinent 
shell `commands <docs/howto.rst>`__ will be incorporated. Subsequently, a comprehensive 
evaluation of GALFITools' performance can be conducted through 
automated tests using the following procedure:

To run the tests locally, install and invoke tox:

::
   
   pip install tox


run the tests:

::

    tox 



-----------------

**HOW TO USE**
~~~~~~~~~~~~~~~~~

For comprehensive insights into GALFITools' repertoire 
of routines and their optimal deployment, it is 
recommended to consult the provided documentation on 
usage, accessible via the following link: `How to use <docs/howto.rst>`__.
This comprehensive resource elaborates on the 
practical implementation of each individual routine.

---------------------

**API Reference**
~~~~~~~~~~~~~~~~~~~~


For the customization of these functions 
to align with your specific scripting requirements, 
the API documentation serves as an indispensable resource. 
The detailed instructions for utilizing these 
functions within your own scripts can be found here: 

`API <docs/api.rst>`__


--------------

**License**
--------------

The codebase of GALFITools is governed by the terms of the **MIT** license.


-----------

**Cite as**
-----------

To acknowledge the utility of GALFITools in your research, 
we kindly request that you cite it as follows:

Añorve, Christopher. (2023). canorve/GALFITools: 
GALFITools v0.15.2 (v0.15.2). Zenodo. https://doi.org/10.5281/zenodo.8216473



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



--------------

**Questions?**
--------------


Should any questions or innovative suggestions arise, 
we encourage you to reach out via email to  canorve [at] gmail [dot] com  
or initiate a discussion by opening an  `issue <https://github.com/canorve/GALFITools/issues>`__.  
Your input is invaluable in fostering the continual refinement of 
GALFITools, for the betterment of the *GALFIT* community and beyond.




.. _pyscaffold-notes:


====
Note
====

This project has been set up using PyScaffold 4.2.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.

