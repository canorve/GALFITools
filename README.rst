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


`GALFIT`_  
is a 2 dimensional image fitting algorithm.  
See Peng et al. (2002, AJ, 124, 266) for general description. 

.. _GALFIT: https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html




GALFITools is a collection of Python
scripts that improve the input and 
output analysis of GALFIT.


This is a set of routines which help the GALFIT's  user to make 
masks, cut stars, create multiple initial parameters, simulate 
galaxies, fit MGE, compute sky and compute other photometric variables.


--------------

**Installation**
----------------

The python libraries required are:

-  numpy
-  astropy
-  scipy
-  matplotlib
-  mgefit


Install GALFIT if you haven't done so.

Download the latest release, and installed it via

::

   cd GALFITools 
   pip install . 

or

::

   cd GALFITools 
   python setup.py install


Also, you can install it via pip:

::

   pip install GALFITools 


Run the automated tests:

::

    tox 


A list of commands will be installed and 
run it via shell.


--------------

**HOW TO USE**
~~~~~~~~~~~~~~

To learn how to use every command routine of GALFITools,
check:


`How to use <docs/howto.rst>`__

Where you can find full explanations of each routine.

--------------

**API**
~~~~~~~~~~~~~~

Check how to use the functions for your own scripts:

`API <docs/api.rst>`__



--------------

**License**
--------------

The code is under the license of **MIT**


-----------

**Cite as**
-----------

If you find this code useful, please cite as:

Añorve, Christopher. (2023). canorve/GALFITools: 
GALFITools v0.15.2 (v0.15.2). Zenodo. https://doi.org/10.5281/zenodo.8216473



---------------

**Other Stuff**
---------------

Check EllipSect to create surface brightness profiles
from GALFIT output and estimate other photometric parameters:

`here <https://github.com/canorve/EllipSect>`__


.. _pyscaffold-notes:


--------------

**Questions?**
--------------

Do you have any questions or suggestions? Please send an email to
canorve [at] gmail [dot] com or open an
`issue <https://github.com/canorve/GALFITools/issues>`__

I’m open to new ideas that can benefit the library *GALFITools* and the
*GALFIT* community




====
Note
====

This project has been set up using PyScaffold 4.2.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.

