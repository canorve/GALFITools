
## EllipSectGalfit.py

This is a "quick" (and dirty) substitute for IRAF's ellipse
routine. It creates a "ellipse" profile of galaxy
and GALFIT's model output (peng et al 2002). It also create a
multiple plots for different angles.  

![A85 gal](img/A85mgetwist-gal.png)
![A85 gal](img/A85mgetwist-mod.png)
![A85 gal](img/A85mgetwist.png)
![A85 gal](img/A85mgetwist-mul.png)

To run the code just type in the command line (or in ipython):

```
 ./EllipSectGalfit.py [GALFITOutputFile] [--logx] [--q AxisRatio] [--pa PositionAngle] [--sub] [--pix] [--ranx Value]
 ```

GALFITOutputFile: GALFIT output file  (e.g. galfit.01)
logx: plots X-axis as logarithm

q: introduce axis ratio  where the user introduces the axial symmetric generalized
ellipse where EllipSectGalfit computes the profile. If ignored,
it takes the one of the last component that comes in GALFITOutputFile.

pa: introduce position angle (same as GALFIT).  If ignored,
this code takes the one of the last component that comes in GALFITOutputFile.

sub: plots subcomponents.

pix: plot the x-axis in pixels

ranx: constant that multiplies the range of the x axis

rany: constant that multiplies the range of the y axis

noplot: do not display images

Example:
 EllipSectGalfit.py galfit.01 --logx

or Example:
 EllipSectGalfit.py galfit.02 --q 0.35 --pa 60 --sub


Note: Errors given by EllipSectGalfit will be greater
than those coming from IRAF's ellipse since this program leaves fixed
the axis ratio for the whole galaxy while IRAF's ellipse
can change axis ratio for each isophote.

Note 2: EllipSectGalfit already uses your mask (option "F)" GALFIT) if this
is a FITS image.  


#### Warning
Be sure to run this code in the same path that you run GALFIT.


EllipSectGalfit.py uses the mgefit library which is
described in Cappellari, MNRAS, 333, 400 (2002).
