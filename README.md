# GALFITools


GALFITools is a collection of
programs written in python which helps
the GALFIT's user to make masks, and aid
fitting analysis.

## Installation

Copy or clone this code. These codes are
written for python 3.

The python libraries need to be installed are:
- numpy
- sys
- os
- subprocess
- astropy
- scipy
- matplotlib
- mgefit

The main programs are:
- BorderMask.py
- BorderMask2.py
- ElipSectGalfit.py
- MaskBox.py
- GetK.pl

## BorderMask.py

Program that creates a mask that fills the
edges of an HST image.

To run the code just type in the command line:
```
 ./BorderMask.py [ImageFile] [MaskBorderImage] [Value]
```
Where ImageFile is the original image. MaskBorderImage is the 
output mask image. Value is the number that will have the 
pixel flux in the border mask.



## BorderMask2

This is another version of BorderMask.py. The difference
with the previous one is that BorderMask create the mask for those
regions which pixel flux are below or equal to zero.
On the other hand, BorderMask2 creates the mask on those 
regions where the pixels are equal to zero.

To run the code just type in the command line:
```
 ./BorderMask2.py [ImageFile] [MaskBorderImage] [Value]
```

## ElipSectGalfit.py

This is a "quick" (and dirty) substitute for IRAF's ellipse 
routine. It creates a "ellipse" profile of galaxy 
and GALFIT's model output (peng et al 2002). It also create a 
multiple plots for different angles.  

To run the code just type in the command line (or in ipython):

```
 ./EllipSectGalfit.py [GALFITOutputFile] [AxisRatio OPTIONAL] [PA OPTIONAL]
 ```

Where GALFITOutputFile is the output GALFIT file (e.g. galfit.01). 
AxisRatio is an optional argument where
the user can fix the axis ratio of the axial symmetric generalized 
ellipse where EllipSectGalfit computes the profile. If ignored, 
this code takes the one that comes in GALFITOutputFile.  

Note: Errors given by EllipSectGalfit will be greater
than those from IRAF's ellipse since this program leaves fixed
the axis ratio for the whole galaxy while IRAF's ellipse
can change axis ratio for each isophote.


#### Warning
Be sure to run this code in the same path as you run GALFIT.


## Maskbox.py

Maskbox.py put (or remove) mask patch on
an already existing mask FITS file. It does this
using a box file region. This file is
created by DS9 program. You can create as many 
box regions as you like in a single file. Maskbox will fill
those box regions with the number specified in Value.   

To run the code just type in the command line (or in ipython):

```
./MaskBox.py [ImageFile] [RegFile] [Value]
```
Where ImageFile is the mask image to edit. RegFile
is the DS9 Region file which contains the box region
to patch, and, Value is the flux number you want to
establish for the pixels within the box region.



## GetK.pl

Perl script that gives the K constant for a
determined Sersic index. The one that allows the surface
brightness (Ie) to be at the half of the light radius (Re). 
See Sersic equation.

To run the code just type in the command line:
```
 ./usage: ./GetK.pl [SersicIndex]
```
Where SersicIndex is (obviously) the Sersic index.




#### Questions?
Do you have questions or suggestions?
Please send an email to canorve [at] gmail [dot] com

## License
These codes are under the license of **GNU**

ElipSectGalfit.py uses the mgefit library which is
described in Cappellari, MNRAS, 333, 400 (2002).
