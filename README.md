# GALFITools


GALFITools is a collection of
programs written in python which helps
the GALFIT's user to make masks, and aid
the fitting analysis.

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


The main programs are:
- BorderMask.py
- BorderMask2.py
- ElipSectGalfit.py
- MaskBox.py


## BorderMask.py

Program that makes a mask that fills the
edges of an HST image.

To run the code just type in the command line:
```
 ./BorderMask.py [ImageFile] [MaskBorderImage] [Value]
```
Where ImageFile is the original image. MaskBorderImage is the output mask image. Value is the number that will have the the pixel flux in the border mask.



## BorderMask2

This is another version of BorderMask.py. The difference
is that BorderMask create the mask for those
regions which pixel flux are below or equal to zero.
On the other hand, BorderMask2 creates the mask on those regions where the pixels are equal to zero.

To run the code just type in the command line:
```
 ./BorderMask2.py [ImageFile] [MaskBorderImage] [Value]
```

## ElipSectGalfit.py

This is a "quick" (and dirty) substitute for ellipse routine of IRAF. It creates a "ellipse" profile of the galaxy and model of the GALFIT (peng et al 2002) output. It also create a multiple plots for different angles.  

To run the code just type in the command line (or in ipython):

```
 ./EllipSectGalfit.py [GALFITOutputFile] [AxisRatio OPTIONAL]
 ```

Where GALFITOutputFile is the output GALFIT file (e.g. galfit.01). AxisRatio is an optional argument where
the user can fix the axis ratio of the axial symmetric generalized ellipse where EllipSectGalfit computes the profile. If ignored, this code takes the one that
comes in GALFITOutputFile.  

Note: Errors given by EllipSectGalfit will be greater
than those from IRAF's ellipse since this program fix
the axis ratio for the whole galaxy while IRAF ellipse
can change for each isophote.


#### Warning
Be sure to run this code in the same path as you run GALFIT.


## Maskbox.py

Maskbox.py put (or remove) mask patch on
an already existing mask FITS file. It does this
using a box file region. This file is
made by DS9 program. You can create as many box regions
as you like in a single file. Maskbox will fill
those box regions with the number specified in Value.   

To run the code just type in the command line (or in ipython):

```
./MaskBox.py [ImageFile] [RegFile] [Value]
```
Where ImageFile is the mask image to edit. RegFile
is the DS9 Region file which contains the box region
to patch, and, Value is the flux number you want to
establish for the pixels within the box region.



#### Questions?
Do you have questions or suggestions?
Please send an email to canorve [at] gmail [dot] com

## License
This code is under the license of **GNU**

If you use ElipSectGalfit.py  for your
publications you should cite cappellari et al. (2006) in your paper since I wrote this code using their libraries
