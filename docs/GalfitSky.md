## GalfitSky.py

GalfitSky.py computes the sky using GALFIT

Maskbox.py put (or remove) mask patch on
an already existing mask FITS file. It does this
using a box file region. This file is
created by DS9 program. You can create as many
box regions as you like in a single file. Maskbox will fill
those box regions with the number specified in Value.   

To run the code just type in the command line (or in ipython):

```
./GalfitSky.py [ImageFile] [Magzpt] [scale] [--X x] [--Y y]
```
Where ImageFile is the galaxy image. Magzpt is the
zero point of the image. Scale is the factor which is
increased the ellipse mask of the galaxy. X and Y are
the pixel coordinates of the galaxy to compute its sky.  
