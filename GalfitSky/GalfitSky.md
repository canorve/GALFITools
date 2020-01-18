## GalfitSky.py

GalfitSky.py computes the sky using GALFIT. 

To run the code just type in the command line (or in ipython):

```
./GalfitSky.py [ImageFile] [Magzpt] [scale] [--X x] [--Y y]
```

Where ImageFile is the galaxy image. Magzpt is the
zero point of the image. Scale is the factor which is
increased the ellipse mask of the galaxy. X and Y are
the pixel coordinates of the galaxy center to compute its sky. 
After run it, follow instructions displayed by GalfitSky.py

This program needs some of the Sextractor files that comes
in this directory. Run it along with those.

