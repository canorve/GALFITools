## xy2fits.py

Program that creates a fits mask from a ascii mask.

To run the code just type in the command line:

```
Usage: ./xy2fits.py [ImageFile] [AsciiMask] [--val Value]
```

Where ImageFile is the original image. AsciiMask is the
GALFIT Ascii Mask that only contains the positions (X,Y)
of the bad pixels. Value is the flux that contains the
mask pixels (Default = 1).
