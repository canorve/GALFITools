
## MakeMask.py

MakeMask.py creates a mask from a Sextractor catalog.

To run the code just type in the command line (or in ipython):

```
./MakeMask.py [SexFile] [ImageFile] [MaskFileOut] [scale]
```

SexFile is the sextractor catalog.  ImageFile is the image
file you want to use the mask. MaskFileOut is the mask image
that MakeMask.py will create. Scale is used to increase/decrease 
the size of the ellipses mask for each object in the mask image.

Check the sextractor files in this directory. They will help you
to create the sextractor file this program needs.
