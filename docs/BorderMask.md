
## BorderMask.py

Program that creates a mask that fills the
edges of an HST image.

To run the code just type in the command line:
```
 ./BorderMask.py [ImageFile] [MaskBorderImage] [--val Value] [--le0]
```

Where ImageFile is the original image. MaskBorderImage is the
output mask image. Value is the number that will have the
pixel flux in the border mask. le0 will mask all the pixels below or
equal to zero.
