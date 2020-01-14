
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
