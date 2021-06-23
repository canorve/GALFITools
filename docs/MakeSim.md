
## MakeSim.py

MakeSim creates an artificial galaxy image.
Make your surface brightness model with GALFIT. Run
GALFIT with option "p) 1" and introduce it as a ImageFile.

Add GAIN, sky mean, and skystd as inputs to add to your simulate
galaxy. 

To run the code just type in the command line (or in ipython):

```
 ./MakeSim.py [ImageFile] [GAIN] [skymean] [skystd] [newimage]
```
Where ImageFile is the input ImageFile. GAIN is the
gain of the image, skymean and skystd is the mean and 
standard deviation of the sky and newimage is the name 
of the output image. 




