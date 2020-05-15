
# EllipSectGalfit.py

EllipSectGalfit.py creates surface brightness profiles
for galaxy and galfit model from the galfit output: galfit.XX. 
See peng et al. (2002). It also creates multiple surface brightness 
plots for separate angles.  

**Code: [EllipSectGalfit.py](../EllipSectGalfit.py)**

This is a "quick" substitute for IRAF's ellipse
routine. It creates an Surface brightness profile for the galaxy
and model.

The options to run the code in the terminal (or ipython) are:

```
 ./EllipSectGalfit.py [GALFITOutputFile] [-logx] [-q AxisRatio] [-pa PositionAngle] [-comp] [-pix] [-ranx/y Value] [-grid] [-dpi Value] [-noplot] [-phot] [-sbout] [-noplot] [-minlevel Value] [-sectors Value] [-object Name] [-filter Name] [-snr] [-help] [-checkimg] [-noned] [-distmod Value] [-magcor Value] [-scalekpc Value][-sbdim Value]
 ```

Examples: 

```
 ./EllipSectGalfit.py galfit.01 -logx
```

*or* 

```
 ./EllipSectGalfit.py galfit.02 -q 0.35 -pa 60 -comp
```
### Input

**GALFITOutputFile**: GALFIT output file  (e.g. galfit.01)

## OPTIONS

**help**: Help menu

**logx**: plots X-axis as logarithm

**q**: axis ratio value. If ignored, it takes the one from the last component in GALFITOutputFile.

**pa**: position angle value (same as GALFIT). If ignored, it takes the one 
from the last component in GALFITOutputFile.

**comp**: plots include the individual model components

**pix**: plots the top of x-axis in pixels

**ranx**: constant that multiplies the range of the x axis to increase/decrease *or* 
it can be used as xmin-xmax to change range

**rany**: constant that multiplies the range of the y axis to increase/decrease it
*or* it can be used as ymin-ymax to change range.

**noplot**: the code creates images but do not display them.

**grid**: display a grid in the plot

**dpi**: dots per inch value to increase/decrease resolution.

**sbout**: Creates output file containing the surface brightness profiles.


### Photometry output options

**phot**: Compute photometry. Check the variables created in output file.

The below options are used only if 'phot' is enabled:

**snr**: Creates a signal to noise image. This is created dividing the galaxy image
with the one sigma image created by GALFIT

**object**: used for 'phot' to search in NED. For instance, if you are looking
for photometry data for galaxy m51, then used as "-object m51" the same name 
that you will used to search in NED.

**filter**: used for 'phot' to indicate band for NED. If you need galactic correction
for B  filter then used as "-filter B". Band "R" is the default option.  


Any of the following options disabled the connection to NED

**noned**: avoid to connect to NED. No luminosity nor absolute magnitude is computed.

**distmod**: manual input for Distance Modulus.  

**magcor**: manual input for Galactic Extinction. 

**scalekpc**: manual input for equivalence of ''/kiloparsec. 

**sbdim**: manual input for surface brightness dimming.

### Advanced 

**minlevel**: Parameter given directly to sectors_photometry.
                It stops when it founds this value. Check sectors_photometry manual 

**sectors**: parameter given directly to sectors_photometry. Divide ellipse in 'sectors'
                      Check sectors_photometry manual
                     
**checkimg**: save the images used for sectors_photometry in individual components

## Notes

* You need to install the mgefit library via pip:  

```
pip install mgefit
```


* EllipSectGalfit already uses the mask image (option "*F*" GALFIT) only if this
    is a **FITS** image. In case your mask is an *ASCII* file, you can convert it to **FITS**
    using the [xy2fits.py](xy2fits.md) tool.


* EllipSectGalfit uses axis ratio (*q*) and position angle (*pa*) to create an "ellipse" *grid* 
    using the function *sectors_photometry* from the *mgefit* library. Unlike IRAF's Ellipse, 
    *q* and *pa* are fixed through radius. See the two images below: 

    ![A85 gal](../img/A85.img.png)
    ![A85 mod](../img/A85.img.mod.png)

    For this reason, errors are expected to 
    be greater than those coming from IRAF's ellipse since EllipSectGalfit 
    averages errors for different isophotes, specially for large radius. While, on the other hand, IRAF's ellipse
    can change axis ratio and angular position for each isophote. 
    
    This is how mgefit *sectors_photometry* returns the counts data and, unless I write my own code, I can't change that. 

* Be sure to run this code in the same path that you run GALFIT. 

* Since EllipSectGalfit reads the "*B*" option of galfit.XX file, this
    must be the last GALFIT fit. 

* The angles shown in the multi-plot are measured from the galaxy's major axis.
    They are not measured from the Y-axis as it is the case in GALFIT.

* In order for the program to detect the components, they must share the same 
    center (x,y).

* For the comp option, It could be some small differences between the angle shown 
    in the top right corner and the one from each component. This is because
    *sectors_photometry* is applied different for individual components and the 
    galaxy itself. They are at different angles. To see the real angle which the 
    component is measured check the output file at that angle with the *-sbout* option 

## Examples

See the examples below for an elliptical galaxy that was fitted 
with 7 gaussians (images for this galaxy are displayed above). 

* Example 1: 
    ./EllipSectGalfit.py galfit.46 

    ![A85 ](../img/A85.png)
    ![A85 ](../img/A85.mul.png)

* Example 2:
    ./EllipSectGalfit.py galfit.46 --logx

    ![A85 ](../img/A85.log.png)
    ![A85 ](../img/A85.mul.log.png)

* Example 3 (UPDATE THIS): 
    ./EllipSectGalfit.py galfit.46 --comp
    (displays the 7 gaussians)

    ![A85 ](../img/A85.sub.png)
    ![A85 ](../img/A85.mul.sub.png)

* Example 4: 
    ./EllipSectGalfit.py galfit.46 --pix
    (put pixels marks on top x-axis)

    ![A85 ](../img/A85.pix.png)
    ![A85 ](../img/A85.mul.pix.png)

* Example 5: 
    ./EllipSectGalfit.py galfit.46 --ranx 0.5 
    (range in x-axis is decreased 50%)

    ![A85 ](../img/A85.ranx1.png)
    ![A85 ](../img/A85.mul.ranx1.png)

* Example 6: 
    ./EllipSectGalfit.py galfit.46 --rany 2 
    (range in y-axis doubled)

    ![A85 ](../img/A85.rany1.png)
    ![A85 ](../img/A85.mul.rany1.png)

* Example 7: 
    ./EllipSectGalfit.py galfit.46 --ranx 1-50 
    (x-axis ranges from 1 to 50)

    ![A85 ](../img/A85.ranx2.png)
    ![A85 ](../img/A85.mul.ranx2.png)

* Example 8: 
    ./EllipSectGalfit.py galfit.46 --grid 

    ![A85 ](../img/A85.grid.png)
    ![A85 ](../img/A85.mul.grid.png)


**EllipSectGalfit.py** uses the mgefit library which is
described in Cappellari, MNRAS, 333, 400 (2002).

Check my others GALFIT tools [here](../README.md)
 

