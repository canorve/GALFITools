#!/usr/bin/env python3

import numpy as np
import sys
import os
import stat
import subprocess as sp
import os.path
from astropy.io import fits
import scipy
import scipy.special


def main():

    if len(sys.argv[1:]) == 0:
      print ('Missing arguments')
      print ("Usage: %s [ImageFile] [AsciiMask] [--val Value] " % (sys.argv[0]))
      print ("Example: %s A85.fits mask.txt " % (sys.argv[0]))

      sys.exit()

    ImageFile= sys.argv[1]
    AsciiFile= sys.argv[2]
    Value = 1

    flagvalue=False

    OptionHandleList = ["--val"]
    options = {}
    for OptionHandle in OptionHandleList:
        options[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)] if OptionHandle in sys.argv else None
    if options['val'] != None:
        flagvalue=True

################################
    if flagvalue == True:
        opt={}
        OptionHandle="--val"
        opt[OptionHandle[2:]] = sys.argv[sys.argv.index(OptionHandle)+1]
        Value=np.int(opt['val'])


    print("Value for mask pixels:  ",Value)

######################

    (tmp)=AsciiFile.split(".")

    namefile=tmp[0]

    maskfits=namefile + ".fits"


    (ncol, nrow) = GetAxis(ImageFile)

    MakeImage(maskfits, ncol, nrow)


    X, Y  = np.genfromtxt(AsciiFile, delimiter="", unpack=True)

    X = X.astype(int)
    Y = Y.astype(int)


    PutPix(X,Y,Value,maskfits)



#############################################################################
######################### End of program  ###################################
#     ______________________________________________________________________
#    /___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/___/_/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/|
#   |___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|__/|
#   |_|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|___|/
##############################################################################

def PutPix(X,Y,Value,ImageFits):


      # original file
    hdu=fits.open(ImageFits)
    Image = hdu[0].data


## for some strange reason I have to interchange X and Y
    Image[[Y],[X]]=Value


    hdu[0].data=Image

    hdu.writeto(ImageFits,overwrite=True)

    hdu.close()







def GetAxis(Image):
    # k Check
    "Get number of rows and columns from the image"

    hdu = fits.open(Image)
    ncol = hdu[0].header["NAXIS1"]  # for hubble images
    nrow = hdu[0].header["NAXIS2"]
    hdu.close()
    return ncol, nrow




def MakeImage(newfits, sizex, sizey):
    "create a new blank Image"
# k Check

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((sizey, sizex),dtype=np.float64)
    hdu.writeto(newfits, overwrite=True)

    return True



if __name__ == '__main__':
    main()
