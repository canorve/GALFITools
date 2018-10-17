#! /usr/bin/env python3

# programa que crea una mascara que rellena los bordes con Value
# para imagenes del HST.

import numpy as np
from astropy.io import fits
#import matplotlib.pyplot as plt
import sys
import os
import subprocess as sp



def main():

  if len(sys.argv[1:]) != 3:
    print ('Missing arguments')
    print ("Usage:\n %s [ImageFile] [MaskBorderImage] [Value] " % (sys.argv[0]))
    print ("Example:\n %s image.fits border.fits 100" % (sys.argv[0]))
    sys.exit()


  ImageFile= sys.argv[1]
  BorderFile= sys.argv[2]
  Value = sys.argv[3]



  Value=np.float64(Value)

  BorderMask(ImageFile,BorderFile,Value)

#################################################################
#################################################################

def BorderMask(ImageFile,BorderFile,Value):

  SexFile="sexmask.fits"


  if not os.path.exists(ImageFile):
    print ('%s: image filename does not exist!' %(sys.argv[1]))
    sys.exit()


  (ncol, nrow) = GetAxis(ImageFile)


  MakeImage(BorderFile, ncol, nrow)

  MakeImage1(SexFile, ncol, nrow)


#  BorderImage = np.zeros((ncol, nrow),dtype=np.float64)


  # original file
  hdu=fits.open(ImageFile)
  Image = hdu[0].data
#  Header = hdu[0].header


# output file
  hdu2=fits.open(BorderFile)
  BorderImage = hdu2[0].data

# output file 2
  hdu3=fits.open(SexFile)
  SexImage = hdu3[0].data



  mask = Image <= 0


  if mask.any():

    BorderImage[mask] = Value
    SexImage[mask] = 0


  hdu2[0].data=BorderImage
  hdu3[0].data=SexImage

#  hdu[0].header=Header


  hdu2.writeto(BorderFile,overwrite=True)
  hdu3.writeto(SexFile,overwrite=True)


  hdu.close()
  hdu2.close()
  hdu3.close()




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



def MakeImage1(newfits, sizex, sizey):
    "create a new blank Image"
# k Check

    if os.path.isfile(newfits):
        print("{} deleted; a new one is created \n".format(newfits))

        runcmd = "rm {}".format(newfits)
        errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


    hdu = fits.PrimaryHDU()
    hdu.data = np.ones((sizey, sizex),dtype=np.float64)
    hdu.writeto(newfits, overwrite=True)

    return True





#end of program
if __name__ == '__main__':
  main()
