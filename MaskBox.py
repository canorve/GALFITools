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
import matplotlib.pyplot as plt

# programa que hace parches (mascaras) en una seccion de una imagen FITS

def main(): 

  if len(sys.argv[1:]) != 3:
    print ('Missing arguments')
    print ("Usage:\n %s [ImageFile] [RegFile] [Value] " % (sys.argv[0]))
    print ("Example:\n %s image.fits box.reg 100" % (sys.argv[0]))
    sys.exit()

  ImageFile= sys.argv[1]
  RegFile= sys.argv[2]
  Value = sys.argv[3]

  if not os.path.exists(ImageFile):
    print ('%s: image filename does not exist!' %(sys.argv[1]))
    print ('Creating a new image file ')

    hdu=fits.PrimaryHDU()
    nx = input("enter numbers of pixels in X ") 
    ny = input("enter numbers of pixels in Y ") 

    nx = np.int64(nx)
    ny = np.int64(ny)
    Image = np.zeros([ny,nx])
    hdu.data=Image
    hdu.writeto(ImageFile,overwrite=True) 

  hdu=fits.open(ImageFile)
  #  naxis1=hdu[0].header["NAXIS1"]
  #  naxis2=hdu[0].header["NAXIS2"]

  Image = hdu[0].data



  if not os.path.exists(RegFile):
    print ('%s: reg filename does not exist!' %(sys.argv[2]))
    sys.exit()

  v1=[]
  v2=[]
  v3=[]
  v4=[]


  f1 = open(RegFile,'r')

  lines = f1.readlines()

  f1.close()

  for line in lines:

    b1 =line.split("(")

    p =line.split(",")


    if b1[0] == "box":
      x1=p[0]
      x2=x1[4:]


      v1.append(np.float(x2))
      v2.append(np.float(p[1]))
      v3.append(np.float(p[2]))
      v4.append(np.float(p[3]))


  xpos=np.array(v1)
  ypos=np.array(v2)
  xlong=np.array(v3)
  ylong=np.array(v4)


  xlo = xpos - xlong/2.
  xhi = xpos + xlong/2.

  ylo = ypos - ylong/2.
  yhi = ypos + ylong/2.


  xlo=np.round(xlo)
  xhi=np.round(xhi)

  ylo=np.round(ylo)
  yhi=np.round(yhi)

  xlo=xlo.astype(int)
  xhi=xhi.astype(int)

  ylo=ylo.astype(int)
  yhi=yhi.astype(int)


  TempImage=Image


  for i in range(len(xpos)):
    TempImage[ylo[i]-1:yhi[i],xlo[i]-1:xhi[i]]=Value



  hdu.data=TempImage

  hdu.writeto(ImageFile,overwrite=True) 


  hdu.close()

  
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
if __name__ == '__main__':
  main()


