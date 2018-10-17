#! /usr/bin/env python3


# programa que crea un parche en una seccion de una imagen FITS



import numpy as np
from astropy.io import fits
#import matplotlib.pyplot as plt

#import pyfits as pf
import sys
import os


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
  sys.exit()

if not os.path.exists(RegFile):
  print ('%s: reg filename does not exist!' %(sys.argv[2]))
  sys.exit()


hdu=fits.open(ImageFile)


naxis1=hdu[0].header["NAXIS1"]
naxis2=hdu[0].header["NAXIS2"]

Image = hdu[0].data


v1=[]
v2=[]
v3=[]
v4=[]


# Region file format: DS9 version 4.1
#global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
#image
#box(168.5,1195.5,43,63,0)

#v1,v2,v3,v4 = np.genfromtxt(RegFile,skip_header=3,delimiter=",",usecols = (0,1,2,3),dtype="S10",unpack=True)

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






hdu[0].data=TempImage

hdu.writeto(ImageFile,overwrite=True) # clobber sobreescribe




hdu.close()
