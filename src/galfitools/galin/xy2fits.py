#!/usr/bin/env python3

import numpy as np
import os
import subprocess as sp
import os.path
from astropy.io import fits

import argparse



# code to convert ASCII xy positions to FTIS mask 


def main(): 


    args = ReadParser()

    ImageFile= args.ImageFile
    AsciiFile= args.AsciiMask
    Value = args.val 

    xy2fits().MakeFits(ImageFile, AsciiFile, Value)


def ReadParser():


    #parser argument section

    parser = argparse.ArgumentParser(description="code to convert ASCII x,y positions to FTIS mask")

    parser.add_argument("ImageFile", help="The Image file ")
    parser.add_argument("AsciiMask", help="The ascii file with the x,y positions ")
    parser.add_argument("-c","--val", type=int, help="the value in counts for the masked pixels",default=1)
 
    args = parser.parse_args()


    return args


class xy2fits:


    def MakeFits(self,ImageFile, AsciiFile, Value):
        
        (tmp)=AsciiFile.split(".")

        namefile=tmp[0]

        maskfits=namefile + ".fits"


        (ncol, nrow) = self.GetAxis(ImageFile)

        self.MakeImage(maskfits, ncol, nrow)

        X, Y = np.genfromtxt(AsciiFile, delimiter="", unpack=True)

        X = X.astype(int)
        Y = Y.astype(int)

        X = X-1
        Y = Y-1


        self.PutPix(X,Y,Value,maskfits)

        print("Ascii -> Fits done ")

        return maskfits

    def GetAxis(self,Image):
        "Get number of rows and columns from the image"

        hdu = fits.open(Image)
        ncol = hdu[0].header["NAXIS1"]  # for hubble images
        nrow = hdu[0].header["NAXIS2"]
        hdu.close()

        return ncol, nrow

    def MakeImage(self,newfits, sizex, sizey):
        "create a new blank Image"

        if os.path.isfile(newfits):
            print("{} deleted; a new one is created ".format(newfits))

            runcmd = "rm {}".format(newfits)
            errrm = sp.run([runcmd], shell=True, stdout=sp.PIPE,
                        stderr=sp.PIPE, universal_newlines=True)


        hdu = fits.PrimaryHDU()
        hdu.data = np.zeros((sizey, sizex),dtype=np.float64)
        hdu.writeto(newfits, overwrite=True)

        return True



    def PutPix(self,X,Y,Value,ImageFits):

        # original file
        hdu=fits.open(ImageFits)
        Image = hdu[0].data

        ## for some strange reason I have to interchange X and Y
        Image[[Y],[X]]=Value

        hdu[0].data=Image
        hdu.writeto(ImageFits,overwrite=True)
        hdu.close()




######################



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

