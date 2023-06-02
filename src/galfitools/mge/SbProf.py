#! /usr/bin/env python3

import numpy as np
import sys
import os
import subprocess as sp
from astropy.io import fits
import os.path  
import scipy
import scipy.special
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import mimetypes
import warnings



from astropy.io.votable import parse
from mgefit.sectors_photometry import sectors_photometry
from mgefit.find_galaxy import find_galaxy

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,NullFormatter,
                               AutoMinorLocator,LogLocator,LinearLocator,AutoLocator)



import argparse


def mainSbProf():

    parser = argparse.ArgumentParser(description = "SbProf: creates a surface brightness profile from a ellipse ds9 region")


    # required arguments
    parser.add_argument("Image", help = "image fits file")
    parser.add_argument("Ds9Region", help = "Ds9 ellipse region file")

    parser.add_argument("-mz","--mgzpt", type=float, help="Magnitud zero point", default=25)
    parser.add_argument("-m","--mask", type=float, help="mask fits file" )

    parser.add_argument("-s","--sky", type=float, help="sky value", default=0)
    parser.add_argument("-p","--plate", type=float, help="plate scale ", default=1)
    parser.add_argument("-o","--output", type=str, help="output file", default="sb.png")



    args = parser.parse_args()


    image = args.image
    ds9reg = args.Ds9Region
    mgzpt = args.mgzpt
    mask =  args.mask
    sky = args.sky
    plate = args.plate
    output = args.output

    SbProf(image, ds9reg, mgzpt, mask, sky, plate, output)

    
    print('Done')



def SbProf(image, ds9reg, mgzpt, mask, sky, plate, output):
    'creates the surface brightness profile'




if __name__ == '__main__':

    mainSbProf()
