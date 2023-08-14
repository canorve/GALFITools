

import pytest
import os


import subprocess as sp



from galfitools.mge.mge2galfit import mge2gal
from galfitools.mge.SbProf import sbProf



def test_mge2gal():

    class mge_options:
        pass

    mgeargs = mge_options()


 

    mgeargs.image = "ngc3344J.mge2galfit.fits"
    mgeargs.Ds9regFile ="ds9.mge2galfit.reg"
    path="tests/"
    mgeargs.image=path+mgeargs.image
    mgeargs.Ds9regFile=path+mgeargs.Ds9regFile
 
    mgeargs.magzpt = 20.0820
    mgeargs.twist = False 
    mgeargs.regu = False 
    mgeargs.center = False 
    mgeargs.mask = None 

    mgeargs.psf = 0
    mgeargs.sky = 0 
    mgeargs.plate = 1 
    mgeargs.gauss = False 

    mgeargs.psfile = None 
    mgeargs.sigfile = None 

    mgeargs.freeser = False
    mgeargs.freesky = False 

    mgeargs.regu = False 
    mgeargs.numgauss = None


    mge2gal(mgeargs) 



    parfile="mseGALFIT.txt"
    file2  = "mgegas.txt"
    namefile="sectors.png"
    nameout="tests/ngc3344J.mge2galfit.png"
    namecons="constraints.txt"

    assert os.path.isfile(parfile)
    if os.path.isfile(parfile):
        os.remove(parfile)

    assert os.path.isfile(file2)
    if os.path.isfile(file2):
        os.remove(file2)

    assert os.path.isfile(namefile)
    if os.path.isfile(namefile):
        os.remove(namefile)

    assert os.path.isfile(nameout)
    if os.path.isfile(nameout):
        os.remove(nameout)

    assert os.path.isfile(namecons)
    if os.path.isfile(namecons):
        os.remove(namecons)




    return None




def test_sbProf():

    class sb_options:
        pass


    sbargs = sb_options()

    sbargs.Image  =  "A1656-1-3.sbProf.fits"      
    sbargs.Ds9Region = "ds9.sbProf.reg"
    path="tests/"
    sbargs.Image=path+sbargs.Image
    sbargs.Ds9Region=path+sbargs.Ds9Region
 


    sbargs.mgzpt = 21.817
    sbargs.mask = "mask.sbProf.fits"
    
    sbargs.mask = path+sbargs.mask

    sbargs.axrat = None 
    sbargs.angle = None 
    sbargs.sky = 376 
    sbargs.plate = 0.68 
    sbargs.output = "sb.png"
    sbargs.output= path+sbargs.output
    sbargs.center = False  
    sbargs.ranx = None 
    sbargs.rany = None 
    sbargs.logx = False 
    sbargs.pix = False  
    sbargs.grid = False 

    sbargs.rad = None 
    sbargs.rad2 = None 



    sbProf(sbargs)


    assert os.path.isfile(sbargs.output)

    if os.path.isfile(sbargs.output):
        os.remove(sbargs.output)


    namefile = "gal.png"
    #namefile = path+namefile

    assert os.path.isfile(namefile)

    if os.path.isfile(namefile):
        os.remove(namefile)



    return None


#####################
#####################
#####################
#####################


