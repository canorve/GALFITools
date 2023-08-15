

import pytest
import os


import subprocess as sp



from galfitools.sim.MakeSim import makeSim


def test_makeSim():

    image = "M81-BD.fits" 
    path="tests/"
    image = path+image
    GAIN = 6.5 

    skymean = 0 
    skystd = 1  

    newimage = "M81Sim.fits" 
    newimage= path+newimage



    makeSim(image, GAIN, skymean, skystd, newimage)


    assert os.path.isfile("poissonimg.fits")

    if os.path.isfile("poissonimg.fits"):
        os.remove("poissonimg.fits")

    assert os.path.isfile("skynoise.fits")

    if os.path.isfile("skynoise.fits"):
        os.remove("skynoise.fits")




    assert os.path.isfile(newimage)

    if os.path.isfile(newimage):
        os.remove(newimage)



    return None



#check if galfit is installed
def test_galfit():

    runcmd="galfit -help"
    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  

    assert err.returncode == 0, "is GALFIT installed?"



