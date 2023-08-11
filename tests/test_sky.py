

import pytest
import os


import subprocess as sp



from galfitools.sky.GalfitSky import galfitSky

from galfitools.sky.Sky import sky



def test_galfitSky():

    imgname =   "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits" 

    mgzpt =  21.817
    scale = 0.68 

    X = 169 
    Y = 544 

    initsky = args.initsky


    galfitSky(imgname, maskfile, mgzpt, scale, X, Y, initsky)



    assert os.path.isfile("sky.txt")

    if os.path.isfile("sky.txt"):
        os.remove("sky.txt")





    return None


def test_sky():



    imgname =   "A1656-1-3.sbProf.fits"
    maskfile = "mask.sbProf.fits" 
    filereg = "ds9.sky.reg"

    mean, sig = sky(imgname, maskimage, filereg)

    print("Sky within 3 sigma:") 

    print("mean sky: {:.3f} ".format(mean))
    print("std sky: {:.3f} ".format(sig))

    tol = 1e-3

    result1 =  369.717 
    result2 = 4.407 



    diffsky1 = abs(mean -result1)

    assert diffsky1 < tol

    diffsky2 = abs(sig -result2)

    assert diffsky2 < tol



    return None



#####################
#####################
#####################
#####################
#####################



# simple run test
def test_exit():
    with pytest.raises(SystemExit) as e:
        ellipsectors.run()
    assert e.type == SystemExit 
    assert e.value.code == 2 


#check if galfit is installed
def test_galfit():

    runcmd="galfit -help"
    err = sp.run([runcmd],shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)  

    assert err.returncode == 0, "is GALFIT installed?"




# checking the creation of files
def test_files():

    arg=['tests/galfit.01', '-np']

    path="tests/"

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-01-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)



    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)



    assert os.path.isfile(filepng)
    assert os.path.isfile(filemulpng)


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)




# checking the creation of the components file
def test_comp():

    arg=['tests/galfit.01','--comp', '--noplot']

    path="tests/"

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-01-mul.png"

    filepng = path+filepng
    filemulpng = path+filemulpng



    filecomp = "imgblock-01-comp.fits"

    filecomp= path+filecomp

    if os.path.isfile(filecomp):
        os.remove(filecomp)


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    assert os.path.isfile(filecomp),"is GALFIT installed?"


    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)


    if os.path.isfile(filecomp):
        os.remove(filecomp)






# checking the creation of the sbout files
def test_phot():


    arg=['tests/galfit.01','--phot', '--noned', '--noplot']

    path="tests/"

    filephot = "imgblock-01-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-01-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-01-mul.png"

    filesig = "imgblock-01-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-01-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-01-cub.png"
    filecube = path+filecube



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)




    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 1e-3

    bt= 1.000 
    tidal = 0.688 
    lchinu = 1.038
    bump = 0.064 
    snr = 2.871 
    std_snr = 5.324 
    aic= 629.555 
    bic = 660.334
    effrad3 = 1.192
    effrad = 2.623 
    effrad9 = 14.914

    diffbt =abs(bt-photapi.BulgeToTotal ) 
    difftidal = abs(tidal-photapi.tidal)
    difflchinu = abs(lchinu-photapi.objchinu)
    diffbump = abs(bump-photapi.bump)
    diffsnr = abs(snr-photapi.snr)
    diffstd_snr = abs(std_snr-photapi.stdsnr)
    diffaic= abs(aic-photapi.AICrit)
    diffbic = abs(bic-photapi.BICrit)



    diffeffrad = abs(effrad - photapi.EffRad)
    diffeffrad9 = abs(effrad9 - photapi.EffRad9) 
    diffeffrad3 = abs(effrad3 - photapi.EffRad3) 

    assert os.path.isfile(filephot)

    assert diffbt < tol
    assert difftidal < tol
    assert difflchinu < tol
    assert diffbump < tol
    assert diffsnr < tol
    assert diffstd_snr < tol
    assert diffaic < tol
    assert diffbic < tol

    assert diffeffrad < tol
    assert diffeffrad9 < tol
    assert diffeffrad3 < tol


    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)


# checking the computation of sky with the gradient method
def test_gsky():


    arg=['tests/galfit.01','-gsky', '-ri','2', '-skw','3','--noplot']

    path="tests/"

    filephot = "imgblock-01-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-01-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-01-mul.png"

    filesig = "imgblock-01-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-01-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-01-cub.png"
    filecube = path+filecube

    filersky = "imgblock-01-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-01-ringmask.fits"
    filermsky = path + filermsky



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 1e-2

    sky = 0.63


    diffsky = abs(sky - photapi.gradskymean)



    assert diffsky < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)

 # checking the computation of sky with the random method
def test_rsky():


    arg=['tests/galfit.01','-rsky','--noplot']

    path="tests/"

    filephot = "imgblock-01-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-01-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock.01-png"
    filemulpng = "imgblock-01-mul.png"

    filesig = "imgblock-01-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-01-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-01-cub.png"
    filecube = path+filecube

    filersky = "imgblock-01-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-01-ringmask.fits"
    filermsky = path + filermsky



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 


    # read user's input 
    args = ArgParsing(arg)

    # full program:
    photapi = SectorsGalfit(args)




    #tolerance parameter
    tol = 5

    sky = 0.9


    diffsky = abs(sky - photapi.randskymean)



    assert diffsky < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)



# checking the number of free params 
def test_freepar():


    arg=['tests/galfit.01','--noplot']

    path="tests/"

    filephot = "imgblock-01-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-01-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-mul.png"

    filesig = "imgblock-01-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-01-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-01-cub.png"
    filecube = path+filecube

    filersky = "imgblock-01-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-01-ringmask.fits"
    filermsky = path + filermsky



    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)
 

    from ellipsect.inout.galfit  import numParFree
    from ellipsect.inout.galfit  import Galfit 
    from ellipsect.sectors.sect  import PassArgs 
    from ellipsect.inout.galfit import SelectGal 

    # read user's input 
    args = ArgParsing(arg)

    # full program:
    #photapi = SectorsGalfit(args)

    ellconf = PassArgs(args) # from now on, ellconf is used instead of args

    ######################################
    ####### Read Galfit File #############

    galfit = Galfit(ellconf.galfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    galcomps = SelectGal(galcomps,ellconf.distmax,ellconf.numcomp)

    freepar = numParFree(galcomps)




    #tolerance parameter
    tol = 1e-2 

    totpar = 7 

    diffpar = abs(totpar - freepar)


    assert diffpar < tol



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)


# checking the number of free params 
def test_img_size():


    arg=['tests/galfit.01','--noplot']

    path="tests/"

    filephot = "imgblock-01-out.txt"
    filephot= path+filephot

    filecomp = "imgblock-01-comp.fits"
    filecomp= path+filecomp

    filepng = "imgblock-01.png"
    filemulpng = "imgblock-01-mul.png"

    filesig = "imgblock-01-sig.fits"
    filesig= path+filesig

    filecheck = "imgblock-01-check.fits"
    filecheck= path+filecheck

    filecube = "imgblock-01-cub.png"
    filecube = path+filecube

    filersky = "imgblock-01-ring.fits"
    filersky = path + filersky

    filermsky = "imgblock-01-ringmask.fits"
    filermsky = path + filermsky


    tempmask="tempmask.fits"
    filetempmask = tempmask 

    filepng = path+filepng
    filemulpng = path+filemulpng

    if os.path.isfile(filepng):
        os.remove(filepng)
    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)

    if os.path.isfile(filetempmask):
        os.remove(filetempmask)
 


    from ellipsect.inout.galfit  import Galfit 
    from ellipsect.sectors.sect  import PassArgs 
    from ellipsect.inout.galfit import SelectGal 

    # read user's input 
    args = ArgParsing(arg)

    # full program:
    #photapi = SectorsGalfit(args)

    ellconf = PassArgs(args) # from now on, ellconf is used instead of args

    ######################################
    ####### Read Galfit File #############

    galfit = Galfit(ellconf.galfile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    galcomps = SelectGal(galcomps,ellconf.distmax,ellconf.numcomp)



    from ellipsect.inout.read import ReadGALFITout 

    ReadGALFITout(ellconf, galhead, galcomps) 

    from ellipsect.inout.read import prefixNames
    #creates names of the output files based on prefix of galfit output
    prefixNames(ellconf, galhead.outimage)


    from ellipsect.inout.galfit  import readDataImg
    dataimg = readDataImg(ellconf, galhead)


    

    assert dataimg.img.shape == dataimg.model.shape 
    assert dataimg.img.shape == dataimg.mask.shape 



    if os.path.isfile(filepng):
        os.remove(filepng)

    if os.path.isfile(filemulpng):
        os.remove(filemulpng)

    if os.path.isfile(filecomp):
        os.remove(filecomp)

    if os.path.isfile(filephot):
        os.remove(filephot)

    if os.path.isfile(filesig):
        os.remove(filesig)

    if os.path.isfile(filecheck):
        os.remove(filecheck)

    if os.path.isfile(filecube):
        os.remove(filecube)

    if os.path.isfile(filersky):
        os.remove(filersky)

    if os.path.isfile(filermsky):
        os.remove(filermsky)

    if os.path.isfile(filetempmask):
        os.remove(filetempmask)
 






