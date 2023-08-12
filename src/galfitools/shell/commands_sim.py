

import argparse

from galfitools.sim.MakeSim import makeSim

from galfitools.shell.prt import printWelcome

def mainMakeSim():

    printWelcome()

    parser = argparse.ArgumentParser(description="simulates a observed galaxy from a GALFIT model")

    parser.add_argument("image", help="the GALFIT galaxy model")
    parser.add_argument("newimage", help="the name of the new galaxy image")


    parser.add_argument("-s","--sky", type=float, help="the sky background value. default = 0", default=0)
    parser.add_argument("-std","--std", type=float, help="the sky standard deviation. default = 1", default=1)

    parser.add_argument("-g","--gain", type=float, help="the gain value of the image. default =  1", default=1)



    args = parser.parse_args()

    image = args.image 
    GAIN = args.gain 

    skymean = args.sky
    skystd = args.std 

    newimage = args.newimage 




    makeSim(image, GAIN, skymean, skystd, newimage)






