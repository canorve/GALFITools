#! /usr/bin/python3


import numpy as np
import argparse
import glob

from pathlib import Path

# code to make init galfit files 



class InitGal:


    def __init__(self, GalfitFile: str, number: int, param3: float, param4: float, param5: float, param6: float, param7: float, param8: float, param9: float, param10: float, numcomp: int):

        self.GalfitFile = GalfitFile 
        self.number = number
        self.param3 = param3
        self.param4 = param4
        self.param5 = param5
        self.param6 = param6
        self.param7 = param7

        self.param8 = param8
        self.param9 = param9
        self.param10 = param10

        self.numcomp = numcomp

        self.MakeParams()
        self.MakeBash()




    def MakeParams(self):   



        lines = self.ReadFile(self.GalfitFile) #read  galfit input file

        (tmp)=self.GalfitFile.split(".")
        prefixname=tmp[0]


        for idx, item in enumerate(np.arange(1,self.number+1)): # cicle for every galfit file

            self.namefile = prefixname + "-" + str(item) + ".gal"

            fileout = open(self.namefile, "w")

            flagsky3 = False # flag to avoid to write in sky c3

            index = 0
            comp = 0


            while index < len(lines):

                line = lines[index]
                (tmp) = line.split()


                if (tmp[0] == "0)"):
                    comp = comp + 1

                if self.numcomp:

                    if self.numcomp == comp:

                        if (tmp[0] == "3)") and (flagsky3  == False):     # input image
                            if self.param3:
                                c3min, c3max = self.param3[0],self.param3[1]
                                c3 = np.random.default_rng().uniform(c3min,c3max,1)
                                tmp[1] = str(np.round(c3[0],2))
                                line= " ".join(tmp)
                                line= " " + line 

                        if tmp[0] == "4)":     # input image
                            if self.param4:
                                c4min, c4max = self.param4[0],self.param4[1]
                                c4 = np.random.default_rng().uniform(c4min,c4max,1)
                                tmp[1] = str(np.round(c4[0],2))
                                line= " ".join(tmp)
                                line= " " + line 


                        if tmp[0] == "5)":     # input image
                            if self.param5:
                                c5min, c5max = self.param5[0],self.param5[1]
                                c5 = np.random.default_rng().uniform(c5min,c5max,1)
                                tmp[1] = str(np.round(c5[0],2))
                                line= " ".join(tmp)
                                line= " " + line 

                               
                        if tmp[0] == "6)":     # input image
                            if self.param6:
                                c6min, c6max = self.param6[0],self.param6[1]
                                c6 = np.random.default_rng().uniform(c6min,c6max,1)
                                tmp[1] = str(np.round(c6[0],2))
                                line= " ".join(tmp)
                                line= " " + line 


                        if tmp[0] == "7)":     # input image
                            if self.param7:
                                c7min, c7max = self.param7[0],self.param7[1]
                                c7 = np.random.default_rng().uniform(c7min,c7max,1)
                                tmp[1] = str(np.round(c7[0],2))
                                line= " ".join(tmp)
                                line= " " + line 

                        if tmp[0] == "8)":     # input image
                            if self.param8:
                                c8min, c8max = self.param8[0],self.param8[1]
                                c8 = np.random.default_rng().uniform(c8min,c8max,1)
                                tmp[1] = str(np.round(c8[0],2))
                                line= " ".join(tmp)
                                line= " " + line 

                        if tmp[0] == "9)":     # input image
                            if self.param9:
                                c9min, c9max = self.param9[0],self.param9[1]
                                c9 = np.random.default_rng().uniform(c9min,c9max,1)
                                tmp[1] = str(np.round(c9[0],2))
                                line= " ".join(tmp)
                                line= " " + line 

                        if tmp[0] == "10)":     # input image
                            if self.param10:
                                c10min, c10max = self.param10[0],self.param10[1]
                                c10 = np.random.default_rng().uniform(c10min,c10max,1)
                                tmp[1] = str(np.round(c10[0],2))
                                line= " ".join(tmp)

                else:

                    if (tmp[0] == "3)") and (flagsky3  == False):     # input image
                        if self.param3:
                            c3min, c3max = self.param3[0],self.param3[1]
                            c3 = np.random.default_rng().uniform(c3min,c3max,1)
                            tmp[1] = str(np.round(c3[0],2))
                            line= " ".join(tmp)
                            line= " " + line 

                    if tmp[0] == "4)":     # input image
                        if self.param4:
                            c4min, c4max = self.param4[0],self.param4[1]
                            c4 = np.random.default_rng().uniform(c4min,c4max,1)
                            tmp[1] = str(np.round(c4[0],2))
                            line= " ".join(tmp)
                            line= " " + line 


                    if tmp[0] == "5)":     # input image
                        if self.param5:
                            c5min, c5max = self.param5[0],self.param5[1]
                            c5 = np.random.default_rng().uniform(c5min,c5max,1)
                            tmp[1] = str(np.round(c5[0],2))
                            line= " ".join(tmp)
                            line= " " + line 

                           
                    if tmp[0] == "6)":     # input image
                        if self.param6:
                            c6min, c6max = self.param6[0],self.param6[1]
                            c6 = np.random.default_rng().uniform(c6min,c6max,1)
                            tmp[1] = str(np.round(c6[0],2))
                            line= " ".join(tmp)
                            line= " " + line 


                    if tmp[0] == "7)":     # input image
                        if self.param7:
                            c7min, c7max = self.param7[0],self.param7[1]
                            c7 = np.random.default_rng().uniform(c7min,c7max,1)
                            tmp[1] = str(np.round(c7[0],2))
                            line= " ".join(tmp)
                            line= " " + line 

                    if tmp[0] == "8)":     # input image
                        if self.param8:
                            c8min, c8max = self.param8[0],self.param8[1]
                            c8 = np.random.default_rng().uniform(c8min,c8max,1)
                            tmp[1] = str(np.round(c8[0],2))
                            line= " ".join(tmp)
                            line= " " + line 

                    if tmp[0] == "9)":     # input image
                        if self.param9:
                            c9min, c9max = self.param9[0],self.param9[1]
                            c9 = np.random.default_rng().uniform(c9min,c9max,1)
                            tmp[1] = str(np.round(c9[0],2))
                            line= " ".join(tmp)
                            line= " " + line 

                    if tmp[0] == "10)":     # input image
                        if self.param10:
                            c10min, c10max = self.param10[0],self.param10[1]
                            c10 = np.random.default_rng().uniform(c10min,c10max,1)
                            tmp[1] = str(np.round(c10[0],2))
                            line= " ".join(tmp)




                lineout = line + "\n"
                fileout.write(lineout)

                if (tmp[0] == "0)") and (tmp[1] == "sky") :     # avoids to write in sky component
                    flagsky3 = True # sky component must be at the end of the file 
                                    # otherwise it will fail

                index += 1
            fileout.close()


#############################################################

    def MakeBash(self):

        #creating bash script

        listfiles = glob.glob("*.gal")

        bashout = open("rungalfit.sh", "w")

        lineout="#!/bin/bash \n" 
        bashout.write(lineout)

        for gfile in listfiles:

            lineout="galfit " + gfile  + "\n" 
            bashout.write(lineout)


        bashout.close()

        p = Path('rungalfit.sh')
        p.chmod(0o755)




#############################################################

    def ReadFile(self,inputf):

        File = open(inputf,"r")

        # All lines including the blank ones
        lines = (line.rstrip() for line in File)
        #lines = (line.split('#', 1)[0] for line in lines)  # remove comments
        #lines = (line.rstrip() for line in lines) # remove lines containing only comments
        lines = (line for line in lines if line)  # Non-blank lines

        lines = list(lines)

        return lines




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
  mainInitGal()


