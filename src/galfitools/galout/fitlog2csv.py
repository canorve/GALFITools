#! /usr/bin/env python3 


import numpy as np

import re


def log2csv(num:int, fileout: str):


    file = "fit.log"

    #fileout="fitlog.csv"

    fitlogF = open(file,"r")

    f_out = open(fileout, "w")

    # All lines including the blank ones
    lines = (line.rstrip() for line in fitlogF)
    lines = (line.split('#', 1)[0] for line in lines)  # remove comments
    # remove lines containing only comments
    lines = (line.rstrip() for line in lines)
    lines = (line for line in lines if line)  # Non-blank lines

    lines = list(lines)
    index = 0

    linnp = np.array(lines)
    idxfits = np.where(linnp == "-----------------------------------------------------------------------------")



    if num:
        #use user number 
        num=4 #number to be readed from argpaser
        index = idxfits[0][2*num-2]
        cnt=num*2-1
    else:
        #use last one by default
        index = idxfits[0][-2]
        cnt=int(len(idxfits[0])-1)



    line = lines[index]
    (tmp) = line.split()

    if tmp[0] == "-----------------------------------------------------------------------------":

        index += 1
        lineout = lines[index]+"\n"
        f_out.write(lineout)

        index += 1
        lineout = lines[index]+"\n"
        f_out.write(lineout)

        index += 1
        lineout = lines[index]+"\n"
        f_out.write(lineout)


        index += 1
        lineout = lines[index]+"\n"
        f_out.write(lineout)

        #put chi and sky at the beggining

        index2=idxfits[0][cnt] - 4


        line= lines[index2]
        line=re.sub('[^A-Za-z0-9.+.-]+', ' ', line)
        (tmp) = line.split()
        lineout1 = "sky = "+tmp[3]+", "+tmp[4]+", "+tmp[5]

        index2 += 1
        line= lines[index2]
        line=re.sub('[^A-Za-z0-9.+.-]+', ' ', line)
        (tmp) = line.split()


        lineout2 = " skyerr = "+tmp[0]+", "+tmp[1]+", "+tmp[2]

        lineout = lineout1 + "," + lineout2+"\n"

        f_out.write(lineout)


        index2 += 1



        line= lines[index2]
        line=re.sub('[^A-Za-z0-9.+]+', ' ', line)
        (tmp) = line.split()

        lineout1 = "Chi^2 = " + tmp[2] + " ndof = " + tmp[4]

        index2 += 1
        line = lines[index2]
        (tmp) = line.split()

        lineout2 = " ".join(tmp) 

        lineout = lineout1 + " " + lineout2+"\n"

        f_out.write(lineout)


        cnt+=2

        lineout="model parameters, error parameters\n" 
        f_out.write(lineout)
 

        lineout="------------------------------------------------------------\n" 
        f_out.write(lineout)
 
        
        #end of header

        index += 1
        line = lines[index]
        line=re.sub('[^A-Za-z0-9.+.-]+', ' ', line)
        (tmp) = line.split()

        while tmp[0] != "sky":
 
                
            lineout1 = ",".join(tmp) 

            index += 1
            line= lines[index]
            line=re.sub('[^A-Za-z0-9.+.-]+', ' ', line)
            (tmp) = line.split()


            lineout2 = ",".join(tmp) 

            lineout = lineout1 + "," + lineout2+"\n"

            f_out.write(lineout)


            index += 1
            line= lines[index]
            line=re.sub('[^A-Za-z0-9.+.-]+', ' ', line)
            (tmp) = line.split()



    fitlogF.close()
    f_out.close()


