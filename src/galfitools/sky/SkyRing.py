#! /usr/bin/env python

import numpy as np 

from  galfitools.mge.mge2galfit import Ds9ell2Kronell
from  galfitools.mge.mge2galfit import GetInfoEllip
from  galfitools.mge.mge2galfit import GetSize
from  galfitools.mge.mge2galfit import GetPmax
from  galfitools.mge.mge2galfit import GetAxis

from astropy.io import fits


def SkyRing(image, mask, ds9regfile, width, center):
    """Computes the sky using gradient over rings """


    ##end input
    obj, xpos, ypos, rx, ry, angle = GetInfoEllip(ds9regfile)
    xx, yy, Rkron, theta, eps = Ds9ell2Kronell(xpos,ypos,rx,ry,angle)

    q = 1 - eps

    hdu = fits.open(image)
    datimg = hdu[0].data
    hdu.close()
 
    hdumask = fits.open(mask)
    maskimg = hdumask[0].data
    hdumask.close()



    (ncol, nrow) = GetAxis(image)

    if center:
        print('center of ds9 ellipse region will be used')
        xpeak, ypeak = xpos, ypos
    else:        
        (xmin, xmax, ymin, ymax) = GetSize(xx, yy, Rkron, theta, eps, ncol, nrow)
        xpeak, ypeak = GetPmax(datimg, maskimg, xmin, xmax, ymin, ymax)
 


    mean, std, median, rad = SkyCal().GetEllipSky(datimg,maskimg,xpeak,ypeak,
                                                        theta,q,Rkron,width,
                                                        "ring.fits","ringmask.fits")



    return mean, std, median, rad



class SkyCal:
    "computes sky background "

    def GetEllipSky(self, ImageFile, MaskFile, xx, yy, thetadeg, q, Rinit, width,namering,ringmask,outliers=True):
        "Gradient sky method"

        self.xx = xx 
        self.yy = yy

        self.thetadeg = 90 + thetadeg
        self.q = q
        self.e = (1 - self.q)
        self.Rinit = Rinit
        self.width = width 

        self.NumRings = 5  # number of rings per loop # read this from function?

        self.ringmask = ringmask 

        self.outliers = outliers

        ###

        #hdumask = fits.open(MaskFile)
        #self.maskimg = hdumask[0].data
        #hdumask.close()
        self.maskimg = MaskFile.copy()



        #hdu = fits.open(ImageFile)
        #self.img = hdu[0].data
        #hdu.close()
        self.img = ImageFile.copy()

        ####
       
        (self.nrow,self.ncol) = self.img.shape


        xmin,xmax,ymin,ymax,Rkron = self.GetXYRBorder()

        self.R = Rkron

        if self.Rinit > Rkron: # avoid radius greater than the border

            self.Rinit= Rkron/3 #is this number ok?
            print("Rinit is greater than image size")
            print("using Rinit = {:0.2f} ".format(self.Rinit))
            print("change this value with --radinit ")





        (xmin, xmax, ymin, ymax) = self.GetSize(self.xx, self.yy, Rkron, self.thetadeg, self.q, self.ncol, self.nrow) # obtain corners of R

        theta = self.thetadeg * np.pi / 180  # Rads!!!

        patch = self.maskimg[ymin - 1:ymax, xmin - 1:xmax] # logical patch mask image

        self.invpatch=np.logical_not(patch)


        Rings=np.arange(self.Rinit,self.R,self.width) # anillos de tamaño width

        sky = np.array([])
        skymed = np.array([])
        skystd = np.array([])
        radius = np.array([])


        count = 0
        idx=0


        ########################################
        #creating ring masks file

        masksky = np.zeros(np.shape(self.img))

        val = 1 
        for ridx, ritem in enumerate(Rings):


            bring = (Rings[ridx] + self.width) * self.q
            aring = Rings[ridx] + self.width 

            #incresing number of points per rad 
            points = 2*np.pi * np.sqrt(0.5*(aring**2 + bring**2)) #Aprox. 
            points = 2*points #doubling the number of points
            points = int(round(points))
           
            alpha = np.linspace(0,2*np.pi,points)
            for tidx, item in enumerate(range(self.width)):

                bim=(Rings[ridx]+item)*self.q

                tempxell = self.xx + (Rings[ridx]+item) * np.cos(alpha) * np.cos(theta) - bim * \
                  np.sin(alpha) * np.sin(theta)

                tempyell = self.yy + (Rings[ridx]+item) * np.cos(alpha) * np.sin(theta) + bim * \
                  np.sin(alpha) * np.cos(theta)

                tempxell = tempxell.round().astype("int")
                tempyell = tempyell.round().astype("int")

                tempxell,tempyell = self.CorSize(tempxell,tempyell)
                masksky[tempyell - 1, tempxell - 1] = val
                
            val += 1 


        #hdu[0].data=masksky
        #hdu.writeto(self.ringmask,overwrite=True) 

        ########################################
        ########################################



        #computing sky in every ring
        for ind, item in enumerate(Rings):


            maskring,idx=self.GetRingMask(masksky[ymin - 1:ymax, xmin - 1:xmax],idx)


            flatimg=self.img[ymin - 1:ymax, xmin - 1:xmax][maskring].flatten()  
            flatimg.sort()

            tot=len(flatimg)

            top=round(.8*tot)
            bot=round(.2*tot)


            if self.outliers:   # eliminate top 80% and bottom 20%
                imgpatch=flatimg[bot:top]
            else:
                imgpatch=flatimg

            mean=np.mean(imgpatch)
            std=np.std(imgpatch)

            median=np.median(imgpatch)

            sky=np.append(sky,mean)
            skymed=np.append(skymed,median)
            skystd=np.append(skystd,std)
            radius=np.append(radius,Rings[idx] + self.width/2)

            print("Ring = {}; rad = {:.2f}; sky mean = {:.2f}; sky std = {:.2f}; median: {:.2f} ".format(idx+1,Rings[idx] + self.width/2,mean,std, median))


            # calcular gradiente
            if (count >= self.NumRings):

                # [1:-1] avoiding the first and last element for gradient 
                gradmask = np.gradient(sky[1:-1]) >= 0 
               
                count = 0
                tempidx=np.where(np.gradient(sky[1:-1]) >= 0) 
                
                if (sky[1:-1][gradmask].any()): 
                    
                    savidx=tempidx[0][0]
                    maskring,none =self.GetRingMask(masksky[ymin - 1:ymax, xmin - 1:xmax],savidx)


                    print("sky computed in ring {} ".format(savidx+2))
                    
                    print("Ring radius = {:.2f} marked in {} ".format(radius[1:-1][savidx],namering))
                    print("the counts value within ring represent the long axis") 
                    self.img[ymin - 1:ymax, xmin - 1:xmax][maskring] = radius[1:-1][savidx] 
                    break

            count += 1
            idx +=1

            if idx == (len(Rings)-1): 
                print("The edge of image has been reached. Sky can not be computed")
                return 0,0,0,0


        #hdu[0].data=self.img
        #hdu.writeto(namering,overwrite=True) 

        finmean,finmedian,finstd,finRad = sky[1:-1][gradmask],skymed[1:-1][gradmask],skystd[1:-1][gradmask],radius[1:-1][gradmask]


        return finmean[0],finstd[0],finmedian[0],finRad[0]


    def GetRingMask(self,masksky,idx):
        ''' obtains the ring selected by index idx'''

        ring= idx + 1

        maskring = masksky == ring 

        maskring=maskring*self.invpatch

        ringcont=0

        while( not(maskring.any()) and (ringcont < 10)):

            if (ringcont == 0):
                    print("Selecting next ring ")

            idx += 1
            ring= idx + 1

            maskring = masksky == ring 

            maskring=maskring*self.invpatch


            ringcont+=1 # avoid eternal loop

        if (ringcont == 10):
            print("max. iteration reached. I couldn't find a ring") 
            return 0,0 # It couldn't found any ring ending 
 




        return maskring,idx


    def CorSize(self,xell,yell):
        '''Correct size for image borders'''
        masksx = xell < 1
        masksy = yell < 1

        xell[masksx] = 1
        yell[masksy] = 1

        masksx = xell > self.ncol
        masksy = yell > self.nrow 

        xell[masksx] = self.ncol 
        yell[masksy] = self.nrow

        return xell,yell




    def GetXYRBorder(self):
        "this subroutine get the coordinates of the border"

        q =  self.q

        theta = self.thetadeg * (np.pi / 180)  # rads!!

        thetax=np.sqrt((np.cos(theta))**2 + (q**2)*(np.sin(theta))**2 )
        thetay=np.sqrt((q**2)*(np.cos(theta))**2 + (np.sin(theta))**2 )


        if (self.thetadeg >-45 and self.thetadeg <= 45):

            xmax=self.ncol
            xmin =1
            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >45 and self.thetadeg <= 135):
            ymax=self.nrow
            ymin =1
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

        elif (self.thetadeg >135 and self.thetadeg <= 225):
            xmax=1
            xmin =self.ncol

            R1 = (xmax - self.xx)/thetax
            R2 = (self.xx - xmin)/thetax

            if (np.abs(R1)>=np.abs(R2)):
                R = R1
            else:
                R = R2

     
        elif (self.thetadeg >225 and self.thetadeg <= 315):
            ymax=1
            ymin =self.nrow
            R1 = (ymax - self.yy)/thetay
            R2 = (self.yy - ymin)/thetay

            if (np.abs(R1)>=np.abs(R2)):
               R = R1
            else:
               R = R2

        R = np.abs(R) #avoids negative numbers
        bim = q * R

        # getting size

        xmin = self.xx - np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        xmax = self.xx + np.sqrt((R**2) * (np.cos(theta))**2 +
                         (bim**2) * (np.sin(theta))**2)

        ymin = self.yy - np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        ymax = self.yy + np.sqrt((R**2) * (np.sin(theta))**2 +
                         (bim**2) * (np.cos(theta))**2)

        mask = xmin < 1
        if mask.any():
            if isinstance(xmin,np.ndarray):
                xmin[mask] = 1
            else:
                xmin = 1

        mask = xmax > self.ncol
        if mask.any():
            if isinstance(xmax,np.ndarray):
                xmax[mask] = self.ncol
            else:
                xmax = self.ncol

        mask = ymin < 1
        if mask.any():
            if isinstance(ymin,np.ndarray):
                ymin[mask] = 1
            else:
                ymin = 1

        mask = ymax > self.nrow
        if mask.any():
            if isinstance(ymax,np.ndarray):
                ymax[mask] = self.nrow
            else:
                ymax =self.nrow

        xmin=np.int32(np.round(xmin))
        ymin=np.int32(np.round(ymin))
        xmax=np.int32(np.round(xmax))
        ymax=np.int32(np.round(ymax))


        return (xmin,xmax,ymin,ymax,np.int32(R))


    def GetSize(self,x, y, R, theta, q, ncol, nrow):
        '''this subroutine get the maximun
        and minimim pixels for Kron and sky ellipse'''

        bim = q * R

        theta = theta * (np.pi / 180)  # rads!!

        # getting size
        xmin = x - np.sqrt((R**2) * (np.cos(theta))**2 +
                           (bim**2) * (np.sin(theta))**2)

        xmax = x + np.sqrt((R**2) * (np.cos(theta))**2 +
                           (bim**2) * (np.sin(theta))**2)

        ymin = y - np.sqrt((R**2) * (np.sin(theta))**2 +
                           (bim**2) * (np.cos(theta))**2)

        ymax = y + np.sqrt((R**2) * (np.sin(theta))**2 +
                           (bim**2) * (np.cos(theta))**2)

        mask = xmin < 1
        if mask.any():
            if isinstance(xmin,np.ndarray):
                xmin[mask] = 1
            else:
                xmin = 1

        mask = xmax > ncol

        if mask.any():
            if isinstance(xmax,np.ndarray):
                xmax[mask] = ncol
            else:
                xmax = ncol

        mask = ymin < 1
        if mask.any():
            if isinstance(ymin,np.ndarray):
                ymin[mask] = 1
            else:
                ymin = 1

        mask = ymax > nrow
        if mask.any():
            if isinstance(ymax,np.ndarray):
                ymax[mask] = nrow
            else:
                ymax = nrow

        xmin=np.int32(np.round(xmin))
        ymin=np.int32(np.round(ymin))
        xmax=np.int32(np.round(xmax))
        ymax=np.int32(np.round(ymax))


        return (xmin, xmax, ymin, ymax)


    ##### End of sky class ##################
    #########################################
    #########################################

