#! /usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from astropy.io import fits
from scipy import stats
from matplotlib.patches import Ellipse

from astropy import units as U

import argparse




def displayCube(cubeimage: str, namecube: str, dpival: int, brightness: float, contrast:float, 
                    cmap: str, scale: float, noplot: bool) -> None:


    ######################
    #shows the image cube#
    ######################{


    ell=[]


    ShowCube(cubeimage, namepng = namecube, dpival 
             = dpival, bri = brightness, con = contrast, 
              cmap = cmap, ellipse = ell, plate = scale)


    if noplot is False:
        plt.pause(1.5)
 
    plt.close()

    #}
    #####################


def Comp2Ellip(galhead, galcomps, N, lw=1):
    ''' converts galfit component parameter into an Ellipse object''' 


    ellipses = [] 

    #color value
    values = range(N)
    jet = cm = plt.get_cmap('jet') 
    cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    
    for idx, item in enumerate(galcomps.N):

        if galcomps.Active[idx] == True:
            # correcting coordinates
            xc = galcomps.PosX[idx] - galhead.xmin + 1
            yc = galcomps.PosY[idx] - galhead.ymin + 1


            pa = galcomps.PosAng[idx] + 90

            w = galcomps.Rad[idx]
            h = galcomps.Rad[idx]*galcomps.AxRat[idx]


            colorVal = scalarMap.to_rgba(values[idx])


            ell=Ellipse((xc, yc), width = w, height = h, angle = pa,
                         edgecolor = colorVal,
                         facecolor = 'none',
                         linewidth = lw)

            ellipses.append(ell)


    return ellipses






class ShowCube:

    def __init__(self, cubeimg: str, namepng="cubeout.png", dpival=100, 
                bri = 0, con = 1, cmap='viridis', ellipse=[], plate = 1):
        """
        This routine shows the GALFIT output cube image: galaxy, model and residual    
        """

        
        hdu = fits.open(cubeimg)
        data = (hdu[1].data.copy()).astype(float)
        model = (hdu[2].data.copy()).astype(float)
        residual = (hdu[3].data.copy()).astype(float)
        hdu.close()

        flatmodimg = model.flatten()  
        flatresimg = residual.flatten()  

        flatmodimg.sort()
        flatresimg.sort()

        restot = len(flatresimg)

        restop = round(.9*restot)
        resbot = round(.1*restot)

        modtot = len(flatmodimg)

        modtop = round(.9*modtot)
        modbot = round(.1*modtot)



        modimgpatch = flatmodimg#[modbot:modtop]
        resimgpatch = flatresimg[resbot:restop]

        resmin = np.min(resimgpatch)
        resmax = np.max(resimgpatch)


        modmin = np.min(modimgpatch)
        modmax = np.max(modimgpatch)


        data = data.clip(modmax/1e4,modmax)
        model = model.clip(modmax/1e4)

        modmin = modmax/1e4


        middle = (modmax - modmin)/2


        #brightness auto-adjust according to the contrast value 

        Autobri = middle*(con -1) + modmin*(1-con) 


        #user can re-adjust brightness in addition to Autobri
        newdata = con*(data - middle) + middle + Autobri + bri*(modmax-middle)
        newmodel = con*(model - middle) + middle + Autobri + bri*(modmax-middle)



        mask=data < 0 
        data[mask] = 1 # avoids problems in log
     
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows = 1, ncols = 3)

        im = ax1.imshow(newdata, origin ='lower', interpolation='nearest', norm 
                    = colors.LogNorm(vmin=modmin, vmax=modmax), cmap = cmap)


        ax1.set_title('Data')

        y,x = data.shape

        xt = .02*x
        yt = .02*y

        lxline = round(.1*x)
        lyline = round(.1*y)
        x1 = [xt, xt+lxline]
        y1 = [lyline, lyline]

        arcsec = lxline*plate*U.arcsec 


        if arcsec.value >= 60:
            lxlinearc = arcsec.to("arcmin").value
            s = "{}\'".format(round(lxlinearc))
        else:
            lxlinearc = arcsec.value
            s = "{}\'\'".format(round(lxlinearc))
 
        ax1.plot(x1, y1, color="white", linewidth=3)

        ax1.text(xt+round(lxline/5),lyline+yt,s,color='white',fontsize=14)
    
        #ax1.set_xlabel(r'$\circ$')
        #ax2.set_ylabel('23\"')

        for ell in ellipse:
            ax1.add_patch(ell)


        ax2.imshow(newmodel, origin='lower', interpolation='nearest', norm 
                    = colors.LogNorm(vmin = modmin, vmax = modmax), cmap = cmap)


        ax2.set_title('GALFIT Model')

        ax3.imshow(residual, origin='lower', vmin = resmin, vmax = resmax, cmap = cmap)
        ax3.set_title('Residual')

        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)
        plt.savefig(namepng, dpi = dpival)
    


#end of program
if __name__ == '__main__':
    mainShowCube()
