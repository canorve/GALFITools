#! /usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from astropy.io import fits
from scipy import stats
from matplotlib.patches import Ellipse

def main():



    # uncomment to see the examples run with your own GALFIT files

    #ShowCube("A1656-509-bar.fits",namepng="A1656-509-bar.png")
    #ShowCube("A85gal.fits",namepng="A85gal.png",frac=1)
    #ShowCube("A2256  -576.fits",namepng="A2256-576.png")

    #ell1=Ellipse((399, 385), width=25, height=14,angle=113,
    #                 edgecolor='red',
    #                 facecolor='none',
    #                 linewidth=2)

    #ell2=Ellipse((399, 385), width=55, height=24,angle=113,
    #                 edgecolor='blue',
    #                 facecolor='none',
    #                 linewidth=2)
 
    #ell=[ell1,ell2]


    #ShowCube("A2029cD.fits",namepng="A2029cD.png",frac=0.3,ellipse=ell)




class ShowCube:

    def __init__(self, cubeimg: str,namepng="cubeout.png",dpival=100,frac= 0.2,cmap='viridis',ellipse=[]):
        """
        This routine shows the GALFIT output cube image: galaxy, model and residual    
        """

        
        hdu = fits.open(cubeimg)
        data = (hdu[1].data.copy()).astype(float)
        model = (hdu[2].data.copy()).astype(float)
        residual = (hdu[3].data.copy()).astype(float)
        hdu.close()

        flatmodimg=model.flatten()  
        flatresimg=residual.flatten()  

        flatmodimg.sort()
        flatresimg.sort()

        restot=len(flatresimg)

        restop=round(.9*restot)
        resbot=round(.1*restot)

        modimgpatch=flatmodimg#[modbot:modtop]
        resimgpatch=flatresimg[resbot:restop]

        modmin = np.min(modimgpatch)
        modmax = np.max(modimgpatch)

        if frac  < 1:
            modmin = (1-frac)*modmin 
            modmax = frac*modmax


        if (modmin > modmax):
            modmin, modmax = modmax, modmin


        resmin = np.min(resimgpatch)
        resmax = np.max(resimgpatch)


        mask=data < 0 
        data[mask] = 1 # avoids problems in log
     
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(14, 5), nrows=1, ncols=3)
        fig.subplots_adjust(left=0.04, right=0.98, bottom=0.02, top=0.98)

        #ax1.imshow(data, origin='lower',vmin=modmin, vmax=modmax,cmap=cmap)
        ax1.imshow(data, origin='lower',norm=colors.LogNorm(vmin=modmin, vmax=modmax),cmap=cmap)
        ax1.set_title('Data')


        for ell in ellipse:
            ax1.add_patch(ell)


        ax2.imshow(model, origin='lower',norm=colors.LogNorm(vmin=modmin, vmax=modmax),cmap=cmap)
        #ax2.imshow(model, origin='lower',vmin=modmin, vmax=modmax,cmap=cmap)
        ax2.set_title('GALFIT Model')

        ax3.imshow(residual, origin='lower',vmin=resmin, vmax=resmax,cmap=cmap)
        ax3.set_title('Residual')

        plt.savefig(namepng,dpi=dpival)
     
        plt.show()



#end of program
if __name__ == '__main__':
    main()
