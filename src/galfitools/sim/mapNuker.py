#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from os import path

import mgefit
from mgefit.mge_fit_1d import mge_fit_1d
from scipy.special import gammaincinv


from scipy.interpolate import UnivariateSpline

def main():



    rb = 10 

    print("computing mean of Rbreak/Rgamma")

    mean, dev, mi,ma = Rgamean(rb)

    print("mean = {:.2f} ; std = {:.2f} ".format(mean,dev))
    print("min value = {:.2f} ; max = {:.2f}".format(mi,ma))


def Rgamma(Rb,alpha,beta,gamma):


    rgam = Rb*(((0.5-gamma)/(beta-0.5))**(1/alpha)) #formula

    if np.isscalar(rgam):
        Rg  = rgam
    else:
        mask = np.isfinite(rgam)
        Rg  = rgam[mask]

    return Rg 



def Rgamean(Rb):

    alpha = np.arange(2,4.1,0.01) #range of alpha to search
    beta = np.arange(0.5,3.1,0.01) #range of beta to search
    gamma = np.arange(0,0.5,0.01) #range of gamma to search

    alpha = np.arange(1,4.1,0.01) #range of alpha to search
    beta = np.arange(0.5,2.1,0.01) #range of beta to search
    gamma = np.arange(0,0.5,0.01) #range of gamma to search



    albega = cartesian((alpha,beta,gamma))
    
    alpha= albega[:,0]
    beta = albega[:,1]
    gamma = albega[:,2]

    Rg = Rgamma(Rb,alpha,beta,gamma)

    #mask = np.isfinite(Rg)
    #output[~np.isfinite(output)] = 0

    kons=Rb/Rg

    kons.sort()

    tot=len(kons)
    
    print("removing top 80% and bottom 20%")

    top=round(.8*tot)
    bot=round(.2*tot)
    #top=np.int(top)
    #bot=np.int(bot)

    kon=kons[bot:top]

    n = len(kon)

    print("the total number of parameters combinations is",n)

    if n%2 == 0:
        idmed1 = int(n/2 - 1)
        idmed2 = int(n/2) 
        median = (kon[idmed1] + kon[idmed2])/2 
    else:
        idmed = int((n+1)/2 - 1)
        median = kon[idmed]

    print("the median is {:.2f} ".format(median))

    m = np.mean(kon)
    s = np.std(kon) 
    mi = np.min(kon)
    ma = np.max(kon)

    plt.close()
    plt.clf()

    counts, bins = np.histogram(kon)
    plt.stairs(counts, bins)
    plt.hist(bins[:-1], bins, weights=counts)
    #counts,bins,patches=plt.hist(kon)

    idmod = np.where(counts==max(counts))[0][0] 

    print("the mode is {:.2f} ".format(bins[idmod]))

    #plt.hist(kon)


    plt.grid(True)
    plt.minorticks_on()
    plt.savefig("histn.png")
 


    return m,s,mi,ma







def cartesian(arrays, out=None):
    """
    Generate a Cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the Cartesian product of.
    out : ndarray
        Array to place the Cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing Cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    #m = n / arrays[0].size
    m = int(n / arrays[0].size)
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
        #for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m, 1:] = out[0:m, 1:]
    return out



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
    main()
