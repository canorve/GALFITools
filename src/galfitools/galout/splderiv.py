
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import numpy as np

#check modify

#rad,me=np.genfromtxt("Blackout.txt",unpack=True)
rad,me=np.genfromtxt("Blueout.txt",unpack=True)
#rad,me=np.genfromtxt("Redout.txt",unpack=True)
#rad,me=np.genfromtxt("Greenout.txt",unpack=True)
#rad,me=np.genfromtxt("Magentaout.txt",unpack=True)

lrad= np.log10(rad)
lme= np.log10(me)
 
yspl = UnivariateSpline(lrad,lme,s=0,k=4)

plt.plot(lrad,lme,'ro',label = 'data')
xran = np.linspace(lrad[0],lrad[-1],1000)


plt.plot(xran, yspl(xran))

plt.grid(True)
plt.minorticks_on()
plt.savefig("Data.png")

plt.pause(1.5)

plt.clf()

yspl2d = yspl.derivative(n=2)
yspl1d = yspl.derivative(n=1)

plt.plot(xran,yspl2d(xran))

plt.grid(True)
plt.minorticks_on()
plt.savefig("BreakSpl.png")



plt.pause(1.5)

plt.clf()



dev1 = yspl1d(xran)
dev2 = yspl2d(xran)

kappa = np.abs(dev2)/(1 + dev1**2)**(3/2)

plt.plot(xran,kappa)


plt.grid(True)
plt.minorticks_on()
plt.savefig("KappaSpl.png")



idx = np.where(yspl2d(xran) == min(yspl2d(xran)))

brad = xran[idx[0][0]]

rad = 10**(brad)

print("break radius = ",rad)





