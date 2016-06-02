import os
import time
from astropy.cosmology import FlatLambdaCDM
import h5py
import numpy as np
import matplotlib.pyplot as plt
# from dV import dV
# from dV import Nz

cosmo = FlatLambdaCDM(H0=100, Om0=0.3)

a = h5py.File('/gpfs/data/icc-sum1/Aardvark0.hdf5', 'r+')
za = a["z"][...]

y0,binEdges0=np.histogram(za,bins=60)
z = 0.5*(binEdges0[1:]+binEdges0[:-1])
dz = binEdges0[2] - binEdges0[1]
rc = cosmo.comoving_distance(z)

om0=1
H0=100
c = 299792.458
dV = rc**2 * c * 1.738 * dz / ( H0 * (om0 * (1 + z)**3 + 1 - om0 )**0.5 )
Nz = y0 / dV

yfit = np.poly1d(np.polyfit(z,Nz,15))
nis = yfit(za)
wg = 1 / ( 1 + 3000 * nis)

h5z = a['w']
h5z.attrs['Definition']="weight"

fig = plt.figure()
plt.plot(za,wg,',')
plt.xlabel('z')
plt.ylabel(r'n(z) $h^3$ $Mpc^{-3}$ ')
fig.suptitle('Number density of all galaxies in Aardvark at r = 19.3')
plt.show()