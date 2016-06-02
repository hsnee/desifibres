import numpy as np
import h5py
import matplotlib.pyplot as plt

f = h5py.File('./j3_weights/Aardvark.hdf5','r')
z = f['z'][...]
w = f['weight'][...]

y0, binEdges0 = np.histogram(z,bins=45,normed=True)

bincentres0 = 0.5*(binEdges0[1:]+binEdges0[:-1])

y1,binEdges1 = np.histogram(z,weights=w,bins=45,normed=True)

bincentres1 = 0.5*(binEdges1[1:]+binEdges1[:-1])

plt.figure()
plt.plot(bincentres0, y0, '-', color='green', label='unweighted')
plt.plot(bincentres1, y1, '-', color='red',   label=  'weighted')
plt.xlabel('z')
plt.ylabel('N(z) - normalised')
plt.title('Redshift distributions for complete at 19.7 before and after applying j3 weighting')
plt.legend(fancybox=True)
plt.show()
