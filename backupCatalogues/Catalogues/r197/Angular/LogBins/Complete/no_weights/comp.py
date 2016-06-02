import numpy as np
import matplotlib.pyplot as plt

folder = '/gpfs/data/icc-sum1/Catalogues/r197/Angular/LogBins/Complete/'

path1 = 'ang.dat'
path2 = 'true_no_weights/ang.dat'
path3 = 'sample/ang.dat'
path4 = 'redshift_slices/point1_point2/ang.dat'
path5 = 'redshift_slices/point2_point4/ang.dat'
path6 = 'redshift_slices/point6_1point4/ang.dat'
path7 = 'redshift_slices/point4_point6/ang.dat'
path8 = 'redshift_slices/point0_poin1/ang.dat'

path = [path1, path2, path3, path4, path5, path6, path7, path8]

x, y = [[[0] for _ in range(8)] for _ in range(2)]

for i in range(len(path)):
    cor_array1 = np.loadtxt(path[i],dtype='f8')
    cor_array  = np.array(cor_array1)
    x[i] = cor_array[:,0]
    y[i] = cor_array[:,1]
    
plt.figure(1)
plt.plot(x[0],y[0],'gv', markersize=3, label='Unity Weights')
plt.plot(x[1],y[1],'b2', markersize=7, label='No Weights')
plt.plot(x[2],y[2],'r+', markersize=6, label='Random Sample')
plt.xscale('log')
plt.yscale('log')
plt.title('Angular correlation function for unweighted complete at r = 19.7')
plt.legend(fancybox=True)
plt.show()

plt.figure(2)
plt.plot(x[2],y[2]/y[1],'r+', label='Sample/No weights')
plt.plot(x[2],y[0]/y[1],'gv', label ='Unity weights/No weights')
plt.xscale('log')
plt.yscale('log')
plt.legend(fancybox=True)
plt.show()

plt.figure(3)
plt.plot(x[7],y[7],'kd', markersize=5, label='0.0<z<0.1')
plt.plot(x[3],y[3],'gv', markersize=3, label='0.1<z<0.2')
plt.plot(x[4],y[4],'b2', markersize=7, label='0.2<z<0.4')
plt.plot(x[6],y[6],'yo', markersize=3, label='0.4<z<0.6')
plt.plot(x[5],y[5],'r+', markersize=6, label='0.6<z<1.4')
#plt.plot(x[7],y[7],'kd', markersize=5, label='0.0<z<0.1')
plt.title('Angular correlation function for different redshift slices in unweighted complete at r = 19.7')
plt.xscale('log')
plt.yscale('log')
plt.legend(fancybox=True)
plt.show()

