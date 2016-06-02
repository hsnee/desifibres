
###################### Arguments ###################
# magnitude = 3,5,7 or 9.                          #
# cortype = 'Angular' or 'Monopole'                #
# binning = 'LogBins' or 'LinBins'                 #
# complete, targ_raw and targ_corrected:           #
# 'no_weights', 'j3_weights', 'cell_weights'       #
# 'woftheta_weights', 'j3_cell_weights',           #
# 'woftheta_cell_weights', 'all_weights'           #
# correlation: 'Complete/any_of_the_above' or      #
# 'Targeted/any_of_the_above'                      #
####################################################


def juxtapose(magnitude,cortype,binning,correlation1,correlation2):
    import numpy as np
    import matplotlib.pyplot as plt
    if cortype == 'Angular':
        cor = 'ang'
    else:
        cor = 'mono'
    path1 = './r19'+str(magnitude)+'/'+str(cortype)+'/'+str(binning)+'/'+str(correlation1)+'/'+cor+'.dat'
    path2 = './r19'+str(magnitude)+'/'+str(cortype)+'/'+str(binning)+'/'+str(correlation2)+'/'+cor+'.dat'
    cor_array1 = np.loadtxt(path1,dtype='f8')
    cor_array  = np.array(cor_array1)
    x = cor_array[:,0]
    y = cor_array[:,1]
    
    cor_array2a = np.loadtxt(path2,dtype='f8')
    cor_array2 = np.array(cor_array2a)
    x1 = cor_array2[:,0]
    y1 = cor_array2[:,1]
    plt.figure()
    plt.plot(x,y,'gv',label=str(correlation1))
    plt.plot(x1,y1,'r+',label=str(correlation2))
    plt.xscale('log')
    plt.yscale('log')
    if cortype=='Angular':
        plt.xlabel(r'$\theta^o$')
        plt.ylabel(r'$\omega(\theta)$')
    else:
        plt.xlabel(r'$r$')
        plt.ylabel(r'$\xi(r)$')
    plt.title(cortype+'correlation function at 19.'+str(magnitude))
    plt.legend(fancybox=True)
    plt.show()

def triple(magnitude,cortype,binning,complete,targ_raw,targ_corrected):
    import numpy as np
    import matplotlib.pyplot as plt
    if cortype == 'Angular':
        cor = 'ang'
    else:
        cor = 'mono'
    path1 = './r19'+str(magnitude)+'/'+str(cortype)+'/'+str(binning)+'/'+str(complete)+'/'+cor+'.dat'
    path2 = './r19'+str(magnitude)+'/'+str(cortype)+'/'+str(binning)+'/'+str(targ_raw)+'/'+cor+'.dat'
    path3 = './r19'+str(magnitude)+'/'+str(cortype)+'/'+str(binning)+'/'+str(targ_corrected)+'/'+cor+'.dat'
    path = [path1, path2, path3]
    x, y = [[[0] for _ in range(3)] for _ in range(2)]
    for i in range(len(path)):
        cor_array1 = np.loadtxt(path[i],dtype='f8')
        cor_array  = np.array(cor_array1)
        x[i] = cor_array[:,0]
        y[i] = cor_array[:,1]
        
    plt.figure(1)
    plt.plot(x[0],y[0],'gv', markersize=3, label=str(complete))
    plt.plot(x[1],y[1],'b2', markersize=7, label=str(targ_raw))
    plt.plot(x[2],y[2],'r+', markersize=6, label=str(targ_corrected))
    plt.xscale('log')
    if cortype=='Angular':
        plt.xlabel(r'$\theta^o$')
        plt.ylabel(r'$\omega(\theta)$')
    else:
        plt.xlabel(r'$r$ (Mpc/h)')
        plt.ylabel(r'$\xi(r)$')
    plt.yscale('log')
    plt.legend(fancybox=True)
    if targ_corrected=='Targeted/cell_weights':
        plt.title(cortype+' correlation function at 19.'+str(magnitude)+' for 44 by 24 cells')
    else:
        plt.title(cortype+' correlation function at 19.'+str(magnitude))
    if binning=='LogBins':
        xline=np.logspace(-3,1,1000)
    else:
        xline=np.logspace(0,3,1000)
    yline=np.ones(1000)
    plt.figure(2)
    plt.plot(x[0],y[2]/y[0],'r+',label='corrected/complete')
    plt.plot(x[0],y[1]/y[0],'b2',label='raw/complete')
    plt.plot(xline,yline,'b-')
    if cortype== 'Monopole':
        plt.ylabel(r'$\xi(r)|targ/\xi(r)/comp$')
        plt.xlabel('r (Mpc/h)')
    else:
        plt.ylabel(r'$\omega(\theta)|targ/\omega(\theta)|comp$')
        plt.xlabel(r'$\theta^o$')
    plt.xscale('log')
    plt.title(cortype+' correlation function (r = 19.'+str(magnitude)+')')
    plt.yscale('log')
    plt.ylim(0.1,10)
    plt.legend()
    plt.show()
    
def compare(cortype,binning):
	import numpy as np
	import matplotlib.pyplot as plt
	if cortype == 'Angular':
		cor = 'ang'
	else:
		cor = 'mono'
	path0a = './r193/'+str(cortype)+'/'+str(binning)+'/Complete/no_weights/'+cor+'.dat'
	path1a = './r195/'+str(cortype)+'/'+str(binning)+'/Complete/no_weights/'+cor+'.dat'
	path2a = './r197/'+str(cortype)+'/'+str(binning)+'/Complete/no_weights/'+cor+'.dat'
	path3a = './r199/'+str(cortype)+'/'+str(binning)+'/Complete/no_weights/'+cor+'.dat'
	path0b = './r193/'+str(cortype)+'/'+str(binning)+'/Targeted/no_weights/'+cor+'.dat'
	path1b = './r195/'+str(cortype)+'/'+str(binning)+'/Targeted/no_weights/'+cor+'.dat'
	path2b = './r197/'+str(cortype)+'/'+str(binning)+'/Targeted/no_weights/'+cor+'.dat'
	path3b = './r199/'+str(cortype)+'/'+str(binning)+'/Targeted/no_weights/'+cor+'.dat'
	path = [path0a, path0b, path1a, path1b, path2a, path2b, path3a, path3b]
	x, y = [[[0] for _ in range(8)] for _ in range(2)]
	for i in range(len(path)):
		cor_array1 = np.loadtxt(path[i],dtype='f8')
		cor_array  = np.array(cor_array1)
		x[i] = cor_array[:,0]
		y[i] = cor_array[:,1]
	if binning=='LogBins':
		xline=np.logspace(-3,1,1000)
	else:
		xline=np.logspace(0,3,1000)
	yline=np.ones(1000)
	fig, ax = plt.subplots()
	ax.plot(x[0],y[1]/y[0],'+', markersize=5, color='b', label='r = 19.3')
	ax.plot(x[2],y[3]/y[2],'2', markersize=7, color='m', label='r = 19.5')
	ax.plot(x[4],y[5]/y[4],'o', markersize=2, color='g', label='r = 19.7')
	ax.plot(x[6],y[7]/y[6],'^', markersize=3, color='r', label='r = 19.9')
	ax.plot(xline,yline,'k-')
	ax.set_xlabel('r (Mpc/h)')
	ax.set_ylabel(r'$\xi(r)|targ/\xi(r)/comp$')
	ax.set_xscale('log')
	ax.axvspan(100, 110, ymin=-1, ymax=1.4, alpha=0.3, facecolor='red',edgecolor='white')
	ax.set_title('Ratios of targeted/complete, unweighted, uncorrected monopole correlation for different magnitudes')
	ax.set_ylim([0.8,1.4])
	ax.legend(fancybox=True)
	plt.show()
