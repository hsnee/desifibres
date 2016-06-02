
def grid_output(numra, numdec,magnitude=7,interp='none',bins=20,comparable=False,really_comparable=False):
    
    ######################################
    # Function inputs:
    #
    # This function only provides an easy way to run weights from grid.py and plots the results
    # numra is the number of cells in ra
    # numdec is the number of cells in sin(dec)
    # magnitude assumes that the folders are saved under:
    # /gpfs/data/icc-sum1/Angular/MAGNITUDE/AardvarkNUMb.hdft
    # where MAGNITUDE is 193, 195, 197 or 199
    # and NUM is between 0 and 10 (see lines 24 to 35 below)
    # change interp if you would like python to plot them with an interpolation
    # change bins if you would like a different number of bins in the histogram (figure4)
    # comparable sets and minimum and maximum values for the colors so that the weights are comparable at different magnitude limits.
    # really_comparable plots weight distributions for all magnitude limits, this will take a very long time to run. 
    #####################################

    import numpy as np
    import matplotlib.pyplot as plt
    from grid import cell_weights as weights

    if magnitude==3:
        complete_catalogue='0'
        targeted_catalogue='1'
    elif magnitude==5:
        complete_catalogue='3'
        targeted_catalogue='4'
    elif magnitude==7:
        complete_catalogue='6'
        targeted_catalogue='7'
    else:
        complete_catalogue='9'
        targeted_catalogue='10'
    ws = weights('/gpfs/data/icc-sum1/Angular/19'+str(magnitude)+'/Aardvark'+str(complete_catalogue)+'b.hdf5','/gpfs/data/icc-sum1/Angular/19'+str(magnitude)+'/Aardvark'+str(targeted_catalogue)+'b.hdf5',numra,numdec)
    w = ws[0]
    xedges = ws[1]
    yedges = ws[2]
    comp_array_normed = ws[3]
    targ_array_normed = ws[4]
    
    # Must know what the max of comp_array is to force it as colorbar max.
    vmx = np.amax(map(max,comp_array_normed))
    vmn = np.amin(map(min,targ_array_normed))

    # Must define limits for imshow.
    extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
    
    # Firstly, a plot of the weights, i.e. ratio of number of galaxies
    # in the complete catalogue over the number of galaxies that 
    # were targeted by the fibres.
    
    fig = plt.figure(1)
    imgplot = plt.imshow(w,extent=extent,aspect='auto',interpolation=interp)
    plt.xlabel('ra')
    if comparable==True:
        plt.clim(1.02,1.2)
    else:
        pass
    plt.ylabel(r'sin ($\delta$)')
    plt.title('Ratio of number of galaxies in complete/targeted catalogues ('+str(numra)+','+str(numdec)+'), magnitude=19.'+str(magnitude)+'')
    plt.colorbar()
    plt.show()
    
    # Secondly, a plot of the number density per cell per square degrees for the complete survey
    
    fig = plt.figure(2)
    imgplot = plt.imshow(comp_array_normed,extent=extent,aspect='auto',interpolation=interp)
    plt.xlabel('ra')
    plt.ylabel(r'sin ($\delta$)')
    cbar = plt.colorbar() 
    cbar.set_label(r'Number of Galaxies per $degrees^2$')
    plt.clim(vmn,vmx)
    plt.title("Number of galaxies per cell per square degrees (complete) ["+str(numra)+", "+str(numdec)+"],magnitude=19."+str(magnitude)+"")
    plt.show()
    
    # Thirdly, a similar plot for the targeted survey
    
    fig = plt.figure(3)
    imgplot = plt.imshow(targ_array_normed,extent=extent,aspect='auto',interpolation=interp)
    plt.xlabel('ra')
    plt.ylabel(r'sin ($\delta$)')
    cbar = plt.colorbar() 
    cbar.set_label(r'Number of Galaxies per $degrees^2$')
    plt.clim(vmn,vmx)
    plt.title('Number of galaxies per cell per square degrees (targeted) ['+str(numra)+', '+str(numdec)+'],magnitude=19.'+str(magnitude)+'')
    plt.show()
    
    # Lastly, weight distributions
    
    weights_vector = np.ravel(w)
    plt.figure(4)
    y0,binEdges0=np.histogram(weights_vector,bins=bins, normed=True)
    bincenters0 = 0.5*(binEdges0[1:]+binEdges0[:-1])
    plt.plot(bincenters0,y0,linestyle="--")
    if comparable==True:
        plt.xlim([1, 1.2])
    else:
        pass
    plt.xlabel('Weights')
    plt.ylabel('Counts')
    plt.title('Weight distribution for a '+str(numra)+' by '+str(numdec)+' cells grid, magnitude=19.'+str(magnitude)+'')
    plt.show()

    if really_comparable==True:
        ws1 = weights('/gpfs/data/icc-sum1/Angular/193/Aardvark0b.hdf5','/gpfs/data/icc-sum1/Angular/193/Aardvark1b.hdf5',numra,numdec)
        ws2 = weights('/gpfs/data/icc-sum1/Angular/195/Aardvark3b.hdf5','/gpfs/data/icc-sum1/Angular/195/Aardvark4b.hdf5',numra,numdec)
        ws3 = weights('/gpfs/data/icc-sum1/Angular/197/Aardvark6b.hdf5','/gpfs/data/icc-sum1/Angular/197/Aardvark7b.hdf5',numra,numdec)
        ws4 = weights('/gpfs/data/icc-sum1/Angular/199/Aardvark9b.hdf5','/gpfs/data/icc-sum1/Angular/199/Aardvark10b.hdf5',numra,numdec)
        w1 = ws1[0]
        w2 = ws2[0]
        w3 = ws3[0]
        w4 = ws4[0]
        weights_vector1 = np.ravel(w1)
        weights_vector2 = np.ravel(w2)
        weights_vector3 = np.ravel(w3)
        weights_vector4 = np.ravel(w4)

        plt.figure(5)
        y0,binEdges0=np.histogram(weights_vector1,bins=bins, normed=True)
        y1,binEdges1=np.histogram(weights_vector2,bins=bins, normed=True)
        y2,binEdges2=np.histogram(weights_vector3,bins=bins, normed=True)
        y3,binEdges3=np.histogram(weights_vector4,bins=bins, normed=True)

        bincenters0 = 0.5*(binEdges0[1:]+binEdges0[:-1])
        plt.plot(bincenters0,y0,linestyle="-",color='blue', label="r = 19.3")
        bincenters1 = 0.5*(binEdges1[1:]+binEdges1[:-1])
        plt.plot(bincenters1,y1,linestyle="-", color='red', label="r = 19.5")
        bincenters2 = 0.5*(binEdges2[1:]+binEdges2[:-1])
        plt.plot(bincenters2,y2,linestyle="-", color='yellow', label="r = 19.7")
        bincenters3 = 0.5*(binEdges3[1:]+binEdges3[:-1])
        plt.plot(bincenters3,y3,linestyle="-", color='green', label="r = 19.9")
        plt.xlim([1, 1.2])
        plt.xlabel('Weights')
        plt.ylabel('Counts')
        plt.legend()
        plt.title('Weight distribution for a '+str(numra)+' by '+str(numdec)+' cells grid')
        plt.show()

    else: 
        pass
