def ang_cor(magnitude, cell_weights=True, j3_weights=False, debug=False):
    import numpy as np
    import matplotlib.pyplot as plt
    ##############
    # magnitude: #
    # 3 --> 19.3 #
    # 5 --> 19.5 #
    # 7 --> 19.7 #
    # 9 --> 19.9 #
    ##############

    # Directories:
    main_folder = '/gpfs/data/icc-sum1/Catalogues/r19'+str(magnitude)+'/Angular/LogBins/'
    if j3_weights == True:
        comp_subfolder = '/j3_weights/'
    else:
        comp_subfolder = '/no_weights/'

    if cell_weights == True:
        targ_subfolder = '/cell_weights/'
    elif j3_weights == True:
        targ_subfolder = '/j3_weights/'
    else:
	targ_subfolder = '/no_weights/'        
    comp_folder = main_folder+'Complete'+comp_subfolder
    targ_folder = main_folder+'Targeted'+targ_subfolder

    comp_cat = comp_folder+'ang.dat'
    targ_cat = targ_folder+'ang.dat'
    
    cordata = [comp_cat,targ_cat] 
    
    # Read in the correlations
    x,y  = [[[0] for _ in range(2)] for _ in range(2)]
    for i in range(len(cordata)):
        cor_array1 = np.loadtxt(cordata[i],dtype='f8')
        cor_array  = np.array(cor_array1)
        x[i] = cor_array[:,0]
        y[i] = cor_array[:,1]
    
    x=x[0]
    y_comp = y[0]
    y_targ = y[1]
    y_ratio = (1.0 + y_comp)/(1 + y_targ)

    if debug==True:
        plt.figure()
        plt.plot(x,y_ratio,'+')
        plt.xscale('log')
        plt.xlabel(r'$\theta$')
        plt.ylabel(r'$(1 + \omega(\theta)|all)/(1 + \omega(\theta)|targ)$')
        plt.show()
    else:
        pass
    return y_ratio
