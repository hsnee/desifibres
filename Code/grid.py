##########################################################################################
## This script is to be used to "lay a grid" in cartesian ra, sin(dec) to 2 catalogues 	##
## and take the ratio of the number of particles in each grid between the 2 catalogues  ##
## to apply this ratio as an extra weight to make large scale correction for redshift   ##
## incompleteness. 									##
##########################################################################################
from __future__ import division
def weights(complete,targeted,numra,numdec):
	import h5py
	import numpy as np
	import time
	import matplotlib.pyplot as plt
	###################################
	# Things to input (if not using as function)
	# import the two catalogues
	# complete =
	# targeted =
	#
	# number of cells in ra
	# numra = 
	#
	# number of cells in sin(dec)
	# numdec = 
        ####################################
        ####################################
	print 'Remember to have used future division'
	print 'Reading Complete survey as: ',complete
	time.sleep(1)
	print 'Reading Targeted survey as: ',targeted
	time.sleep(1)
	print 'Number of cells in ra = ', numra
	time.sleep(1)
	print 'Number of cells in sin (dec) = ',numdec
	
	comp = h5py.File(complete, "r")
	targ = h5py.File(targeted, "r")
	
	numprod = numra*numdec
	
	rac  = comp["ra"][...]
	decc = comp["dec"][...]
	decc = np.sin(decc*np.pi/180.0)
	rat  = targ["ra"][...]
	dect = targ["dec"][...]
	dect = np.sin(dect*np.pi/180.0)
	
	
	# get the range of ra and dec, it's enough to do this for either catalogues
	rangera  = np.amax(rac)  - np.amin(rac)
	rangedec = np.amax(decc) - np.amin(decc)
	
	
	# make up an array of cell boundaries using the information above 
	x = [[0] for _ in range(numra)]
	y = [[0] for _ in range(numdec)]
	weights    = [[2. for _ in range(numra)] for _ in range(numdec)]
	comp_array = [[2. for _ in range(numra)] for _ in range(numdec)]
	targ_array = [[2. for _ in range(numra)] for _ in range(numdec)]
	
	for i in range(len(x)):
		x[i] = np.amin(rac) + rangera*(i+1)/numra
		
       	for i in range(len(y)):
       		y[i] = np.amax(decc) - rangedec*(i+1)/numdec

	# Let's calculate the (average) size of a cell in square degrees
	cell_size = (70/numdec) * (130/numra)
	
       	# Now, x and y are the upper boundaries of the cells (i.e. x[0] is upper boundary of cell 0.
	
               
	if numra==numdec:
		comp_array, xedges0, yedges0 = np.histogram2d(rac, decc, numra )
		targ_array, xedges1, yedges1 = np.histogram2d(rat, dect, numra )
		weights = comp_array/targ_array
		comp_array_normed = comp_array/cell_size
		targ_array_normed = targ_array/cell_size
		return (weights,xedges0,yedges0,comp_array_normed,targ_array_normed)
	else:
		for i in range(len(x)):
			for j in range(len(y)):
				if i+j==0:
					comp_array[j][i] = len(rac[(rac<x[i])&(decc>y[j])])
					targ_array[j][i] = len(rat[(rat<x[i])&(dect>y[j])])
				elif i==0:
					comp_array[j][i] = len(rac[(rac<x[i])&(decc>y[j])&(decc<=y[j-1])])
					targ_array[j][i] = len(rat[(rat<x[i])&(dect>y[j])&(dect<=y[j-1])])
				elif j==0:
					comp_array[j][i] = len(rac[(rac<x[i])&(rac>=x[i-1])&(decc>y[j])])
					targ_array[j][i] = len(rat[(rat<x[i])&(rat>=x[i-1])&(dect>y[j])])
				else:
					comp_array[j][i] = len(rac[(rac<x[i])&(decc>y[j])&(rac>=x[i-1])&(decc<=y[j-1])])
					targ_array[j][i] = len(rat[(rat<x[i])&(dect>y[j])&(rat>=x[i-1])&(dect<=y[j-1])])
		for i in range(len(x)):
			for j in range(len(y)):
				weights[j][i] = comp_array[j][i]/targ_array[j][i]
		xxmin = [np.amin(rac)]
		yymax = [np.amax(decc)]
		xedges0 = np.concatenate([xxmin,x])
		yedges0 = np.concatenate([yymax,y])
		yedges0 = yedges0[::-1]
		comp_array = np.array(comp_array)
		targ_array = np.array(targ_array)
		comp_array_normed = comp_array / cell_size
		targ_array_normed = targ_array / cell_size
		# Now let's return not only the weights but also the galaxy numbers per cell per square degrees;
		# edges are the still the same. 
		return (weights,xedges0,yedges0,comp_array_normed,targ_array_normed)

