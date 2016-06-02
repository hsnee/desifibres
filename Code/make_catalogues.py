##########################################################################################
## This script is created to make appropriate ra, dec and r selections on the Aardvark	##
## mock survey and run a batch script (which in turn, at least originally intended, to 	##
## submit the "CUTE" code to get the pair correlation function of Aard.			##
##########################################################################################

def make_cat(magnitude, weights=0, debug=False, skip_randoms=False, superdebug=False, same_weights=False, numra=22, numdec=12):
    import time
    import h5py
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import legend
    from astropy.cosmology import FlatLambdaCDM
    from astropy import units as u
    from grid import cell_weights as weigh
    ##############################################
    # Weights: 0 ---> no weight                  #
    #          1 ---> unit weights               #
    #          2 ---> redshift-dependant weights #
    #          3 ---> cell-corrected weights     #
    #          4 ---> 2 and 3                    #
    ##############################################
    if weights==0:
        print 'Not using any weights'   
    elif weights==1:
        print 'Using unity weights'
    elif weights==2:
        print 'Using redshift-dependant weights'
    elif weights==3:
        print 'Using cell weights'
    elif weights==4:
        print 'Using cell and redshift-dependant weights'
    else:
        pass
    if skip_randoms==True:
        print 'Skipping randoms'
    else:
        pass
    
    # Get Cosmology
    cosmo = FlatLambdaCDM(H0=100, Om0=0.3)
    
    # read in the original mock survey and import all the data
    Aardvark  = "/gpfs/data/icc-sum1/aardvark_mockcat.hdf5"
    print 'Reading ',Aardvark
    a = h5py.File(Aardvark, "r")
    za=a["z"][...]
    ra=a["ra"][...]
    dec=a["dec"][...]
    which_magnitude = np.array([0,0])
    if magnitude==3:
        which_magnitude[0] = 0
        which_magnitude[1] = 1
        ffolder = "/Angular/193"
        zx3 = a["istat19.3"][...]
    elif magnitude==5:
        zx3 = a["istat19.5"][...]
        which_magnitude[0] = 3
        which_magnitude[1] = 4
        ffolder = "/Angular/195"
    elif magnitude==7:
        zx3 = a["istat19.7"][...]
        which_magnitude[0] = 6
        which_magnitude[1] = 7
        ffolder = "/Angular/197"
    else:
        zx3 = a["istat19.9"][...]
        which_magnitude[0] = 9
        which_magnitude[1] = 10
        ffolder = "/Angular/199"
    if weights==1: # convention
        which_targeting = 'b'
    elif weights==0:
        which_targeting = 'b'
    else:
        which_targeting = 'a'
   
    # Names of catalogues:
    # Data catalogue is above
    
    oldrname = '/gpfs/data/icc-sum1/rancat.hdf5' # original randoms catalogue
    newtname = "./Aardvark"+str(which_magnitude[1])+which_targeting+".hdf5" # targeted catalogue
    newfname = "./Aardvark"+str(which_magnitude[0])+which_targeting+".hdf5" # complete catalogue
    newrname = "./randomcat"+str(which_magnitude[0])+which_targeting+".hdf5" # randoms corresponding to complete
    newrname_targ='./randomcat'+str(which_magnitude[1])+which_targeting+'.hdf5' # randoms corresponding to targeted
    if superdebug==True:
        print 'D = ', Aardvark
        print 'R = ', oldrname
        print 'D_comp = ', newfname
        print 'R_comp = ', newrname
        print 'D_targ = ', newtname
        print 'R_targ = ', newrname_targ
    else:
        pass
    # cut the edges off ra, dec and z and z3
    ra1=ra[(ra<250)&(ra>120)&(dec<-10)&(dec>-70)]
    dec1=dec[(ra<250)&(ra>120)&(dec<-10)&(dec>-70)]
    z3 = zx3[(ra<250)&(ra>120)&(dec<-10)&(dec>-70)]
    zx = za[(ra<250)&(ra>120)&(dec<-10)&(dec>-70)]
    ra_complete = ra1[z3>0]
    dec_complete = dec1[z3>0]
    redshift_complete = zx[z3>0]
    ra_targeted = ra1[z3==1]
    dec_targeted = dec1[z3==1]
    redshift_targeted = zx[z3==1]
    a.close()
    
    # Write two new Aardvark file, first up: Complete. 
    print 'Writing ',newfname
    newf = h5py.File(newfname, "w")
    h5ra = newf.create_dataset('ra',data=ra_complete)
    h5ra.attrs["Definition"] = "Right Ascension in degrees [0,360]"
    h5dec = newf.create_dataset('dec',data=dec_complete)
    h5dec.attrs['Definition'] = "Declination in degrees [-90,90]" 
    h5z = newf.create_dataset('z',data=redshift_complete)
    h5z.attrs['Definition']="Redshift"
    zx2 = newf['z'][...]
    
    # Determining Weights
    y0,binEdges0=np.histogram(zx2,bins=60)
    z = 0.5*(binEdges0[1:]+binEdges0[:-1])
    dz = binEdges0[2] - binEdges0[1]
    rc = cosmo.comoving_distance(z)
    rc = rc / u.Mpc
    om0=1
    H0=100
    c = 299792.458
    dV = rc**2 * c * 1.738 * dz / ( H0 * (om0 * (1.0 + z)**3 + 1.0 - om0 )**0.5 )
    Nz = y0 / dV
    Nz[Nz==0]=3*10**(-9)
    yfit = np.poly1d(np.polyfit(z,np.log10(Nz),9))
    nis = np.power(10.0,yfit(zx2))
    if (weights==1)|(weights==3):
        print 'Writing unity weights'
        wg = np.ones(len(zx2))
        h5z = newf.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"
        
    elif (weights==2)|(weights==4):
        print 'Writing redshift-dependant weights'
        wg = 1.0 / ( 1.0 + 3000.0 * nis)
        h5z = newf.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"        
    else:
        print 'Not writing any weights'
    if debug==True:
        print 'Weight distribution for complete_data'
        xp = np.linspace(0,1.6,100)
        idxran = np.random.random_integers(0,high=len(wg)-1,size=2000000)
        xran = redshift_complete[idxran]
        yran = wg[idxran]
        fig = plt.figure()
        plt.plot(xran,yran,'+')
        plt.xlabel("Redshift")
        plt.ylabel("Weight per galaxy")
        fig.suptitle('Galaxy weights as a function of redshift for all galaxies in Aardvark at r = 19.'+str(magnitude))
        plt.show()
    else:
        pass
    newf.close() 
    if weights>0:
        del wg
    else:
        pass
    
    ####################################################
    # Now add the weights to the corresponding randoms #
    ####################################################
    if skip_randoms==True:
        pass
    else:
        if weights==0:
            pass
        else:
            print 'Reading ',oldrname
            ran_cat = h5py.File(oldrname, "r")
            ra_ran = ran_cat["ra"][...]
            dec_ran = ran_cat["dec"][...]
            which_magnitude = np.array([0,0])
            if magnitude==3:
                zx3_ran = ran_cat["z19.3"][...]
            elif magnitude==5:
                zx3_ran = ran_cat["z19.5"][...]
            elif magnitude==7:
                zx3_ran = ran_cat["z19.7"][...]
            else:
                zx3_ran = ran_cat["z19.9"][...]

            # Make cutoffs
            ra1_ran = ra_ran[(ra_ran<250)&(ra_ran>120)&(dec_ran<-10)&(dec_ran>-70)]
            dec1_ran=dec_ran[(ra_ran<250)&(ra_ran>120)&(dec_ran<-10)&(dec_ran>-70)]
            z3_ran = zx3_ran[(ra_ran<250)&(ra_ran>120)&(dec_ran<-10)&(dec_ran>-70)]
            ran_cat.close()
            
            # Write new random file
            print 'writing ',newrname
            newr = h5py.File(newrname, "w")
            h5ra_ran = newr.create_dataset('ra',data=ra1_ran)
            h5ra_ran.attrs["Definition"] = "Right Ascension in degrees [0,360]"
            h5dec_ran = newr.create_dataset('dec',data=dec1_ran)
            h5dec_ran.attrs['Definition'] = "Declination in degrees [-90,90]" 
            h5z_ran = newr.create_dataset('z',data=z3_ran)
            h5z_ran.attrs['Definition']="Redshift"
            zz = newr['z']
            newr.close()
            # Now since the randoms for the complete and targeted are the same, let's write it here now, reopen it later on
            print 'writing ',newrname_targ
            newr_targ = h5py.File(newrname_targ,'w')
            h5ra_ran_targ = newr_targ.create_dataset('ra',data=ra1_ran)
            h5ra_ran_targ.attrs['Definition'] = 'Right Ascension in degrees [0,360]'
            h5dec_ran_targ = newr_targ.create_dataset('dec',data=dec1_ran)
            h5dec_ran_targ.attrs['Definition'] = 'Declination in degrees [-90,90]'
            h5z_ran_targ = newr_targ.create_dataset('z',data=z3_ran)
            h5z_ran_targ.attrs['Definition'] = 'Redshift'
            newr_targ.close()
            
            # Read complete in again to apply weights
            newr = h5py.File(newrname,'r+')
            zz = newr['z'][...]
            nis2 = np.power(10.0,yfit(zz))
            if (weights==1)|(weights==3):
                w2 = np.ones(len(zz))
            elif (weights==2)|(weights==4):
                w2 = 1.0 / ( 1.0 + 3000.0 * nis2 )
            else:
                pass
            time.sleep(2)
            h5w2 = newr.create_dataset('weight',data=w2)
            h5w2.attrs['Definition']="weight"
            if debug==True:
                print 'Weight distribution for randoms_complete'
                xp = np.linspace(0,1.6,100)
                idxran = np.random.random_integers(0,high=len(w2)-1,size=2000000)
                xran = zz[idxran]
                yran = w2[idxran]
                fig = plt.figure()
                plt.plot(xran,yran,'+')
                plt.xlabel("Redshift")
                plt.ylabel("Weight per galaxy")
                fig.suptitle('Galaxy weights as a function of redshift for all randoms in Aardvark at r = 19.'+str(magnitude))
                plt.show()
            else:
                pass
            newr.close()
            if weights>0:
                del w2,y0,binEdges0,zx2,z,rc,dV,Nz
            else:
                pass

    ##########################################################
    ##########################################################
    ##            Now write targeted files                  ##
    ##########################################################
    ##########################################################

    print 'writing', newtname
    newt = h5py.File(newtname, "w")
    h5ra = newt.create_dataset('ra',data=ra_targeted)
    h5ra.attrs["Definition"] = "Right Ascension in degrees [0,360]"
    h5dec = newt.create_dataset('dec',data=dec_targeted)
    h5dec.attrs['Definition'] = "Declination in degrees [-90,90]" 
    h5z = newt.create_dataset('z',data=redshift_targeted)
    h5z.attrs['Definition']="Redshift"
    zx2 = newt['z'][...]
    
    # Determining Weights
    y0,binEdges0=np.histogram(zx2,bins=60)
    z = 0.5*(binEdges0[1:]+binEdges0[:-1])
    dz = binEdges0[2] - binEdges0[1]
    rc = cosmo.comoving_distance(z)
    rc = rc / u.Mpc
    om0=1
    H0=100
    c = 299792.458
    dV = rc**2 * c * 1.738 * dz / ( H0 * (om0 * (1.0 + z)**3 + 1.0 - om0 )**0.5 )
    Nz = y0 / dV
    Nz[Nz==0]=3*10**(-9)
    yfit2 = np.poly1d(np.polyfit(z,np.log10(Nz),9))
    
    if same_weights==True:
        nis = np.power(10.0,yfit(zx2))
    else:
        nis = np.power(10.0,yfit2(zx2))
    if weights==1:
        print 'Writing unity weights'
        wg = np.ones(len(zx2))
        h5z = newt.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"
    elif weights==2:
        print 'Writing redshift-dependant weights'
        wg = 1.0 / ( 1.0 + 3000.0 * nis)
        h5z = newt.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"
    elif weights==3:
        W = weigh(newfname,newtname,numra,numdec)
        ws = W[0]
	wg = np.zeros(len(ra_targeted))
	for k in range(len(ra_targeted)):
            indx = (ra_targeted[k] - 120)/(130.0/numra)
            indy = (np.sin(dec_targeted[k]*3.1415/180) + 0.17364817766693033)/-(0.76604444311/numdec)
	    indx = int(indx)
	    indy = int(indy)
	    wg[k] = ws[indy][indx]
        h5z = newt.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"
    elif weights==4:
        W = weigh(newfname,newtname,numra,numdec)
        ws = W[0]
	wg1 = np.zeros(len(ra_targeted))
	for k in range(len(ra_targeted)):
            indx = (ra_targeted[k] - 120)/(130.0/numra)
            indy = (np.sin(dec_targeted[k]*3.1415/180) + 0.17364817766693033)/-(0.76604444311/numdec) 
            indx = int(indx)
            indy = int(indy)
            wg1[k] = ws[indy][indx]
        wg2 = 1.0 / ( 1.0 + 3000.0 * nis)
        wg = wg1*wg2
        h5z = newt.create_dataset('weight',data=wg)
        h5z.attrs['Definition']="Weight"
    else:
        pass
    if debug==True:
        print 'Weight distribution for targeted'
        xp = np.linspace(0,1.6,100)
        idxran = np.random.random_integers(0,high=len(wg)-1,size=3000000)
        xran = redshift_targeted[idxran]
        yran = wg[idxran]
        fig = plt.figure()
        plt.plot(xran,yran,'+')
        plt.xlabel("Redshift")
        plt.ylabel("Weight per galaxy")
        fig.suptitle('Galaxy weights as a function of redshift for targeted galaxies in Aardvark at r = 19.'+str(magnitude))
        plt.show()
    else:
        pass
    newt.close() 

    ####################################################
    # Now add the weights to the corresponding randoms #
    ####################################################
    if skip_randoms==True:
        pass
    else:
        if weights==0:
            pass
        else:
            newr_targ = h5py.File(newrname_targ, "r+")
            print 'Adding weights to targeted randoms'
            zz2 = newr_targ['z'][...]
            if same_weights==True:
                nis2 = np.power(10.0,yfit(zz2))
            else:
                nis2 = np.power(10.0,yfit2(zz2))
            if (weights==1)|(weights==3):
                w2 = np.ones(len(zz2))
            elif (weights==2)|(weights==4):
                w2 = 1.0 / ( 1.0 + 3000.0 * nis2 )
            else:
                pass
            if debug==True:
                print 'Weight distibution for randoms - Targeted'
                xp = np.linspace(0,1.6,100)
                idxran = np.random.random_integers(0,high=len(w2)-1,size=2000000)
                xran = zz2[idxran]
                yran = w2[idxran]
                fig = plt.figure()
                plt.plot(xran,yran,'+')
                plt.xlabel("Redshift")
                plt.ylabel("Weight per galaxy")
                fig.suptitle('Galaxy weights as a function of redshift for all galaxies in Aardvark at r = 19.'+str(magnitude))
                plt.show()
            else:
                pass
            if skip_randoms==True:
                pass
            else:
                if (weights==1)|(weights==2)|(weights==3)|(weights==4):
                    h5w2 = newr_targ.create_dataset('weight',data=w2)
                    h5w2.attrs['Definition']="weight"
                    time.sleep(5)
                else: 
                    pass
            newr_targ.close()
