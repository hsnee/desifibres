///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//                      Correlators with OpenMP                      //
/*********************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "define.h"
#include "woftheta.h"

void auto_full_bf(int npix_full,int *indices,
		  RadialPixel *pixrad,histo_t hh[])
{
  //////
  // Radial cross-correlator
  int i;

  for(i=0;i<NB_RED*NB_DZ*NB_THETA;i++) 
    hh[i]=0;

#pragma omp parallel default(none)			\
  shared(npix_full,indices,pixrad,hh,n_side_phi)
  {
    int j;
    histo_t hthread[NB_RED*NB_DZ*NB_THETA];
    double cth_aperture=cos(1./I_THETA_MAX);
    for(j=0;j<NB_RED*NB_DZ*NB_THETA;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad[ip1].np;
      RadialPixelInfo *pi1=pixrad[ip1].pi;
      int *bounds=pi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];

	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(pi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_aperture) {
	    double z_mean=0.5*(redshift1+pi1->redshifts[jj]);
	    int iz_mean=(int)((z_mean-RED_0)*I_RED_INTERVAL*NB_RED);
	    if((iz_mean>=0)&&(iz_mean<NB_RED)) {
	      double dz=fabs(redshift1-pi1->redshifts[jj]);
	      int idz=(int)(dz*I_DZ_MAX*NB_DZ);
	      if((idz<NB_DZ)&&(idz>=0)) {
		int ith;
#ifdef _LOGBIN
		if(prod!=1) {
#ifdef _TRUE_ACOS
		  prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
		  prod=1-MIN(1,prod);
		  prod=0.5*log10(2*prod+0.3333333*prod*prod+
				 0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		  ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
		}
		else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
		prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
		prod=1-MIN(1,prod);
		prod=sqrt(2*prod+0.333333333*prod*prod+
			  0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
		if((ith<NB_THETA)&&(ith>=0)) {
		  int index=ith+NB_THETA*(idz+NB_DZ*iz_mean);
#ifdef _WITH_WEIGHTS
		  hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		  hthread[index]++;
#endif //_WITH_WEIGHTS
		}
	      }
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=pixrad[ip2].np;
		RadialPixelInfo *pi2=pixrad[ip2].pi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(pi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_aperture) {
		    double z_mean=0.5*(redshift1+pi2->redshifts[jj]);
		    int iz_mean=(int)((z_mean-RED_0)*I_RED_INTERVAL*NB_RED);
		    if((iz_mean>=0)&&(iz_mean<NB_RED)) {
		      double dz=fabs(redshift1-pi2->redshifts[jj]);
		      int idz=(int)(dz*I_DZ_MAX*NB_DZ);
		      if((idz<NB_DZ)&&(idz>=0)) {
			int ith;
#ifdef _LOGBIN
			if(prod!=1) {
#ifdef _TRUE_ACOS
			  prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
			  prod=1-MIN(1,prod);
			  prod=0.5*log10(2*prod+0.3333333*prod*prod+
					 0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
			  ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
			}
			else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
			prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
			prod=1-MIN(1,prod);
			prod=sqrt(2*prod+0.333333333*prod*prod+
				  0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
			ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
			if((ith<NB_THETA)&&(ith>=0)) {
			  int index=ith+NB_THETA*(idz+NB_DZ*iz_mean);
#ifdef _WITH_WEIGHTS
			  hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[index]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_RED*NB_DZ*NB_THETA;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_full_bf(int npix_full,int *indices,
		   RadialPixel *pixrad1,RadialPixel *pixrad2,
		   histo_t hh[])
{
  //////
  // Radial cross-correlator
  int i;

  for(i=0;i<NB_RED*NB_DZ*NB_THETA;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,pixrad1,pixrad2,hh,n_side_phi)
  {
    int j;
    histo_t hthread[NB_RED*NB_DZ*NB_THETA];
    double cth_aperture=cos(1./I_THETA_MAX);
    for(j=0;j<NB_RED*NB_DZ*NB_THETA;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad1[ip1].np;
      RadialPixelInfo *pi1=pixrad1[ip1].pi;
      int *bounds=pi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad2[ip2].np>0) {
	      int jj;
	      int np2=pixrad2[ip2].np;
	      RadialPixelInfo *pi2=pixrad2[ip2].pi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(pi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_aperture) {
		  double z_mean=0.5*(redshift1+pi2->redshifts[jj]);
		  int iz_mean=(int)((z_mean-RED_0)*I_RED_INTERVAL*NB_RED);
		  if((iz_mean>=0)&&(iz_mean<NB_RED)) {
		    double dz=fabs(redshift1-pi2->redshifts[jj]);
		    int idz=(int)(dz*I_DZ_MAX*NB_DZ);
		    if((idz<NB_DZ)&&(idz>=0)) {
		      int ith;
#ifdef _LOGBIN
		      if(prod!=1) {
#ifdef _TRUE_ACOS
			prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
			prod=1-MIN(1,prod);
			prod=0.5*log10(2*prod+0.3333333*prod*prod+
				       0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
			ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
		      }
		      else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
		      prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
		      prod=1-MIN(1,prod);
		      prod=sqrt(2*prod+0.333333333*prod*prod+
				0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		      ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
		      if((ith<NB_THETA)&&(ith>=0)) {
			int index=ith+NB_THETA*(idz+NB_DZ*iz_mean);
#ifdef _WITH_WEIGHTS
			hthread[index]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[index]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_RED*NB_DZ*NB_THETA;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_rad_bf(int npix_full,int *indices,RadialPixel *pixrad,
		 histo_t hh[])
{
  //////
  // Radial auto-correlator

  int i;
  for(i=0;i<NB_DZ;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,pixrad,hh,n_side_phi,aperture_los)
  {
    int j;
    histo_t hthread[NB_DZ];
    double cth_aperture=cos(aperture_los);

    for(j=0;j<NB_DZ;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad[ip1].np;
      RadialPixelInfo *pi1=pixrad[ip1].pi;
      int *bounds=pi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];

	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(pi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_aperture) {
	    double dz=fabs(redshift1-pi1->redshifts[jj]);
	    int ired=(int)(dz*I_DZ_MAX*NB_DZ);
	    if((ired<NB_DZ)&&(ired>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ired]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[ired]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=pixrad[ip2].np;
		RadialPixelInfo *pi2=pixrad[ip2].pi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(pi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_aperture) {
		    double dz=fabs(redshift1-pi2->redshifts[jj]);
		    int ired=(int)(dz*I_DZ_MAX*NB_DZ);
		    if((ired<NB_DZ)&&(ired>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ired]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[ired]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_DZ;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_rad_bf(int npix_full,int *indices,
		  RadialPixel *pixrad1,RadialPixel *pixrad2,
		  histo_t hh[])
{
  //////
  // Radial cross-correlator
  int i;

  for(i=0;i<NB_DZ;i++) 
    hh[i]=0;

#pragma omp parallel default(none)					\
  shared(npix_full,indices,pixrad1,pixrad2,hh,n_side_phi,aperture_los)
  {
    int j;
    histo_t hthread[NB_DZ];
    double cth_aperture=cos(aperture_los);

    for(j=0;j<NB_DZ;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=pixrad1[ip1].np;
      RadialPixelInfo *pi1=pixrad1[ip1].pi;
      int *bounds=pi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(pi1->pos[N_POS*ii]);
	double redshift1=pi1->redshifts[ii];
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(pixrad2[ip2].np>0) {
	      int jj;
	      int np2=pixrad2[ip2].np;
	      RadialPixelInfo *pi2=pixrad2[ip2].pi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(pi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_aperture) {
		  double dz=fabs(redshift1-pi2->redshifts[jj]);
		  int ired=(int)(dz*I_DZ_MAX*NB_DZ);
		  if((ired<NB_DZ)&&(ired>=0)) {
#ifdef _WITH_WEIGHTS
		    hthread[ired]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		    hthread[ired]++;
#endif //_WITH_WEIGHTS
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_DZ;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_ang_bfdd(int npix_full,int *indices,Box2D *boxes,
		 histo_t hh[])
{
  //////
  // Angular auto-correlator

  int i;
  for(i=0;i<NB_THETA;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,boxes,hh,n_side_phi)
  {
    int j;
    histo_t hthread[NB_THETA];
    double cth_max=cos(1/I_THETA_MAX);

    for(j=0;j<NB_THETA;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes[ip1].np;
      Box2DInfo *bi1=boxes[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(bi1->pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(bi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_max) {
	    int ith;
#ifdef _LOGBIN
	    if(prod!=1) {
#ifdef _TRUE_ACOS
	      prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
	      prod=1-MIN(1,prod);
	      prod=0.5*log10(2*prod+0.3333333*prod*prod+
			     0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	      ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
	    }
	    else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
	    prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
	    prod=1-MIN(1,prod);
	    prod=sqrt(2*prod+0.333333333*prod*prod+
		      0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	    ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
	    if((ith<NB_THETA)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ith]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
	      hthread[ith]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=boxes[ip2].np;
		Box2DInfo *bi2=boxes[ip2].bi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(bi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_max) {
		    int ith;
#ifdef _LOGBIN
		    if(prod!=1) {
#ifdef _TRUE_ACOS
		      prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
		      prod=1-MIN(1,prod);
		      prod=0.5*log10(2*prod+0.3333333*prod*prod+
				     0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		      ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
		    }
		    else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
		    prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
		    prod=1-MIN(1,prod);
		    prod=sqrt(2*prod+0.333333333*prod*prod+
			      0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		    ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
		    if((ith<NB_THETA)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ith]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
		      hthread[ith]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_THETA;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_ang_bfrr(int npix_full,int *indices,Box2D *boxes,
		 histo_t hh[])
{
  //////
  // Angular auto-correlator

  int i;
  for(i=0;i<NB_THETA;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(npix_full,indices,boxes,hh,n_side_phi)
  {
    int j;
    histo_t hthread[NB_THETA];
    double cth_max=cos(1/I_THETA_MAX);

    for(j=0;j<NB_THETA;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes[ip1].np;
      Box2DInfo *bi1=boxes[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	double *pos1=&(bi1->pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double *pos2=&(bi1->pos[N_POS*jj]);
	  double prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  if(prod>cth_max) {
	    int ith;
#ifdef _LOGBIN
	    if(prod!=1) {
#ifdef _TRUE_ACOS
	      prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
	      prod=1-MIN(1,prod);
	      prod=0.5*log10(2*prod+0.3333333*prod*prod+
			     0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	      ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
	    }
	    else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
	    prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
	    prod=1-MIN(1,prod);
	    prod=sqrt(2*prod+0.333333333*prod*prod+
		      0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	    ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
	    if((ith<NB_THETA)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[ith]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int icth;
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes[ip2].np>0) {
	      if(ip2>ip1) {
		int np2=boxes[ip2].np;
		Box2DInfo *bi2=boxes[ip2].bi;
		for(jj=0;jj<np2;jj++) {
		  double *pos2=&(bi2->pos[N_POS*jj]);
		  double prod=pos1[0]*pos2[0]+
		    pos1[1]*pos2[1]+pos1[2]*pos2[2];
		  if(prod>cth_max) {
		    int ith;
#ifdef _LOGBIN
		    if(prod!=1) {
#ifdef _TRUE_ACOS
		      prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
		      prod=1-MIN(1,prod);
		      prod=0.5*log10(2*prod+0.3333333*prod*prod+
				     0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		      ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
		    }
		    else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
		    prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
		    prod=1-MIN(1,prod);
		    prod=sqrt(2*prod+0.333333333*prod*prod+
			      0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		    ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
		    if((ith<NB_THETA)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[ith]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_THETA;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_ang_bf(int npix_full,int *indices,
		  Box2D *boxes1,Box2D *boxes2,
		  histo_t hh[])
{
  //////
  // Angular cross-correlator

  int i;
  for(i=0;i<NB_THETA;i++) 
    hh[i]=0;

#pragma omp parallel default(none)			\
  shared(npix_full,indices,boxes1,boxes2,hh,n_side_phi)
  {
    int j;
    histo_t hthread[NB_THETA];
    double cth_max=cos(1/I_THETA_MAX);

    for(j=0;j<NB_THETA;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<npix_full;j++) {
      int ii;
      int ip1=indices[j];
      int np1=boxes1[ip1].np;
      Box2DInfo *bi1=boxes1[ip1].bi;
      int *bounds=bi1->bounds;
      for(ii=0;ii<np1;ii++) {
	int icth;
	double *pos1=&(bi1->pos[N_POS*ii]);
	for(icth=bounds[0];icth<=bounds[1];icth++) {
	  int iphi;
	  int icth_n=icth*n_side_phi;
	  for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	    int iphi_true=(iphi+n_side_phi)%n_side_phi;
	    int ip2=iphi_true+icth_n;
	    if(boxes2[ip2].np>0) {
	      int jj;
	      int np2=boxes2[ip2].np;
	      Box2DInfo *bi2=boxes2[ip2].bi;
	      for(jj=0;jj<np2;jj++) {
		double *pos2=&(bi2->pos[N_POS*jj]);
		double prod=pos1[0]*pos2[0]+
		  pos1[1]*pos2[1]+pos1[2]*pos2[2];
		if(prod>cth_max) {
		  int ith;
#ifdef _LOGBIN
		  if(prod!=1) {
#ifdef _TRUE_ACOS
		    prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
		    prod=1-MIN(1,prod);
		    prod=0.5*log10(2*prod+0.3333333*prod*prod+
				   0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		    ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
		  }
		  else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
		  prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
		  prod=1-MIN(1,prod);
		  prod=sqrt(2*prod+0.333333333*prod*prod+
			    0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
		  ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN
		  if((ith<NB_THETA)&&(ith>=0)) {
#ifdef _WITH_WEIGHTS
		    hthread[ith]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		    hthread[ith]++;
#endif //_WITH_WEIGHTS
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_THETA;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void corr_ang_pm(Cell2D *cellsD,Cell2D *cellsR,
		 histo_t DD[],histo_t DR[],histo_t RR[])
{
  //////
  // PM angular correlator

  int i;
  for(i=0;i<NB_THETA;i++) {
    DD[i]=0;
    DR[i]=0;
    RR[i]=0;
  }

#pragma omp parallel default(none)					\
  shared(cellsD,cellsR,DD,DR,RR,n_side_phi,n_boxes2D)
  {
    int ip1;
    histo_t DDthread[NB_THETA];
    histo_t DRthread[NB_THETA];
    histo_t RRthread[NB_THETA];
    double cth_max=cos(1/I_THETA_MAX);

    for(ip1=0;ip1<NB_THETA;ip1++) {
      DDthread[ip1]=0;
      DRthread[ip1]=0;
      RRthread[ip1]=0;
    }

#pragma omp for nowait schedule(dynamic)
    for(ip1=0;ip1<n_boxes2D;ip1++) {
      np_t nD1=cellsD[ip1].np;
      np_t nR1=cellsR[ip1].np;
      Cell2DInfo *ci1;
      int *bounds;
      double *pos1;
      int icth;
      if(nD1>0)
	ci1=cellsD[ip1].ci;
      else if(nR1>0)
	ci1=cellsR[ip1].ci;
      else continue;
      bounds=ci1->bounds;
      pos1=ci1->pos;
      for(icth=bounds[0];icth<=bounds[1];icth++) {
	int iphi;
	int icth_n=icth*n_side_phi;
	for(iphi=bounds[2];iphi<=bounds[3];iphi++) {
	  double *pos2;
	  double prod;
	  int iphi_true=(iphi+n_side_phi)%n_side_phi;
	  int ip2=iphi_true+icth_n;
	  
	  if(cellsD[ip2].np>0)
	    pos2=(cellsD[ip2].ci)->pos;
	  else if(cellsR[ip2].np>0)
	    pos2=(cellsR[ip2].ci)->pos;
	  else continue;

	  prod=pos1[0]*pos2[0]+
	    pos1[1]*pos2[1]+pos1[2]*pos2[2];
	  
	  if(prod>cth_max) {
	    int ith;
#ifdef _LOGBIN
	    if(prod!=1) {
#ifdef _TRUE_ACOS
	      prod=log10(acos((MIN(1,prod))));
#else //_TRUE_ACOS
	      prod=1-MIN(1,prod);
	      prod=0.5*log10(2*prod+0.3333333*prod*prod+
			     0.0888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	      ith=(int)(N_LOGINT*(prod-LOG_TH_MAX)+NB_THETA);
	    }
	    else ith=-1;
#else //_LOGBIN
#ifdef _TRUE_ACOS
	    prod=acos((MIN(1,prod)));
#else //_TRUE_ACOS
	    prod=1-MIN(1,prod);
	    prod=sqrt(2*prod+0.333333333*prod*prod+
		      0.08888888889*prod*prod*prod);
#endif //_TRUE_ACOS
	    ith=(int)(prod*NB_THETA*I_THETA_MAX);
#endif //_LOGBIN

	    if((ith<NB_THETA)&&(ith>=0)) {
	      DDthread[ith]+=nD1*cellsD[ip2].np;
	      DRthread[ith]+=nD1*cellsR[ip2].np;
	      RRthread[ith]+=nR1*cellsR[ip2].np;
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(ip1=0;ip1<NB_THETA;ip1++) {
	DD[ip1]+=DDthread[ip1];
	DR[ip1]+=DRthread[ip1];
	RR[ip1]+=RRthread[ip1];
      }
    }

  } //end omp parallel

  for(i=0;i<NB_THETA;i++) {
    DD[i]/=2;
    RR[i]/=2;
  }
}
void auto_mono_bfrr(int nbox_full,int *indices,Box3D *boxes,
		  histo_t hh[])
{
  //////
  // Monopole auto-correlator

  int i;
  for(i=0;i<NB_R;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_R];
    double r2_max=1./(I_R_MAX*I_R_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(1/(I_R_MAX*dx))+1;
    }

    for(j=0;j<NB_R;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    int ir;
#ifdef _LOGBIN
	    if(r2>0)
	      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
	    else
	      ir=-1;
#else //_LOGBIN
	    ir=(int)(sqrt(r2)*I_R_MAX*NB_R);
#endif //_LOGBIN
	    if((ir<NB_R)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
	      hthread[ir]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      int ir;
#ifdef _LOGBIN
		      if(r2>0)
			ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
		      else
			ir=-1;
#else //_LOGBIN
		      ir=(int)(sqrt(r2)*I_R_MAX*NB_R);
#endif //_LOGBIN
		      if((ir<NB_R)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
			hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[ir]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_R;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_mono_bfdd(int nbox_full,int *indices,Box3D *boxes,
		  histo_t hh[])
{
  //////
  // Monopole auto-correlator for data-data pairs

  int i;
  for(i=0;i<NB_R;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_R];
    double r2_max=1./(I_R_MAX*I_R_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(1/(I_R_MAX*dx))+1;
    }

    for(j=0;j<NB_R;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    int ir;
#ifdef _LOGBIN
	    if(r2>0)
	      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
	    else
	      ir=-1;
#else //_LOGBIN
	    ir=(int)(sqrt(r2)*I_R_MAX*NB_R);
#endif //_LOGBIN
	    if((ir<NB_R)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
	      hthread[ir]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
	      hthread[ir]++;
#endif //_WITH_WEIGHTS
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      int ir;
#ifdef _LOGBIN
		      if(r2>0)
			ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
		      else
			ir=-1;
#else //_LOGBIN
		      ir=(int)(sqrt(r2)*I_R_MAX*NB_R);
#endif //_LOGBIN
		      if((ir<NB_R)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
			hthread[ir]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
			hthread[ir]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_R;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_mono_bf(int nbox_full,int *indices,
		   Box3D *boxes1,Box3D *boxes2,
		   histo_t hh[])
{
  //////
  // Monopole cross-correlator

  int i;
  for(i=0;i<NB_R;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_R];
    double r2_max=1./(I_R_MAX*I_R_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(1/(I_R_MAX*dx))+1;
    }

    for(j=0;j<NB_R;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    int ir;
#ifdef _LOGBIN
		    if(r2>0)
		      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
		    else
		      ir=-1;
#else //_LOGBIN
		    ir=(int)(sqrt(r2)*I_R_MAX*NB_R);
#endif //_LOGBIN
		    if((ir<NB_R)&&(ir>=0)) {
#ifdef _WITH_WEIGHTS
		      hthread[ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		      hthread[ir]++;
#endif //_WITH_WEIGHTS
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_R;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_3d_ps_bfdd(int nbox_full,int *indices,Box3D *boxes,
		   histo_t hh[])
{
  //////
  // 3D Sigma-Pi Data-Data Correlator

  int i;
  for(i=0;i<NB_RL*NB_RT;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_RL*NB_RT];
    double r2_max=1./(I_RT_MAX*I_RT_MAX)+1./(I_RL_MAX*I_RL_MAX);
    double rt2_max=1./(I_RT_MAX*I_RT_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

    for(j=0;j<NB_RL*NB_RT;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3],xcm[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  xcm[0]=0.5*(pos1[0]+pos2[0]);
	  xcm[1]=0.5*(pos1[1]+pos2[1]);
	  xcm[2]=0.5*(pos1[2]+pos2[2]);
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
	      sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
	    int irl=(int)(rl*I_RL_MAX*NB_RL);
	    if((irl<NB_RL)&&(irl>=0)) {
	      double rt2=r2-rl*rl;
	      if(rt2<rt2_max) {
		int irt=(int)(sqrt(rt2)*I_RT_MAX*NB_RT);
		if((irt<NB_RT)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
		  hthread[irl+NB_RL*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		  hthread[irl+NB_RL*irt]++;
#endif //_WITH_WEIGHTS
		}
	      }
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3],xcm[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    xcm[0]=0.5*(pos1[0]+pos2[0]);
		    xcm[1]=0.5*(pos1[1]+pos2[1]);
		    xcm[2]=0.5*(pos1[2]+pos2[2]);
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
		      int irl=(int)(rl*I_RL_MAX*NB_RL);
		      if((irl<NB_RL)&&(irl>=0)) {
			double rt2=r2-rl*rl;
			if(rt2<rt2_max) {
			  int irt=(int)(sqrt(rt2)*I_RT_MAX*NB_RT);
			  if((irt<NB_RT)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
			    hthread[irl+NB_RL*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			    hthread[irl+NB_RL*irt]++;
#endif //_WITH_WEIGHTS
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_RL*NB_RT;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}


void auto_3d_ps_bfrr(int nbox_full,int *indices,Box3D *boxes,
		   histo_t hh[])
{
  //////
  // 3D Sigma Pi Random-Random Correlator

  int i;
  for(i=0;i<NB_RL*NB_RT;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_RL*NB_RT];
    double r2_max=1./(I_RT_MAX*I_RT_MAX)+1./(I_RL_MAX*I_RL_MAX);
    double rt2_max=1./(I_RT_MAX*I_RT_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

    for(j=0;j<NB_RL*NB_RT;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3],xcm[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  xcm[0]=0.5*(pos1[0]+pos2[0]);
	  xcm[1]=0.5*(pos1[1]+pos2[1]);
	  xcm[2]=0.5*(pos1[2]+pos2[2]);
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
	      sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
	    int irl=(int)(rl*I_RL_MAX*NB_RL);
	    if((irl<NB_RL)&&(irl>=0)) {
	      double rt2=r2-rl*rl;
	      if(rt2<rt2_max) {
		int irt=(int)(sqrt(rt2)*I_RT_MAX*NB_RT);
		if((irt<NB_RT)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
		  hthread[irl+NB_RL*irt]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
		  hthread[irl+NB_RL*irt]++;
#endif //_WITH_WEIGHTS
		}
	      }
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3],xcm[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    xcm[0]=0.5*(pos1[0]+pos2[0]);
		    xcm[1]=0.5*(pos1[1]+pos2[1]);
		    xcm[2]=0.5*(pos1[2]+pos2[2]);
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
		      int irl=(int)(rl*I_RL_MAX*NB_RL);
		      if((irl<NB_RL)&&(irl>=0)) {
			double rt2=r2-rl*rl;
			if(rt2<rt2_max) {
			  int irt=(int)(sqrt(rt2)*I_RT_MAX*NB_RT);
			  if((irt<NB_RT)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
			    hthread[irl+NB_RL*irt]+=pos1[3]*pos2[3]*woftheta(pos1,pos2);
#else //_WITH_WEIGHTS
			    hthread[irl+NB_RL*irt]++;
#endif //_WITH_WEIGHTS
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_RL*NB_RT;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_3d_ps_bf(int nbox_full,int *indices,
		    Box3D *boxes1,Box3D *boxes2,
		    histo_t hh[])
{
  //////
  // 3D Sigma Pi Cross Correlator.

  int i;
  for(i=0;i<NB_RL*NB_RT;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_RL*NB_RT];
    double r2_max=1./(I_RT_MAX*I_RT_MAX)+1./(I_RL_MAX*I_RL_MAX);
    double rt2_max=1./(I_RT_MAX*I_RT_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

    for(j=0;j<NB_RL*NB_RT;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3],xcm[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  xcm[0]=0.5*(pos1[0]+pos2[0]);
		  xcm[1]=0.5*(pos1[1]+pos2[1]);
		  xcm[2]=0.5*(pos1[2]+pos2[2]);
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    double rl=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
		      sqrt(xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2]);
		    int irl=(int)(rl*I_RL_MAX*NB_RL);
		    if((irl<NB_RL)&&(irl>=0)) {
		      double rt2=r2-rl*rl;
		      if(rt2<rt2_max) {
			int irt=(int)(sqrt(rt2)*I_RT_MAX*NB_RT);
			if((irt<NB_RT)&&(irt>=0)) {
#ifdef _WITH_WEIGHTS
			  hthread[irl+NB_RL*irt]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[irl+NB_RL*irt]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_RL*NB_RT;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void auto_3d_rm_bf(int nbox_full,int *indices,Box3D *boxes,
		   histo_t hh[])
{
  //////
  // Monopole auto-correlator

  int i;
  for(i=0;i<NB_R3D*NB_CTH;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_R3D*NB_CTH];
    double r2_max=1./(I_R3D_MAX*I_R3D_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

    for(j=0;j<NB_R3D*NB_CTH;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	double *pos1=&(boxes[ip1].pos[N_POS*ii]);
	
	int jj;
	for(jj=ii+1;jj<np1;jj++) {
	  double r2;
	  double *pos2=&(boxes[ip1].pos[N_POS*jj]);
	  double xr[3],xcm[3];
	  xr[0]=pos1[0]-pos2[0];
	  xr[1]=pos1[1]-pos2[1];
	  xr[2]=pos1[2]-pos2[2];
	  xcm[0]=0.5*(pos1[0]+pos2[0]);
	  xcm[1]=0.5*(pos1[1]+pos2[1]);
	  xcm[2]=0.5*(pos1[2]+pos2[2]);
	  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	  if(r2<r2_max) {
	    int ir;
#ifdef _LOGBIN
	    if(r2>0)
	      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R3D_MAX)+NB_R3D);
	    else
	      ir=-1;
#else //_LOGBIN
	    ir=(int)(sqrt(r2)*I_R3D_MAX*NB_R3D);
#endif //_LOGBIN
	    if((ir<NB_R3D)&&(ir>=0)) {
	      int icth;
	      if(r2==0) icth=0;
	      else {
		double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
		  sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
		icth=(int)(cth*NB_CTH);
	      }
	      if((icth<NB_CTH)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
		hthread[icth+NB_CTH*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
		hthread[icth+NB_CTH*ir]++;
#endif //_WITH_WEIGHTS
	      }
	    }
	  }
	}

	int iz;
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes[ip2].np>0) {
		if(ip2>ip1) {
		  int np2=boxes[ip2].np;
		  for(jj=0;jj<np2;jj++) {
		    double r2;
		    double *pos2=&(boxes[ip2].pos[N_POS*jj]);
		    double xr[3],xcm[3];
		    xr[0]=pos1[0]-pos2[0];
		    xr[1]=pos1[1]-pos2[1];
		    xr[2]=pos1[2]-pos2[2];
		    xcm[0]=0.5*(pos1[0]+pos2[0]);
		    xcm[1]=0.5*(pos1[1]+pos2[1]);
		    xcm[2]=0.5*(pos1[2]+pos2[2]);
		    r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		    if(r2<r2_max) {
		      int ir;
#ifdef _LOGBIN
		      if(r2>0)
			ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R3D_MAX)+NB_R3D);
		      else
			ir=-1;
#else //_LOGBIN
		      ir=(int)(sqrt(r2)*I_R3D_MAX*NB_R3D);
#endif //_LOGBIN
		      if((ir<NB_R3D)&&(ir>=0)) {
			int icth;
			if(r2==0) icth=0;
			else {
			  double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			    sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
			  icth=(int)(cth*NB_CTH);
			}
			if((icth<NB_CTH)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
			  hthread[icth+NB_CTH*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			  hthread[icth+NB_CTH*ir]++;
#endif //_WITH_WEIGHTS
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_R3D*NB_CTH;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}

void cross_3d_rm_bf(int nbox_full,int *indices,
		    Box3D *boxes1,Box3D *boxes2,
		    histo_t hh[])
{
  //////
  // Monopole auto-correlator

  int i;
  for(i=0;i<NB_R3D*NB_CTH;i++) 
    hh[i]=0;

#pragma omp parallel default(none)				\
  shared(nbox_full,indices,boxes1,boxes2,hh,n_side,l_box)
  {
    int j;
    histo_t hthread[NB_R3D*NB_CTH];
    double r2_max=1./(I_R3D_MAX*I_R3D_MAX);
    int irange[3];
    
    for(j=0;j<3;j++) {
      double dx=l_box[j]/n_side[j];
      irange[j]=(int)(sqrt(r2_max)/dx)+1;
    }

    for(j=0;j<NB_R3D*NB_CTH;j++)
      hthread[j]=0;

#pragma omp for nowait schedule(dynamic)
    for(j=0;j<nbox_full;j++) {
      int ii;
      int ip1=indices[j];

      int np1=boxes1[ip1].np;

      int ix1=ip1%n_side[0];
      int iz1=ip1/(n_side[0]*n_side[1]);
      int iy1=(ip1-ix1-iz1*n_side[0]*n_side[1])/n_side[0];

      int ixmin=MAX(ix1-irange[0],0);
      int ixmax=MIN(ix1+irange[0],n_side[0]-1);
      int iymin=MAX(iy1-irange[1],0);
      int iymax=MIN(iy1+irange[1],n_side[1]-1);
      int izmin=MAX(iz1-irange[2],0);
      int izmax=MIN(iz1+irange[2],n_side[2]-1);

      for(ii=0;ii<np1;ii++) {
	int iz;
	double *pos1=&(boxes1[ip1].pos[N_POS*ii]);
	for(iz=izmin;iz<=izmax;iz++) {
	  int iy;
	  int iz_n=iz*n_side[0]*n_side[1];
	  for(iy=iymin;iy<=iymax;iy++) {
	    int ix;
	    int iy_n=iy*n_side[0];
	    for(ix=ixmin;ix<=ixmax;ix++) {
	      int ip2=ix+iy_n+iz_n;
	      if(boxes2[ip2].np>0) {
		int jj;
		int np2=boxes2[ip2].np;
		for(jj=0;jj<np2;jj++) {
		  double r2;
		  double *pos2=&(boxes2[ip2].pos[N_POS*jj]);
		  double xr[3],xcm[3];
		  xr[0]=pos1[0]-pos2[0];
		  xr[1]=pos1[1]-pos2[1];
		  xr[2]=pos1[2]-pos2[2];
		  xcm[0]=0.5*(pos1[0]+pos2[0]);
		  xcm[1]=0.5*(pos1[1]+pos2[1]);
		  xcm[2]=0.5*(pos1[2]+pos2[2]);
		  r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
		  if(r2<r2_max) {
		    int ir;
#ifdef _LOGBIN
		    if(r2>0)
		      ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R3D_MAX)+NB_R3D);
		    else
		      ir=-1;
#else //_LOGBIN
		    ir=(int)(sqrt(r2)*I_R3D_MAX*NB_R3D);
#endif //_LOGBIN
		    if((ir<NB_R3D)&&(ir>=0)) {
		      int icth;
		      if(r2==0) icth=0;
		      else {
			double cth=fabs(xr[0]*xcm[0]+xr[1]*xcm[1]+xr[2]*xcm[2])/
			  sqrt((xcm[0]*xcm[0]+xcm[1]*xcm[1]+xcm[2]*xcm[2])*r2);
			icth=(int)(cth*NB_CTH);
		      }
		      if((icth<NB_CTH)&&(icth>=0)) {
#ifdef _WITH_WEIGHTS
			hthread[icth+NB_CTH*ir]+=pos1[3]*pos2[3];
#else //_WITH_WEIGHTS
			hthread[icth+NB_CTH*ir]++;
#endif //_WITH_WEIGHTS
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } // end omp for

#pragma omp critical
    {
      for(j=0;j<NB_R3D*NB_CTH;j++)
	hh[j]+=hthread[j];
    }

  } //end omp parallel
}
