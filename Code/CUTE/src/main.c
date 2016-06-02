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
//                               Main                                //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "define.h"
#include "common.h"
#include "woftheta.h"

void read_dr_catalogs(Catalog *cat_d,Catalog *cat_r,
		      np_t *sum_wd,np_t *sum_wd2,
		      np_t *sum_wr,np_t *sum_wr2)
{
  //////
  // Reads or creates random and data catalogs

  Catalog cat_dat,cat_ran;

  cat_dat=read_catalog(fnameData,sum_wd,sum_wd2);
  if(gen_ran) {
    read_mask();
    if(corr_type!=1)
      read_red_dist();
    timer(0);
    cat_ran=mk_random_cat(fact_n_rand*cat_dat.np);
    timer(1);
    end_mask();
    *sum_wr=(np_t)(fact_n_rand*cat_dat.np);
    *sum_wr2=(np_t)(fact_n_rand*cat_dat.np);
  }
  else
    cat_ran=read_catalog(fnameRandom,sum_wr,sum_wr2);

#ifdef _DEBUG
  write_Catalog(cat_dat,"debug_DatCat.dat");
  write_Catalog(cat_ran,"debug_RanCat.dat");
#endif //_DEBUG

  *cat_d=cat_dat;
  *cat_r=cat_ran;
}

void run_full_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  RadialPixel *pixrad_dat,*pixrad_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_RED*NB_DZ*NB_THETA];
  histo_t DR[NB_RED*NB_DZ*NB_THETA];
  histo_t RR[NB_RED*NB_DZ*NB_THETA];

  timer(4);

#ifdef _VERBOSE
  printf("*** Full correlation function: \n");
  printf(" - Redshift range: %.3lf < z_mean < %.3lf\n",
	 RED_0,RED_0+1./I_RED_INTERVAL);
  printf(" - # redshift bins: %d\n",NB_RED);
  printf(" - Redshift bin width: %.3lf\n",1/(I_RED_INTERVAL*NB_RED));

  printf(" - Angular range: %.3lf < theta < %.3lf \n",
	 0.,1./I_THETA_MAX);
  printf(" - # angular bins : %d\n",NB_THETA);
  printf(" - Angular resolution: D(theta) = %.3lf \n",
	 1./(I_DZ_MAX*NB_DZ));

  printf(" - Radial range: %.3lf < Dz < %.3lf \n",
	 0.,1/I_DZ_MAX);
  printf(" - # radial bins : %d\n",NB_DZ);
  printf(" - Radial resolution: D(Dz) = %.3lf \n",
	 1./(I_DZ_MAX*NB_DZ));

  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,5);
  pixrad_dat=mk_RadialPixels_from_Catalog(cat_dat,&indices_dat,&nfull_dat,5);
  free_Catalog(cat_dat);
  pixrad_ran=mk_RadialPixels_from_Catalog(cat_ran,&indices_ran,&nfull_ran,5);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_PixRads(n_boxes2D,pixrad_dat,"debug_PixRadDat.dat");
  write_PixRads(n_boxes2D,pixrad_ran,"debug_PixRadRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_full_bf(nfull_dat,indices_dat,pixrad_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_full_bf(nfull_ran,indices_ran,pixrad_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_full_bf(nfull_dat,indices_dat,
	       pixrad_dat,pixrad_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_RadialPixels(n_boxes2D,pixrad_dat);
  free_RadialPixels(n_boxes2D,pixrad_ran);
  free(indices_dat);
  free(indices_ran);
}

void run_radial_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  RadialPixel *pixrad_dat,*pixrad_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_DZ];
  histo_t DR[NB_DZ];
  histo_t RR[NB_DZ];

  timer(4);

#ifdef _VERBOSE
  printf("*** Radial correlation function: \n");
  printf(" - Range: %.3lf < Dz < %.3lf \n",
	 0.,1/I_DZ_MAX);
  printf(" - #bins: %d\n",NB_DZ);
  printf(" - Resolution: D(Dz) = %.3lf \n",
	 1./(I_DZ_MAX*NB_DZ));
  printf(" - Using a brute-force approach \n");
  printf(" - Colinear galaxies within Dtheta = %.3lf (deg) \n",
	 aperture_los/DTORAD);
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,0);
  pixrad_dat=mk_RadialPixels_from_Catalog(cat_dat,&indices_dat,&nfull_dat,0);
  free_Catalog(cat_dat);
  pixrad_ran=mk_RadialPixels_from_Catalog(cat_ran,&indices_ran,&nfull_ran,0);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_PixRads(n_boxes2D,pixrad_dat,"debug_PixRadDat.dat");
  write_PixRads(n_boxes2D,pixrad_ran,"debug_PixRadRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_rad_bf(nfull_dat,indices_dat,pixrad_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_rad_bf(nfull_ran,indices_ran,pixrad_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_rad_bf(nfull_dat,indices_dat,
	       pixrad_dat,pixrad_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_RadialPixels(n_boxes2D,pixrad_dat);
  free_RadialPixels(n_boxes2D,pixrad_ran);
  free(indices_dat);
  free(indices_ran);
}

void run_angular_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  Box2D *boxes_dat,*boxes_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_THETA];
  histo_t DR[NB_THETA];
  histo_t RR[NB_THETA];

  timer(4);

#ifdef _VERBOSE
  printf("*** Angular correlation function: \n");
  printf(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(I_THETA_MAX*DTORAD));
  printf(" - #bins: %d\n",NB_THETA);
#ifdef _LOGBIN
  printf(" - Logarithmic binning with %d bins per decade\n",
	 N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: D(theta) = %.3lf \n",
	 1./(I_THETA_MAX*NB_THETA*DTORAD));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,1);
  boxes_dat=mk_Boxes2D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes2D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_Boxes2D(n_boxes2D,boxes_dat,"debug_Box2DDat.dat");
  write_Boxes2D(n_boxes2D,boxes_ran,"debug_Box2DRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_ang_bfdd(nfull_dat,indices_dat,boxes_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_ang_bfrr(nfull_ran,indices_ran,boxes_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_ang_bf(nfull_dat,indices_dat,
	       boxes_dat,boxes_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_Boxes2D(n_boxes2D,boxes_dat);
  free_Boxes2D(n_boxes2D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
}

void run_angular_corr_pm(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  Cell2D *cells_dat,*cells_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_THETA];
  histo_t DR[NB_THETA];
  histo_t RR[NB_THETA];

  timer(4);

#ifdef _VERBOSE
  printf("*** Angular correlation function: \n");
  printf(" - Range: %.3lf < theta < %.3lf (deg)\n",
	 0.,1/(I_THETA_MAX*DTORAD));
  printf(" - #bins: %d\n",NB_THETA);
#ifdef _LOGBIN
  printf(" - Logarithmic binning with %d bins per decade\n",
	 N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: D(theta) = %.3lf \n",
	 1./(I_THETA_MAX*NB_THETA*DTORAD));
#endif //_LOGBIN
  printf(" - Using a PM approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_2D_params(cat_dat,cat_ran,1);
  cells_dat=mk_Cells2D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  free(indices_dat);
  cells_ran=mk_Cells2D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  free(indices_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_Cells2D(n_boxes2D,cells_dat,"debug_Cell2DDat.dat");
  write_Cells2D(n_boxes2D,cells_ran,"debug_Cell2DRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  timer(0);
  corr_ang_pm(cells_dat,cells_ran,DD,DR,RR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_Cells2D(n_boxes2D,cells_dat);
  free_Cells2D(n_boxes2D,cells_ran);
}

void run_monopole_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  Box3D *boxes_dat,*boxes_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_R];
  histo_t DR[NB_R];
  histo_t RR[NB_R];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** Monopole correlation function: \n");
  printf(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/I_R_MAX);
  printf(" - #bins: %d\n",NB_R);
#ifdef _LOGBIN
  printf(" - Logarithmic binning with %d bins per decade\n",
	 N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: D(r) = %.3lf Mpc/h\n",
	 1./(I_R_MAX*NB_R));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,2);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_mono_bfdd(nfull_dat,indices_dat,boxes_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_mono_bfrr(nfull_ran,indices_ran,boxes_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_mono_bf(nfull_dat,indices_dat,
	       boxes_dat,boxes_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
}

void run_3d_ps_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  Box3D *boxes_dat,*boxes_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_RT*NB_RL];
  histo_t DR[NB_RT*NB_RL];
  histo_t RR[NB_RT*NB_RL];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** 3D correlation function (pi,sigma): \n");
  printf(" - Range: (%.3lf,%.3lf) < (pi,sigma) < (%.3lf,%.3lf) Mpc/h\n",
	 0.,0.,1/I_RL_MAX,1/I_RT_MAX);
  printf(" - #bins: (%d,%d)\n",NB_RL,NB_RT);
  printf(" - Resolution: (d(pi),d(sigma)) = (%.3lf,%.3lf) Mpc/h\n",
	 1./(I_RL_MAX*NB_RL),1./(I_RT_MAX*NB_RT));
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,3);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_3d_ps_bfdd(nfull_dat,indices_dat,boxes_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_3d_ps_bfrr(nfull_ran,indices_ran,boxes_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_3d_ps_bf(nfull_dat,indices_dat,
		 boxes_dat,boxes_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
}

void run_3d_rm_corr_bf(void)
{
  np_t sum_wd,sum_wd2,sum_wr,sum_wr2;
  Catalog cat_dat,cat_ran;

  Box3D *boxes_dat,*boxes_ran;
  int *indices_dat,*indices_ran;
  int nfull_dat,nfull_ran;

  histo_t DD[NB_R3D*NB_CTH];
  histo_t DR[NB_R3D*NB_CTH];
  histo_t RR[NB_R3D*NB_CTH];

  timer(4);

  set_r_z();

#ifdef _VERBOSE
  printf("*** 3D correlation function (r,mu): \n");
  printf(" - Range: %.3lf < r < %.3lf Mpc/h\n",
	 0.,1/I_R3D_MAX);
  printf(" - #bins: %d\n",NB_R3D);
  printf(" - Range: 0.000 < mu < 1.000\n");
  printf(" - #bins: %d\n",NB_CTH);
#ifdef _LOGBIN
  printf(" - Logarithmic binning with %d bins per decade\n",
	 N_LOGINT);
#else //_LOGBIN
  printf(" - Resolution: d(r) = %.3lf Mpc/h\n",
	 1./(I_R3D_MAX*NB_R3D));
#endif //_LOGBIN
  printf(" - Using a brute-force approach \n");
  printf("\n");
#endif

  read_dr_catalogs(&cat_dat,&cat_ran,
		   &sum_wd,&sum_wd2,&sum_wr,&sum_wr2);
  
  printf("*** Boxing catalogs \n");
  init_3D_params(cat_dat,cat_ran,4);
  boxes_dat=mk_Boxes3D_from_Catalog(cat_dat,&indices_dat,&nfull_dat);
  free_Catalog(cat_dat);
  boxes_ran=mk_Boxes3D_from_Catalog(cat_ran,&indices_ran,&nfull_ran);
  free_Catalog(cat_ran);
  printf("\n");
  
#ifdef _DEBUG
  write_Boxes3D(n_boxes3D,boxes_dat,"debug_Box3DDat.dat");
  write_Boxes3D(n_boxes3D,boxes_ran,"debug_Box3DRan.dat");
#endif //_DEBUG

  printf("*** Correlating \n");
  printf(" - Auto-correlating data \n");
  timer(0);
  auto_3d_rm_bf(nfull_dat,indices_dat,boxes_dat,DD);
  timer(2);
  printf(" - Auto-correlating random \n");
  auto_3d_rm_bf(nfull_ran,indices_ran,boxes_ran,RR);
  timer(2);
  printf(" - Cross-correlating \n");
  cross_3d_rm_bf(nfull_dat,indices_dat,
		 boxes_dat,boxes_ran,DR);
  timer(1);

  printf("\n");
  write_CF(fnameOut,DD,DR,RR,
	   sum_wd,sum_wd2,sum_wr,sum_wr2);

  printf("*** Cleaning up\n");
  free_Boxes3D(n_boxes3D,boxes_dat);
  free_Boxes3D(n_boxes3D,boxes_ran);
  free(indices_dat);
  free(indices_ran);
  end_r_z();
}

int main(int argc,char **argv)
{
  //////
  // Main routine
  int ii;
  char fnameIn[128];
  if(argc!=2) {
    printf("Usage ./CUTE <input file>\n");
    exit(1);
  }
  sprintf(fnameIn,"%s",argv[1]);

  setbuf(stdout, NULL);

  printf("\n");
  printf("-----------------------------------------------------------\n");
  printf("|| CUTE - Correlation Utilities and Two-point Estimation ||\n");
  printf("-----------------------------------------------------------\n\n");

  //Initialize random number generator
#ifdef _DEBUG
  srand(1234);
#else
  srand(time(NULL));
#endif
#ifdef _VERBOSE
  printf("Initializing random number generator\n");
  printf("First random number : %d \n",rand());
#endif

#ifdef _VERBOSE
  //Calculate number of threads
  ii=0;
#pragma omp parallel
  {
#pragma omp atomic
    ii++;
  }
  printf("Using %d threads \n",ii);
#endif
  printf("\n");

  read_run_params(fnameIn);
  
  if(corr_type==0)
    run_radial_corr_bf();
  else if(corr_type==1) {
    if(use_pm==1)
      run_angular_corr_pm();
    else
      run_angular_corr_bf();
  }
  else if(corr_type==2)
    run_monopole_corr_bf();
  else if(corr_type==3)
    run_3d_ps_corr_bf();
  else if(corr_type==4)
    run_3d_rm_corr_bf();
  else if(corr_type==5)
    run_full_corr_bf();
  else {
    fprintf(stderr,"CUTE: wrong correlation type.\n");
    exit(0);
  }
  printf("             Done !!!             \n");

  return 0;
}
