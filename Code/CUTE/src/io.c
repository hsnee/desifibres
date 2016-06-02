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
//                         Read-write routines                       //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"
#include <hdf5.h>

static int n_objects=-1;
static int estimator=-1;
static int input_format=-1;

static int read_line(FILE *fi,double *zz,double *cth,
		     double *phi,double *weight)
{
  double x0,x1,x2;
  char s0[1024];
  int sr;
  if(fgets(s0,sizeof(s0),fi)==NULL) return 1;

#ifdef _WITH_WEIGHTS
  double x3;
  sr=sscanf(s0,"%lf %lf %lf %lf",&x0,&x1,&x2,&x3);
  if(sr!=4) return 1;
  *weight=x3;
#else //_WITH_WEIGHTS
  sr=sscanf(s0,"%lf %lf %lf",&x0,&x1,&x2);
  if(sr!=3) return 1;
  *weight=1;
#endif //_WITH_WEIGHTS
  
  //////
  // Modify here to add other formats
  // x0, x1, x2 are the first columns
  // in the data file
  if(input_format==0) {
    // z  cos(theta)  phi
    if((x1>1)||(x1<-1)) {
      fprintf(stderr,"CUTE: wrong cos(theta) = %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=x1;
    *phi=x2;
  }
  else if(input_format==1) {
    // z  dec  ra
    if((x1<-90)||(x1>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=cos(DTORAD*(90-x1));
    *phi=DTORAD*x2;
  }
  else if(input_format==2) {
    // ra  dec  z
    if((x1<-90)||(x1>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x2;
    *cth=cos(DTORAD*(90-x1));
    *phi=DTORAD*x0;
  }
  else if(input_format==3) {
    // z ra  dec
    if((x2<-90)||(x2>90)) {
      fprintf(stderr,"CUTE: wrong declination: %lf \n",x1);
      return 1;
    }
    *zz=x0;
    *cth=cos(DTORAD*(90-x2));
    *phi=DTORAD*x1;
  }
  else {
    fprintf(stderr,"CUTE: wrong input format %d \n",
	    input_format);
    exit(1);
  }

  if((*zz)<0) {
    fprintf(stderr,"Wrong redshift = %lf \n",(*zz));
    return 1;
  }

  (*phi)=wrap_phi((*phi));
  
  return 0;
}

static void make_CF(histo_t DD,histo_t DR,histo_t RR,
		    np_t sum_wd,np_t sum_wd2,
		    np_t sum_wr,np_t sum_wr2,
		    double *corr,double *ercorr)
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double edd,edr,err;
  double ddd,ddr,drr;
  double norm_dd=0.5*((double)sum_wd*sum_wd-sum_wd2);
  double norm_rr=0.5*((double)sum_wr*sum_wr-sum_wr2);
  double norm_dr=((double)sum_wd)*sum_wr;

  edd=1./sqrt((double)DD);
  edr=1./sqrt((double)DR);
  err=1./sqrt((double)RR);
  ddd=(double)(DD);
  ddr=(double)(DR);
  drr=(double)(RR);

  double c,ec;
  if((DD==0)||(DR==0)||(RR==0)) {
    c=0;
    ec=0;
  }
  else {
    if(estimator==0) { //PH
      c=(norm_rr/norm_dd)*(ddd/drr)-1;
      ec=(1+c)*(edd+err);
    }
    else if(estimator==1) { //DP
      c=(norm_dr/norm_dd)*(ddd/ddr)-1;
      ec=(1+c)*(edd+edr);
    }
    else if(estimator==2) { //HAM
      c=4*(norm_dr/norm_dd)*(norm_dr/norm_rr)*
	(ddd/ddr)*(drr/ddr)-1;
      ec=(1+c)*(edd+edr+err);
    }
    else if(estimator==3) { //LS
      c=(norm_rr/norm_dd)*(ddd/drr)-2*(norm_rr/norm_dr)*(ddr/drr)+1;
      ec=(1+c)*(edd+edr+err);
    }
    else if(estimator==4) { //HEW
      c=(norm_rr/norm_dd)*(ddd/drr)-(norm_rr/norm_dr)*(ddr/drr);
      ec=(1+c)*(edd+edr+err);
    }
    else {
      fprintf(stderr,"WTF\n");
      exit(1);
    }
  }

  *corr=c;
  *ercorr=ec;
}

void write_CF(char *fname,
	      histo_t *DD,histo_t *DR,histo_t *RR,
	      np_t sum_wd,np_t sum_wd2,
	      np_t sum_wr,np_t sum_wr2)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  printf("*** Writing output file ");
#ifdef _VERBOSE
  printf("%s ",fname);
#endif
  printf("\n");

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"Error opening output file %s",fname);
    fprintf(stderr," using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==0) {
    for(ii=0;ii<NB_DZ;ii++) {
      double dz;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
      dz=(ii+0.5)/(NB_DZ*I_DZ_MAX);
      fprintf(fo,"%lE %lE %lE ",dz,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==1) {
    for(ii=0;ii<NB_THETA;ii++) {
      double th;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
#ifdef _LOGBIN
      th=pow(10,((ii+0.5)-NB_THETA)/N_LOGINT+LOG_TH_MAX)/DTORAD;
#else //_LOGBIN
      th=(ii+0.5)/(NB_THETA*I_THETA_MAX*DTORAD);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE ",th,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==2) {
    for(ii=0;ii<NB_R;ii++) {
      double rr;
      double corr,ercorr;
      make_CF(DD[ii],DR[ii],RR[ii],sum_wd,sum_wd2,
	      sum_wr,sum_wr2,&corr,&ercorr);
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE ",rr,corr,ercorr);
#ifdef _WITH_WEIGHTS
      fprintf(fo,"%lE %lE %lE\n",DD[ii],DR[ii],RR[ii]);
#else //_WITH_WEIGHTS
      fprintf(fo,"%llu %llu %llu\n",DD[ii],DR[ii],RR[ii]);
#endif //_WITH_WEIGHTS
    }
  }
  else if(corr_type==3) {
    for(ii=0;ii<NB_RT;ii++) {
      int jj;
      double rt=(ii+0.5)/(NB_RT*I_RT_MAX);
      for(jj=0;jj<NB_RL;jj++) {
	double corr,ercorr;
	double rl=(jj+0.5)/(NB_RL*I_RL_MAX);
	int ind=jj+NB_RL*ii;
	make_CF(DD[ind],DR[ind],RR[ind],sum_wd,sum_wd2,
		sum_wr,sum_wr2,&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE ",rl,rt,corr,ercorr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE\n",DD[ind],DR[ind],RR[ind]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu\n",DD[ind],DR[ind],RR[ind]);
#endif //_WITH_WEIGHTS
      }
    }
  }
  else if(corr_type==4) {
    for(ii=0;ii<NB_R3D;ii++) {
      int jj;
      double rr;
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_R3D)/N_LOGINT+LOG_R3D_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_R3D*I_R3D_MAX);
#endif //_LOGBIN
      for(jj=0;jj<NB_CTH;jj++) {
	double corr,ercorr;
	double mu=(jj+0.5)/(NB_CTH);
	int ind=jj+NB_CTH*ii;
	make_CF(DD[ind],DR[ind],RR[ind],sum_wd,sum_wd2,
		sum_wr,sum_wr2,&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE ",mu,rr,corr,ercorr);
#ifdef _WITH_WEIGHTS
	fprintf(fo,"%lE %lE %lE\n",DD[ind],DR[ind],RR[ind]);
#else //_WITH_WEIGHTS
	fprintf(fo,"%llu %llu %llu\n",DD[ind],DR[ind],RR[ind]);
#endif //_WITH_WEIGHTS
      }
    }
  }
  else if(corr_type==5) {
    for(ii=0;ii<NB_RED;ii++) {
      int jj;
      double z_mean=RED_0+(ii+0.5)/(I_RED_INTERVAL*NB_RED);
      for(jj=0;jj<NB_DZ;jj++) {
	int kk;
	double dz=(jj+0.5)/(NB_DZ*I_DZ_MAX);
	for(kk=0;kk<NB_THETA;kk++) {
	  double theta=(kk+0.5)/(NB_THETA*I_THETA_MAX*DTORAD);
	  int index=kk+NB_THETA*(jj+NB_DZ*ii);
	  double corr,ercorr;
	  make_CF(DD[index],DR[index],RR[index],sum_wd,sum_wd2,
		  sum_wr,sum_wr2,&corr,&ercorr);
	  fprintf(fo,"%lE %lE %lE %lE %lE ",z_mean,dz,theta,corr,ercorr);
#ifdef _WITH_WEIGHTS
	  fprintf(fo,"%lE %lE %lE\n",DD[index],DR[index],RR[index]);
#else //_WITH_WEIGHTS
	  fprintf(fo,"%llu %llu %llu\n",DD[index],DR[index],RR[index]);
#endif //_WITH_WEIGHTS
	}
      }
    }
  }
  fclose(fo);
  
  printf("\n");
}

void write_CF_cuda(char *fname,unsigned long long *DD,
		   unsigned long long *DR,unsigned long long *RR,
		   int nD,int nR)
{
  //////
  // Writes correlation function to file fname
  FILE *fo;
  int ii;

  printf("*** Writing output file ");
#ifdef _VERBOSE
  printf("%s ",fname);
#endif
  printf("\n");

  fo=fopen(fname,"w");
  if(fo==NULL) {
    char oname[64]="output_CUTE.dat";
    fprintf(stderr,"Error opening output file %s",fname);
    fprintf(stderr," using ./output_CUTE.dat");
    fo=fopen(oname,"w");
    if(fo==NULL) error_open_file(oname);
  }

  if(corr_type==1) {
    for(ii=0;ii<NB_HISTO_1D;ii++) {
      double th;
      double corr,ercorr;
      make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),
	      (histo_t)(RR[ii]),(np_t)(nD),(np_t)(nD),
	      (np_t)(nR),(np_t)(nR),&corr,&ercorr);
#ifdef _LOGBIN
      th=pow(10,((ii+0.5)-NB_HISTO_1D)/N_LOGINT+LOG_TH_MAX)/DTORAD;
#else //_LOGBIN
      th=(ii+0.5)/(NB_HISTO_1D*I_THETA_MAX*DTORAD);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
	      th,corr,ercorr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==2) {
    for(ii=0;ii<NB_HISTO_1D;ii++) {
      double rr;
      double corr,ercorr;
      make_CF((histo_t)(DD[ii]),(histo_t)(DR[ii]),
	      (histo_t)(RR[ii]),(np_t)(nD),(np_t)(nD),
	      (np_t)(nR),(np_t)(nR),&corr,&ercorr);
#ifdef _LOGBIN
      rr=pow(10,((ii+0.5)-NB_HISTO_1D)/N_LOGINT+LOG_R_MAX);
#else //_LOGBIN
      rr=(ii+0.5)/(NB_HISTO_1D*I_R_MAX);
#endif //_LOGBIN
      fprintf(fo,"%lE %lE %lE %llu %llu %llu\n",
	      rr,corr,ercorr,DD[ii],DR[ii],RR[ii]);
    }
  }
  else if(corr_type==3) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double rt=(ii+0.5)/(NB_HISTO_2D*I_RT_MAX);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double corr,ercorr;
	double rl=(jj+0.5)/(NB_HISTO_2D*I_RL_MAX);
	int ind=jj+NB_HISTO_2D*ii;
	make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),
		(histo_t)(RR[ind]),(np_t)(nD),(np_t)(nD),
		(np_t)(nR),(np_t)(nR),&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE %llu %llu %llu\n",
		rl,rt,corr,ercorr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  else if(corr_type==4) {
    for(ii=0;ii<NB_HISTO_2D;ii++) {
      int jj;
      double mu=(ii+0.5)/(NB_HISTO_2D);
      for(jj=0;jj<NB_HISTO_2D;jj++) {
	double corr,ercorr;
	double rr;
#ifdef _LOGBIN
	rr=pow(10,((jj+0.5)-NB_HISTO_2D)/N_LOGINT+LOG_R3D_MAX);
#else //_LOGBIN
	rr=(jj+0.5)/(NB_HISTO_2D*I_R3D_MAX);
#endif //_LOGBIN
	int ind=jj+NB_HISTO_2D*ii;
	make_CF((histo_t)(DD[ind]),(histo_t)(DR[ind]),
		(histo_t)(RR[ind]),(np_t)(nD),(np_t)(nD),
		(np_t)(nR),(np_t)(nR),&corr,&ercorr);
	fprintf(fo,"%lE %lE %lE %lE %llu %llu %llu\n",
		mu,rr,corr,ercorr,DD[ind],DR[ind],RR[ind]);
      }
    }
  }
  fclose(fo);
  
  printf("\n");
}

static void check_params(void)
{
  //////
  // Check all parameters are there and are sensible

  //input format
  if(input_format==-1) {
    fprintf(stderr,"CUTE: input format was not provided \n");
    exit(1);
  }

  //vital files
  if(!strcmp(fnameData,"default")) {
    fprintf(stderr,"CUTE: Data catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameOut,"default")) {
    fprintf(stderr,"CUTE: Output filename was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameRandom,"default")) {
    fprintf(stderr,"CUTE: Random catalog was not provided \n");
    exit(1);
  }
  if(!strcmp(fnameMask,"default")) {
    fprintf(stderr,"CUTE: Mask was not provided \n");
    exit(1);
  }

  //numbers of objects from data and random
  if((n_objects!=-1)&&(n_objects<0)) {
    fprintf(stderr,"CUTE: Wrong wumber of lines from the data file. ");
    fprintf(stderr,"All objects in the catalog will be read.");
    n_objects=-1;
  }

  //Random generated?
  if(!strcmp(fnameRandom,"none"))
    gen_ran=1;
  else
    gen_ran=0;

  //# random particles
  if(gen_ran) {
    if(fact_n_rand<1) {
      fprintf(stderr,"CUTE: Number of random particles was not provided, ");
      fprintf(stderr,"using the same as data \n");
      fact_n_rand=1;
    }
  }

  //z distribution
  if(corr_type!=1) {
    if(!strcmp(fnamedNdz,"default")) {
      fprintf(stderr,"CUTE: Redshift distribution was not provided \n");
      exit(1);
    }
  }

  //Cosmological parameters
  if((corr_type!=0)&&(corr_type!=1)&&(corr_type!=5)) {
    if((omega_M<0)||(omega_L<0)||(weos<-5))
      fprintf(stderr,"CUTE: Wrong (or inexistent) cosmological parameters \n");
  }

  //Pixels for radial correlation
  if(corr_type==0) {
    if(aperture_los<=0) {
      fprintf(stderr,"CUTE: wrong aperture for radial 2PCF %lf deg\n",
	      aperture_los/DTORAD);
      fprintf(stderr," using 1 deg \n");
      aperture_los=1.*DTORAD;
    }
  }

  //PM options for angular correlation function
  if(corr_type==1) {
    //PM option
    if((use_pm!=0)&&(use_pm!=1)) {
      fprintf(stderr,"CUTE: No pixel option for angular correlations was given\n");
      exit(1);
    }
    if(use_pm) {
      if(n_side_cth<0) {
	fprintf(stderr,"CUTE: #pixels in spherical coords for angular correlation was not provided.");
	fprintf(stderr," Using 1024 \n");
	n_side_cth=1024;
	n_side_phi=2048;
      }
    }
  }

}

void read_run_params(char *fname)
{
  //////
  // Reads and checks the parameter file
  FILE *fi;
  int n_lin,ii;
  char estim[64]="none";
  
  printf("*** Reading run parameters \n");

  //Read parameters from file
  fi=fopen(fname,"r");
  if(fi==NULL) error_open_file(fname);
  n_lin=linecount(fi);
  rewind(fi);
  for(ii=0;ii<n_lin;ii++) {
    char s0[512],s1[64],s2[256];
    if(fgets(s0,sizeof(s0),fi)==NULL)
      error_read_line(fname,ii+1);
    if((s0[0]=='#')||(s0[0]=='\n')) continue;
    int sr=sscanf(s0,"%s %s",s1,s2);
    if(sr!=2)
      error_read_line(fname,ii+1);

    if(!strcmp(s1,"data_filename="))
      sprintf(fnameData,"%s",s2);
    else if(!strcmp(s1,"random_filename="))
      sprintf(fnameRandom,"%s",s2);
    else if(!strcmp(s1,"num_lines=")) {
      if(!strcmp(s2,"all"))
	n_objects=-1;
      else
	n_objects=atoi(s2);
    }
    else if(!strcmp(s1,"input_format="))
      input_format=atoi(s2);
    else if(!strcmp(s1,"output_filename="))
      sprintf(fnameOut,"%s",s2);
    else if(!strcmp(s1,"mask_filename="))
      sprintf(fnameMask,"%s",s2);
    else if(!strcmp(s1,"z_dist_filename="))
      sprintf(fnamedNdz,"%s",s2);
    else if(!strcmp(s1,"corr_estimator=")) {
      sprintf(estim,"%s",s2);
      if(!strcmp(estim,"PH"))
	estimator=0;
      else if(!strcmp(estim,"DP"))
	estimator=1;
      else if(!strcmp(estim,"HAM"))
	estimator=2;
      else if(!strcmp(estim,"LS"))
	estimator=3;
      else if(!strcmp(estim,"HEW"))
	estimator=4;
      else {
	fprintf(stderr,"CUTE: Unknown estimator %s, using 'LS'\n",estim);
	estimator=3;
      }
    }
    else if(!strcmp(s1,"corr_type=")) {
      if(!strcmp(s2,"radial")) corr_type=0;
      else if(!strcmp(s2,"angular")) corr_type=1;
      else if(!strcmp(s2,"monopole")) corr_type=2;
      else if(!strcmp(s2,"3D_ps")) corr_type=3;
      else if(!strcmp(s2,"3D_rm")) corr_type=4;
      else if(!strcmp(s2,"full")) corr_type=5;
      else {
	fprintf(stderr,"CUTE: wrong corr type %s.",s2);
	fprintf(stderr," Possible types are \"radial\", \"angular\", \"full\",");
	fprintf(stderr," \"monopole\", \"3D_ps\" and \"3D_rm\".\n");
      }
    }
    else if(!strcmp(s1,"np_rand_fact="))
      fact_n_rand=atoi(s2);
    else if(!strcmp(s1,"omega_M="))
      omega_M=atof(s2);
    else if(!strcmp(s1,"omega_L="))
      omega_L=atof(s2);
    else if(!strcmp(s1,"w="))
      weos=atof(s2);
    else if(!strcmp(s1,"radial_aperture="))
      aperture_los=atof(s2)*DTORAD;
    else if(!strcmp(s1,"use_pm="))
      use_pm=atoi(s2);
    else if(!strcmp(s1,"n_pix_sph=")) {
      n_side_cth=atoi(s2);
      n_side_phi=2*n_side_cth;
    }
    else
      fprintf(stderr,"CUTE: Unknown parameter %s\n",s1);
  }
  fclose(fi);

  check_params();

#ifdef _VERBOSE
  printf("  Using estimator: %s\n",estim);
  if(gen_ran) {
    printf("  The random catalog will be generated ");
    if(fact_n_rand==1)
      printf("with as many particles as in the data \n");
    else
      printf("with %d times more particles than the data \n",fact_n_rand);
  }
#endif //_VERBOSE

  printf("\n");
}

Catalog read_catalog(char *fname,np_t *sum_w,np_t *sum_w2)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii; 
  hid_t file_id, dset_id, dspace_id; 
  hsize_t dims[10];
  char *dsetname;
  double z_mean=0;
  Catalog cat;

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");
  
// input_format=5 means use hdf5 otherwise stanard code
  if(input_format==5) {
    fprintf(stderr,"CUTE: Reading HDF file \n");
    /* Open the input file */
    if((file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
      error_open_file(fname); 
    /* Open the dataset with the redshifts */
    dsetname="z";
    if((dset_id = H5Dopen(file_id, dsetname, H5P_DEFAULT)) < 0)
      {
	fprintf(stderr,"CUTE: Unable to open HDF5 dataset!\n");
	fprintf(stderr,"      Dataset name: z \n");
	fprintf(stderr,"      File name:    %s\n", fname);
	exit(1); 
      }
    /* Get data set size and allocate memory */
    dspace_id = H5Dget_space(dset_id);
    H5Sget_simple_extent_dims(dspace_id, dims, NULL); 
    cat.np= dims[0];
    cat.red=(double *)malloc(cat.np*sizeof(double));
    if(cat.red==NULL) error_mem_out();
    cat.cth=(double *)malloc(cat.np*sizeof(double));
    if(cat.cth==NULL) error_mem_out();
    cat.phi=(double *)malloc(cat.np*sizeof(double));
    if(cat.phi==NULL) error_mem_out();
#ifdef _WITH_WEIGHTS
    cat.weight=(double *)malloc(cat.np*sizeof(double));
    if(cat.weight==NULL) error_mem_out();
#endif //_WITH_WEIGHTS
    /* Read the redshifts */
    if(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			 cat.red) < 0)
		{
		  fprintf(stderr,"CUTE: failed to read HDF5 dataset\n");
		  fprintf(stderr,"      Dataset name: z \n");
		  fprintf(stderr,"      File name:    %s\n", fname);
		  exit(1);   
		}
    /* Open the and read dataset with the dec values and convert to costheta */
    dsetname="dec";
    if((dset_id = H5Dopen(file_id, dsetname, H5P_DEFAULT)) < 0)
      {
	fprintf(stderr,"CUTE: Unable to open HDF5 dataset!\n");
	fprintf(stderr,"      Dataset name: dec \n");
	fprintf(stderr,"      File name:    %s\n", fname);
	exit(1); 
      }
    if(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			 cat.cth) < 0)
		{
		  fprintf(stderr,"CUTE: failed to read HDF5 dataset\n");
		  fprintf(stderr,"      Dataset name: costheta \n");
		  fprintf(stderr,"      File name:    %s\n", fname);
		  exit(1);   
		}
    /* Open the and read dataset with the phi=ra values */
    dsetname="ra";
    if((dset_id = H5Dopen(file_id, dsetname, H5P_DEFAULT)) < 0)
      {
	fprintf(stderr,"CUTE: Unable to open HDF5 dataset!\n");
	fprintf(stderr,"      Dataset name: ra \n");
	fprintf(stderr,"      File name:    %s\n", fname);
	exit(1); 
      }
    if(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			 cat.phi) < 0)
		{
		  fprintf(stderr,"CUTE: failed to read HDF5 dataset\n");
		  fprintf(stderr,"      Dataset name: ra \n");
		  fprintf(stderr,"      File name:    %s\n", fname);
		  exit(1);   
		}
#ifdef _WITH_WEIGHTS
    /* Open the and read dataset with the phi values */
    dsetname="weight";
    if((dset_id = H5Dopen(file_id, dsetname, H5P_DEFAULT)) < 0)
      {
	fprintf(stderr,"CUTE: Unable to open HDF5 dataset!\n");
	fprintf(stderr,"      Dataset name: weight \n");
	fprintf(stderr,"      File name:    %s\n", fname);
	exit(1); 
      }
    if(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			 cat.weight) < 0)
		{
		  fprintf(stderr,"CUTE: failed to read HDF5 dataset\n");
		  fprintf(stderr,"      Dataset name: weight \n");
		  fprintf(stderr,"      File name:    %s\n", fname);
		  exit(1);   
		}
#endif //_WITH_WEIGHTS	
    /* convert to radians wrap the phi angles and accumulate weights */
    ng=cat.np;
    int i_dat=0;
    *sum_w=0;
    *sum_w2=0;
    for(ii=0;ii<ng;ii++) {
      cat.phi[i_dat]=wrap_phi(cat.phi[i_dat]*DTORAD);
      cat.cth[i_dat]=sin(cat.cth[i_dat]*DTORAD);
#ifdef _WITH_WEIGHTS
      (*sum_w)+=cat.weight[i_dat];
      (*sum_w2)+=cat.weight[i_dat]*cat.weight[i_dat];
#else //_WITH_WEIGHTS
      (*sum_w)++;
      (*sum_w2)++;
#endif //_WITH_WEIGHTS
      i_dat++;
    }
  }
  else{ // original code


  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    ng=linecount(fd);
  else
    ng=n_objects;
  rewind(fd);
  printf("  %d lines in the catalog\n",ng);

  //Allocate catalog memory
  cat.np=ng;
  cat.red=(double *)malloc(cat.np*sizeof(double));
  if(cat.red==NULL) error_mem_out();
  cat.cth=(double *)malloc(cat.np*sizeof(double));
  if(cat.cth==NULL) error_mem_out();
  cat.phi=(double *)malloc(cat.np*sizeof(double));
  if(cat.phi==NULL) error_mem_out();
#ifdef _WITH_WEIGHTS
  cat.weight=(double *)malloc(cat.np*sizeof(double));
  if(cat.weight==NULL) error_mem_out();
#endif //_WITH_WEIGHTS

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  *sum_w=0;
  *sum_w2=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,weight;
    int st=read_line(fd,&zz,&cth,&phi,&weight);

    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    if(zz<0) {
      fprintf(stderr,"Wrong redshift = %lf %d\n",zz,ii+1);
      exit(1);
    }
    if((cth>1)||(cth<-1)) {
      fprintf(stderr,"Wrong cos(theta) = %lf %d\n",cth,ii+1);
      exit(1);
    }
    phi=wrap_phi(phi);

    cat.red[i_dat]=zz;
    cat.cth[i_dat]=cth;
    cat.phi[i_dat]=phi;
#ifdef _WITH_WEIGHTS
    cat.weight[i_dat]=weight;
    (*sum_w)+=weight;
    (*sum_w2)+=weight*weight;
#else //_WITH_WEIGHTS
    (*sum_w)++;
    (*sum_w2)++;
#endif //_WITH_WEIGHTS
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  printf("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

#ifdef _WITH_WEIGHTS
  printf("  Effective n. of particles: %lf\n",(*sum_w));
#else //_WITH_WEIGHTS
  printf("  Total n. of particles read: %d\n",(*sum_w));
#endif //_WITH_WEIGHTS
  }
  printf("\n");
  return cat;
}

Catalog_f read_catalog_f(char *fname,int *np)
{
  //////
  // Creates catalog from file fname
  FILE *fd;
  int ng;
  int ii;
  double z_mean=0;
  Catalog_f cat;

  printf("*** Reading catalog ");
#ifdef _VERBOSE
  printf("from file %s",fname);
#endif
  printf("\n");

  //Open file and count lines
  fd=fopen(fname,"r");
  if(fd==NULL) error_open_file(fname);
  if(n_objects==-1) 
    ng=linecount(fd);
  else
    ng=n_objects;
  *np=ng;
  rewind(fd);

  //Allocate catalog memory
  cat.np=ng;
  cat.pos=(float *)malloc(3*cat.np*sizeof(float));
  if(cat.pos==NULL) error_mem_out();

  rewind(fd);
  //Read galaxies in mask
  int i_dat=0;
  for(ii=0;ii<ng;ii++) {
    double zz,cth,phi,rr,sth,dum_weight;
    int st=read_line(fd,&zz,&cth,&phi,&dum_weight);
    if(st) error_read_line(fname,ii+1);
    z_mean+=zz;
    
    sth=sqrt(1-cth*cth);
    if(corr_type!=1)
      rr=z2r(zz);
    else
      rr=1;

    cat.pos[3*i_dat]=(float)(rr*sth*cos(phi));
    cat.pos[3*i_dat+1]=(float)(rr*sth*sin(phi));
    cat.pos[3*i_dat+2]=(float)(rr*cth);
    i_dat++;
  }
  fclose(fd);

  if(i_dat!=ng) {
    fprintf(stderr,"CUTE: Something went wrong !!\n");
    exit(1);
  }

  z_mean/=ng;
#ifdef _VERBOSE
  printf("  The average redshift is %lf\n",z_mean);
#endif //_VERBOSE

  printf("\n");
  return cat;
}
