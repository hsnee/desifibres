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

#ifndef _CUTE_DEFINE_
#define _CUTE_DEFINE_

extern char fnameData[128];
extern char fnameRandom[128];
extern char fnameOut[128];
extern char fnameMask[128];
extern char fnamedNdz[128];

extern int corr_type;

extern double omega_M;
extern double omega_L;
extern double weos;

extern double aperture_los;

extern int use_pm;

extern int fact_n_rand;
extern int gen_ran;

extern int n_side_cth;
extern int n_side_phi;
extern int n_boxes2D;
extern int n_side[3];
extern int n_boxes3D;
extern double l_box[3];

/*                MACROS            */
// Other possible macros
//_DEBUG, _VERBOSE, _TRUE_ACOS, _LOGBIN

#define DTORAD 0.017453292519943295 // x deg = x*DTORAD rad
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) //Maximum of two numbers
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) //Minimum of two numbers
#define ABS(a)   (((a) < 0) ? -(a) : (a)) //Absolute value
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) //min(max(a,low),high)

#ifndef N_LOGINT
#define N_LOGINT 20 //# bins per decade for logarithmic binning
#endif

//Full correlation  (mean redshift range)
#ifndef RED_0
#define RED_0 0.2
#endif
#ifndef I_RED_INTERVAL
#define I_RED_INTERVAL 5.
#endif
#ifndef NB_RED
#define NB_RED 1
#endif 

//Redshift correlations
#ifndef I_DZ_MAX
#define I_DZ_MAX 5. //1/Dz_max
#endif
#ifndef NB_DZ
#define NB_DZ 60 //# bins in Dz for radial correlation
#endif

//Angular correlations
#ifndef I_THETA_MAX
#define I_THETA_MAX 0.9549296585482695 //1/theta_max
#endif
#ifndef LOG_TH_MAX
#define LOG_TH_MAX 0.020028618 //log10(theta_max)
#endif
#ifndef NB_THETA
#define NB_THETA 256 //# bins in theta for angular correlation
#endif

//Monopole
#ifndef I_R_MAX
#define I_R_MAX 0.005 //1/r_max
#endif
#ifndef LOG_R_MAX
#define LOG_R_MAX 2.30103 //log10(r_max)
#endif
#ifndef NB_R
#define NB_R 256 //# bins in r for monopole correlation
#endif

//3-D correlations
#ifndef I_RL_MAX
#define I_RL_MAX 0.005 //1/pi_max
#endif
#ifndef I_RT_MAX
#define I_RT_MAX 0.005 //1/sigma_max
#endif
#ifndef I_R3D_MAX
#define I_R3D_MAX 0.005 //1/r_max
#endif
#ifndef LOG_R3D_MAX
#define LOG_R3D_MAX 2.30103 //log10(r_max)
#endif
#ifndef NB_RL
#define NB_RL 128 //#bins in pi for 3-D anisotropic correlation
#endif
#ifndef NB_RT
#define NB_RT 128 //#bins in sigma for 3-D anisotropic correlation
#endif
#ifndef NB_R3D
#define NB_R3D 128 //#bins in r for 3-D anisotropic correlation
#endif
#ifndef NB_CTH
#define NB_CTH 128 //#bins in mu for 3-D anisotropic correlation
#endif

/////////////////////////////
//Histogram options for CUDA
//<1-D
#define NB_HISTO_1D 256 //#bins for 1-D histograms
//2-D
#ifdef _HISTO_2D_128    //128x128 option for 2D histograms in CUDA
#define NB_HISTO_2D 128 //#bins for 2-D histograms
#define NB_X_BATCH 32   //#columns loaded in each batch (4096/NB_HISTO_2D)
#define NTH_RWS_2D 4    //#thread rows (1024/NB_HISTO_2D)
#endif
#ifdef _HISTO_2D_64     //64x64 option for 2D histograms in CUDA
#define NB_HISTO_2D 64
#define NB_X_BATCH 64
#define NTH_RWS_2D 8
#endif
////////////////////////////

#ifndef RED_COSMO_MAX
#define RED_COSMO_MAX 4.0 //# maximum redshift
#endif

#ifdef _WITH_WEIGHTS
#define N_POS 4
typedef double histo_t;
typedef double np_t;
#else //_WITH_WEIGHTS
#define N_POS 3
typedef unsigned long long histo_t;
typedef int np_t;
#endif //_WITH_WEIGHTS

//Cell for angular 2PCF
typedef struct {
  int bounds[4];
  double pos[3];
} Cell2DInfo;

typedef struct {
  np_t np;
  Cell2DInfo *ci;
} Cell2D;


//Box for angular 2PCF
typedef struct {
  int bounds[4];
  double *pos;
} Box2DInfo;

typedef struct {
  int np;
  Box2DInfo *bi;
} Box2D;


//Pixel for radial 2PCF
typedef struct {
  int bounds[4];
  double *pos;
  double *redshifts;
} RadialPixelInfo;

typedef struct {
  int np;
  RadialPixelInfo *pi;
} RadialPixel;


//Box for 3D 2PCFs
typedef struct {
  int np;
  double *pos;
} Box3D; //3D cell


//Mask cube
typedef struct {
  double z0,zf;
  double cth0,cthf;
  double phi0,phif;
} MaskRegion; //Mask region (cube in z,cos(theta),phi)


//Catalog
typedef struct {
  int np;
  double *red,*cth,*phi;
#ifdef _WITH_WEIGHTS
  double *weight;
#endif //_WITH_WEIGHTS
} Catalog; //Catalog (double precision)

typedef struct {
  int np;
  float *pos;
} Catalog_f; //Catalog (single precision)

#endif //_CUTE_DEFINE_
