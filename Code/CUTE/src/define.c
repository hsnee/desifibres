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
//               Common global variables and macros                  //
/*********************************************************************/
#include <stdio.h>
#include <math.h>

////////// Input parameters //////////
///
//File names
char fnameData[128]="default";   //Data catalog filename
char fnameRandom[128]="default"; //Random catalog filename
char fnameOut[128]="default";    //Output filename
char fnameMask[128]="default";   //Mask filename
char fnamedNdz[128]="default";   //z-distribution filename

//Correlation
int corr_type=-1; //Type of CF

//Random catalogs
int gen_ran=1;             //Do we generate randoms?
int fact_n_rand=-1; //n_random_particles=fact_n_rand*n_data_particles

//Pixels for radial correlation
double aperture_los=0;

//Cosmology
double omega_M=-10;
double omega_L=-10;
double weos=-10;

//PM stuff
int use_pm=-1;
///
//////////////////////////////////////

////////// Internal variables //////////
///
//2D boxing variables
int n_side_cth,n_side_phi,n_boxes2D;

//Boxing variables
int n_side[3],n_boxes3D;
double l_box[3];
///
////////////////////////////////////////
