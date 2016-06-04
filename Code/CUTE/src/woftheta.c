// // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// This is a function to be used in CUTE to weigh the galaxy pairs. The weighting is//
// the 1 + w(theta) ratio from the angular pair correlation function. This function //
// expects a header "woftheta.h"       						    //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include <stdio.h>
#include "woftheta.h"
#include <math.h>
#define PI 3.14159265
double woftheta( double* pos1, double* pos2 )
{
double r1sq = pos1[0]*pos1[0] + pos1[1]*pos1[1] + pos1[2]*pos1[2];
double r2sq = pos2[0]*pos2[0] + pos2[1]*pos2[1] + pos2[2]*pos2[2];
double woftheta1 = 0;
double costheta = (pos1[0]*pos2[0] + pos1[1]*pos2[1] + pos1[2]*pos2[2])/pow((r1sq*r2sq),0.5);
double theta = acos(costheta) * 180.0 / PI;
if(theta>10.0)
{woftheta1 = 1;}
else
{
// first, let's locate the band theta falls in.
  double xx[] = {0,  1.26727800e-3,   1.59540850e-3,   2.00850000e-3,        2.52855200e-3,   3.18325850e-3,   4.00748500e-3,      5.04512450e-3,   6.35143550e-3,   7.99598350e-3,         1.00663445e-2,   1.26727800e-2,   1.59540850e-2,    2.00850000e-2,   2.52855200e-2,   3.18325850e-2,         4.00748500e-2,   5.04512450e-2,   6.35143550e-2,       7.99598350e-2,   1.00663445e-1,   1.26727800e-1,         1.59540850e-1,   2.00850000e-1,   2.52855200e-1,  3.18325850e-1,   4.00748500e-1,   5.04512450e-1,         6.35143550e-1,   7.99598350e-1,   1.00663445,       1.26727800,   1.59540850,   2.00850000,         2.52855200,   3.18325850,   4.00748500,       5.04512450,   6.35143550,   7.99598350,  10.10000000};
  double y[] = {1.33970127,  1.34719965,  1.36031203,  1.3304395 ,  1.33838053,        1.32838375,  1.32194599,  1.31354647,  1.30362769,  1.24764488,        1.21515942,  1.19492406,  1.16846552,  1.14384543,  1.11744933,        1.09309394,  1.07283245,  1.05909543,  1.04650091,  1.03700382,        1.02935922,  1.02339281,  1.01839671,  1.01451162,  1.01154421,        1.00908982,  1.00714546,  1.00552694,  1.00394494,  1.00268202,        1.00212392,  1.00165722,  1.00070611,  1.00009969,  0.99977558,        0.99982882,  0.99976603,  0.9998287 ,  0.99985354,  0.99988024};
int n=39;
int jl = 0;
int ju = n;
int jm = 0;
int ascnd = (xx[n]>xx[1]);
while (ju-jl>1)
	{
	jm = (ju+jl)/2;
	if ((theta>xx[jm])==ascnd)
		jl = jm;
	else 
		ju = jm;
	}
unsigned long j = jl;

//now to calculate from interpolation

woftheta1 = (y[j+1] - y[j])/(xx[j+1]-xx[j]) * theta + y[j+1] - (y[j+1]-y[j])/(xx[j+1]-xx[j])* xx[j+1];

}
// now last jl is the location of this x

//if((theta>5.550)&&(theta<5.552))
//{printf("theta is %f, and w(theta) is %f",theta, woftheta1);} 
return woftheta1;
}


