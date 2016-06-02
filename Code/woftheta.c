// // // // // // // // // // // // // // // // // // // // // // // // // // // //  //
// This is a function to be used in CUTE to weigh the galaxy pairs. The weighting is //
// the 1 + w(theta) ratio from the angular pair correlation function. This function	 //
// expects a header "woftheta.h"													 //
// // // // // // // // // // 								   // // // // // // //  //

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
double xx[] = {  0,   1.41253800e-3,   1.77827900e-3,   2.23872100e-3,   2.81838300e-3,   3.54813400e-3,   4.46683600e-3,   5.62341300e-3,   7.07945800e-3,   8.91250900e-3,   1.12201800e-2,   1.41253800e-2,   1.77827900e-2,   2.23872100e-2,   2.81838300e-2,   3.54813400e-2,   4.46683600e-2,   5.62341300e-2,   7.07945800e-2,   8.91250900e-2,   1.12201800e-1,   1.41253800e-1,   1.77827900e-1,   2.23872100e-1,   2.81838300e-1,   3.54813400e-1,   4.46683600e-1,   5.62341300e-1,   7.07945800e-1,  8.91250900e-1 ,   1.12201800,   1.41253800,   1.77827900,   2.23872100,   2.81838300,   3.54813400,   4.46683600,   5.62341300,   7.07945800,   10.10};
double y[] = { 1.34086794 , 1.34881877,  1.36104718,  1.33141773 , 1.33961938 , 1.32984379,  1.32350764,  1.31510266,  1.30527955,  1.24937412,  1.21707012,  1.1969196,  1.17063995,  1.14601008,  1.11970078,  1.09542368,  1.07519675,  1.06148172,  1.04892469,  1.03941401,  1.0317762 ,  1.02581204,  1.02080322,  1.01690248,  1.01392283,  1.01144776,  1.00947427,  1.00782378,  1.00619241,  1.00487224,  1.00423443,  1.00365341,  1.00253189,  1.00169561,  1.00106032,  1.00073591,  1.00026264,  0.99996944,  0.99977533,  0.99972267};
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


