#define DEFINE_GLOBAL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "global.h"
#include "geom_pbc.c"

double kappa;
double lambda;

// Implementierung fuer S
 double S(double complex *phi, double complex h)
{	
  double complex S=0;
  int k;
  int mu;
  
  for(k=0; k<nvol; k++)
	{
	  //Summe fuer 1. 2. und 4. Terme
	  S+=pow(cabs(phi[k]),2)
	    +lambda*pow(pow(cabs(phi[k]),2)-1,2)
	    -(conj(h)*phi[k]+conj(phi[k])*h);
	    
	    
	    for(mu=1; mu<=ndim; mu++)
		{
					S-=kappa*(conj(phi[k])*phi[nn[mu][k]]
					 +conj(phi[nn[mu][k]])*phi[k]);
		}   
	}
 return S;
 }
  
  
  
  int main (int argc, char **argv)
  {
	  // Parameter fuer 1.1
	  ndim = 2;
	  kappa = 8;
	  lambda =16;
	  int N = 10;
	  double S_phi;
	  
	  
	 // Set fuer geom_pbc
	  int k;
	  lsize = (int *) malloc((ndim+1)*sizeof(int));
	  for (k=1; k <= ndim; k++)
	  {
		  lsize[k] = N;
	  }
	  
	  geom_pbc();
	  
	  double complex *phi;
	  phi = (double complex*) malloc(nvol*sizeof(double complex));
	  for (k=0; k<nvol; k++)
	  {
		  phi[k] = (double complex) 0.5 + 0*I;
	  }
	  
	  printf("S= %f \n", S(phi,1 + 0*I));
	  
	  // 1.2: random Spins
	  
	  double randreal, randim;
	for(k=0; k<nvol; k++)
	{
		randreal = (double)(rand() & 0xFF ) * .1;
		randim = (double)(rand() & 0xFF ) * .1;
		phi[k] = (double complex)randreal + randim*I;
	}
	
	double S_temp;
	
	S_temp = S(phi,0); // Speichern fuer Invarianztest
	
	printf("random Spin S = %f \n", S(phi,0));
	printf("Beispielhafte Werte (1. / letzter) von Phi:\n %f + %f * i \n %f + %f * i \n", creal(phi[0]), cimag(phi[0]), creal(phi[nvol-1]), cimag(phi[nvol-1]));
	
	//Phi geht über zu Phi = exp(i*alpha)*phi
	
	double alpha;
	alpha = 2;
	for(k=0; k<nvol; k++) phi[k] *= cexp(alpha*I);
	
	printf("Phi = Phi * exp(i*alpha= %f): S = %f \n", alpha, S(phi,0));
	printf("Differenz zu random Spin S: diff = %f \n", S(phi,0)-S_temp);
	
	
	free(lsize);
	free(nn);
	free(phi); 
  }
