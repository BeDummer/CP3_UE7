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
  double phi_sqr;
  const double complex h_conj = conj(h);
  double complex phi_conj;
  
  for(k=0; k<nvol; k++)
	{
	  //Summe fuer 1. 2. und 4. Terme
	  phi_sqr = pow(cabs(phi[k]),2);
	  phi_conj = conj(phi[k]);
	  
	  S+=phi_sqr + lambda*pow(phi_sqr - 1,2)
	    -(h_conj*phi[k]+phi_conj*h);
	    
	    
	    for(mu=1; mu<=ndim; mu++)
		{
					S-=kappa*(phi_conj*phi[nn[mu][k]]
					 +conj(phi[nn[mu][k]])*phi[k]);
		}   
	}
 return S;
 }
  
  double complex M(double complex *phi)
  {
    double complex sum = 0.;
    for (int k=0; k<nvol; k++)
      sum += phi[k];
    return (sum/((double) nvol));
  }
  
  double P(double S_tmp)
  {
    return (exp(S_tmp));
  }
  
  void calc_local_dist(double complex *phi, double complex h, double *p)
  {
    int mu;
    double phi_sqr;
    double complex B;
    for (int k=0; k<nvol; k++)
    {
      phi_sqr = pow(cabs(phi[k]),2);
      B = h;
      for(mu=1; mu<=ndim; mu++)
      {
	B += kappa*(phi[nn[mu][k]]+phi[nn[ndim+mu][k]])
      }
      p[k] = exp(conj(B)*phi[k]+B*conj(phi[k])-phi_sqr-lambda*pow((phi_sqr-1),2));
    }
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
