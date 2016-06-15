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
    int k;
    double complex sum = 0.;
    for (k=0; k<nvol; k++)
      sum += phi[k];
    return (sum/((double) nvol));
  }
  
  double P(double S_tmp)
  {
    return (exp(-S_tmp));
  }
  
  void calc_local_dist(double complex *phi, double complex h, double *p)
  {
    int k;
    int mu;
    double phi_sqr;
    double complex B;
    for (k=0; k<nvol; k++)
    {
      phi_sqr = pow(cabs(phi[k]),2);
      B = h;
      for(mu=1; mu<=ndim; mu++)
      {
	B += kappa*(phi[nn[mu][k]]+phi[nn[ndim+mu][k]]);
      }
      p[k] = exp(conj(B)*phi[k]+B*conj(phi[k])-phi_sqr-lambda*pow((phi_sqr-1),2));
    }
  }
  
  int main (int argc, char **argv)
  {
    	  int k;
	  /* Intializes random number generator */
	  time_t t;
	  srand((unsigned) time(&t));
	  
	  // Parameter
	  ndim = 2;
	  kappa = 1.;
	  lambda =1.;
	  int N = 10;  
	  
	 // Set fuer geom_pbc
	  lsize = (int *) malloc((ndim+1)*sizeof(int));
	  for (k=1; k <= ndim; k++)
	  {
		  lsize[k] = N;
	  }
	  
	  geom_pbc();
	  
	  double complex *phi;
	  phi = (double complex*) malloc(nvol*sizeof(double complex));
	  double *p_old;
	  p_old = (double *) malloc(nvol*sizeof(double));
	  double *p_new;
	  p_new = (double *) malloc(nvol*sizeof(double));	  
	  
	  double randreal, randim;
	for(k=0; k<nvol; k++)
	{
		randreal = (double)(rand() & 0xFF ) * .01;
		randim = (double)(rand() & 0xFF ) * .01;
		phi[k] = (double complex)randreal + randim*I;
	}
	
	int randpos = (rand() & 0xFF) % nvol;
	printf("Changed Position: %d\n", randpos);
	
	double S_old = S(phi,0.);
	calc_local_dist(phi,0.,p_old);
	printf("Phi(x)=%f+i*%f  S=%f  p(x)=%e\n", creal(phi[randpos]),cimag(phi[randpos]), S_old, p_old[randpos]);
	randreal = (double)(rand() & 0xFF ) * .01;
	randim = (double)(rand() & 0xFF ) * .01;
	phi[randpos] = (double complex)randreal + randim*I;
	double S_new = S(phi,0.);
	calc_local_dist(phi,0.,p_new);
	printf("Phi'(x)=%f+i*%f  S'=%f  p'(x)=%e\n", creal(phi[randpos]),cimag(phi[randpos]), S_new, p_new[randpos]);
	
	printf("P'/P=%e  p'(x)/p(x)=%e\n",exp(S_old-S_new),(p_new[randpos]/p_old[randpos]));
	double complex M_new = M(phi);
	printf("M=%f + i*%f\n",creal(M_new),cimag(M_new));
	
	free(lsize);
	free(nn);
	free(phi);
	free(p_old);
	free(p_new);
  }
