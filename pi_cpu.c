#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#define NUMBERS_PER_THREAD 64
#define BLOCKSIZE 32
#define GRIDSZE   1
#define NHIST     20

void histogramm(float *x, int N)
{
   float h[NHIST], xr[NHIST], dx;
   int i, k;

   dx=1.0/NHIST;
   for (i=0; i<NHIST; i++)
   {
      h[i]=0.0;
      xr[i]=(float)(i+1) * dx;
   }

   for (k=0; k<N; k++)
   {
      for (i=0; i<NHIST; i++)
      {
         if (x[k]<=xr[i])
         {
            h[i]+=1;
            break;
         }
      }
   }

   for (i=0; i<NHIST; i++)
      h[i]/=((float)N*dx);

    for (i=0; i<NHIST; i++)
       printf("h[%d]: %f\n",i,h[i]);

   printf("N = %d\n",N);
   for (i=2*NHIST; i>=0; i--)
   {
      printf(" ");
      for (k=0; k<NHIST; k++)
      {
         if (h[k]>=(dx*i))
            printf("*");
         else
            printf(" ");
      }
      printf("\n");
   }

   printf(" ");
   for (k=0; k<NHIST; k++)
      printf("-");
   printf("\n");
   printf(" 0");
   for (k=1; k<NHIST-1; k++)
      printf(" ");
   printf("1\n\n");
}

//*********NEU*********
/*    Funktion in_out : nimmt ein zufallvektor von Grösse N der die 
      koordinaten (x,z) ein Schiess repräsentieren. [x_i=z_i, y_i=z_i+N/2-1] 
      Sie nimmt auch ein d_in vektor der sagt ob der schiess i in (1) oder 
      out (0) ist.
      Dann der norm von d_in ist die Zahl von Punkte im Kreis.
*/
//cpu
      float approx_pi_cpu(float *x, float *y, double N)
      {
            int i=0;
            float target=0;
            printf("N=%f\n",N);
            for(i=0;i<N;++i)
            {
                  if(x[i]*x[i]+y[i]*y[i]<=1)
                  {
                      target+=1; 
                  }
            }
            return 4*target/N;
      }



int main(int argc, char **argv)
{
   int N, k, Ngpu;
   float *z;
   float *x;
   float approx_pi=0.0;

   printf("%s Starting...\n\n", argv[0]);

   N=NUMBERS_PER_THREAD;
   if (argc>1)
   {
      N=atoi(argv[1]);
   }

   // Zufallszahlen auf der CPU
   srand((unsigned int)seconds()); // Initialisieren des Zufallszahlen-Generators
   z=(float*)malloc(N*sizeof(float));
   x=(float*)malloc(N*sizeof(float));
   for (k=0; k<N; k++)
   {
      z[k] = (float)(rand()) / (float)RAND_MAX; // Zufallszahlen in (0,1]
      x[k] = (float)(rand()) / (float)RAND_MAX;
   }
   //printf("Histogramm der Zufallszahlen auf der CPU:\n\n");
   //histogramm(z,N); // Histogramm

   approx_pi=approx_pi_cpu(x,z,N);
   printf("Approx von pi : %f \n",approx_pi);


   free(z);
   free(x);
   
}