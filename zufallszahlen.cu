#include <stdio.h>
#include <stdlib.h>
#include <curand_kernel.h> // CURAND Bibliothek header-Datei
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

   // for (i=0; i<NHIST; i++)
   //    printf("h[%d]: %f\n",i,h[i]);

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

__global__ void zufallszahlen_gpu(float *d_z, int N, curandState *states)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;;
   int k,ix;

   // Initialisieren des Zufallszahlen-Generators
   // Der 'state' (Status) wird für jeden Thread unabhängig gespeichert
   curand_init(tid, 0, 0, &states[tid]);

   // gleichverteilte Zufallszahlen in (0,1]
   ix=tid*N;
   for (k=0; k<N; k++)
      d_z[ix+k] = curand_uniform(&states[tid]);
}

int main(int argc, char **argv)
{
   int N, k, Ngpu;
   float *z, *d_z, *h_z;
   curandState *d_states;

   printf("%s Starting...\n\n", argv[0]);

   N=NUMBERS_PER_THREAD;
   if (argc>1)
   {
      N=atoi(argv[1]);
   }

   // Zufallszahlen auf der CPU
   srand((unsigned int)seconds()); // Initialisieren des Zufallszahlen-Generators
   z=(float*)malloc(N*sizeof(float));
   for (k=0; k<N; k++)
   {
      z[k] = (float)(rand()) / (float)RAND_MAX; // Zufallszahlen in (0,1]
   }
   printf("Histogramm der Zufallszahlen auf der CPU:\n\n");
   histogramm(z,N); // Histogramm

   // Zufallszahlen auf der GPU
   // set up device
   int dev = 0;
   cudaDeviceProp deviceProp;
   CHECK(cudaGetDeviceProperties(&deviceProp, dev));
   printf("Using Device %d: %s\n", dev, deviceProp.name);
   CHECK(cudaSetDevice(dev));

   Ngpu=N*BLOCKSIZE*GRIDSZE;
   h_z=(float*)malloc(Ngpu*sizeof(float));
   CHECK(cudaMalloc((void**)&d_z,Ngpu*sizeof(float)));
   CHECK(cudaMalloc((void**)&d_states,BLOCKSIZE*GRIDSZE*sizeof(curandState)));
   zufallszahlen_gpu<<<GRIDSZE,BLOCKSIZE>>>(d_z,N,d_states);
   CHECK(cudaMemcpy(h_z, d_z, Ngpu*sizeof(float), cudaMemcpyDeviceToHost));
   CHECK(cudaGetLastError());
   printf("Histogramm der Zufallszahlen auf der GPU:\n\n");
   histogramm(h_z,Ngpu); // Histogramm

   free(z);
   free(h_z);
   CHECK(cudaFree(d_z));
   CHECK(cudaFree(d_states));
}
