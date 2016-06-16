#include <stdio.h>
#include <stdlib.h>
#include <curand_kernel.h> // CURAND Bibliothek header-Datei
#include "common.h"

#define NUMBERS_PER_THREAD 128
//#define BLOCKSIZE 64
#define GRIDSZE   1
#define NHIST     20


__global__ void zufallszahlen_gpu(float *d_z, int N, curandState *states, float x)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;
   int k,ix;
   
   // Initialisieren des Zufallszahlen-Generators
   // Der 'state' (Status) wird für jeden Thread unabhängig gespeichert
   curand_init(tid, x, 0, &states[tid]);

   // gleichverteilte Zufallszahlen in (0,1]
   ix=tid*N;
   for (k=0; k<N; k++)
      d_z[ix+k] = curand_uniform(&states[tid]);
}

//********* NEU *********
//********* NEU *********
//********* NEU *********

void print_vector_in(int *p, double n)
{
   int j=0;

   printf("In Vector : \n");
   while(j<n)
   {
      printf("  ");
      printf("v(%d) = ",j);
      printf("%d ",p[j]);
      printf("\n");    
      ++j;
   }
}

float approx_pi(int *h_in, int N)
{
  double sum=0;
  int i;
  
  for(i=0; i<N; ++i)
  {
    sum+=h_in[i];
  }
  return 4*sum/N;
}
__global__ void in_out_gpu(float *d_z, float *d_w, int *in, int N)
{
   unsigned int tid  = threadIdx.x + blockDim.x*blockIdx.x;
   in[tid]=0;
   if(tid<N){
      if(d_z[tid]*d_z[tid]+d_w[tid]*d_w[tid]<=1)
      {
      //  printf("IN\n");
        in[tid]=1;
      }
    }
  // else printf("OUT\n");
}

int main(int argc, char **argv)
{
   int N, k, Ngpu;
   float *z, *d_z, *h_z;
   curandState *d_states_z;
   curandState *d_states_w;
   
   N=128;
      if (argc>1)
   {
      N=atoi(argv[1]);
   }
   int BLOCKSIZE=N;

   Ngpu=N*BLOCKSIZE*GRIDSZE;

   
   float *d_w;
   
   int *d_in, *h_in, *d_sum, *h_sum;
   h_in=(int*)malloc(Ngpu*sizeof(int));
   h_sum=(int*)malloc(Ngpu*sizeof(int));
   CHECK(cudaMalloc((void**)&d_in,Ngpu*sizeof(int)));
   CHECK(cudaMalloc((void**)&d_sum,Ngpu*sizeof(int)));

   printf("%s Starting...\n\n", argv[0]);



   // Zufallszahlen auf der CPU
   srand((unsigned int)seconds()); // Initialisieren des Zufallszahlen-Generators
   z=(float*)malloc(N*sizeof(float));
   for (k=0; k<N; k++)
   {
      z[k] = (float)(rand()) / (float)RAND_MAX; // Zufallszahlen in (0,1]
   }
   //printf("Histogramm der Zufallszahlen auf der CPU:\n\n");
   //histogramm(z,N); // Histogramm

   // Zufallszahlen auf der GPU
   // set up device
   int dev = 0;
   cudaDeviceProp deviceProp;
   CHECK(cudaGetDeviceProperties(&deviceProp, dev));
   printf("Using Device %d: %s\n", dev, deviceProp.name);
   CHECK(cudaSetDevice(dev));


   h_z=(float*)malloc(Ngpu*sizeof(float));
   CHECK(cudaMalloc((void**)&d_z,Ngpu*sizeof(float)));
   CHECK(cudaMalloc((void**)&d_states_z,BLOCKSIZE*GRIDSZE*sizeof(curandState)));
   zufallszahlen_gpu<<<GRIDSZE,BLOCKSIZE>>>(d_z,N,d_states_z,(float)(rand()));
   CHECK(cudaMemcpy(h_z, d_z, Ngpu*sizeof(float), cudaMemcpyDeviceToHost));
   CHECK(cudaGetLastError());
   
   CHECK(cudaMalloc((void**)&d_w,Ngpu*sizeof(float)));
   CHECK(cudaMalloc((void**)&d_states_w,BLOCKSIZE*GRIDSZE*sizeof(curandState)));
   zufallszahlen_gpu<<<GRIDSZE,BLOCKSIZE>>>(d_w,N,d_states_w,(float)(rand()));
   CHECK(cudaGetLastError());
   //printf("Histogramm der Zufallszahlen auf der GPU:\n\n");
   //histogramm(h_z,Ngpu); // Histogramm
   
   printf("N = %d\n",N);
   //Anruf von in_out_gpu Funktion
   in_out_gpu<<<GRIDSZE,BLOCKSIZE>>>(d_z,d_w,d_in,N);
   CHECK(cudaDeviceSynchronize());
   CHECK(cudaMemcpy(h_in, d_in, Ngpu*sizeof(int), cudaMemcpyDeviceToHost));
   CHECK(cudaGetLastError());
   //print_vector_in(h_in,N);
   CHECK(cudaMemcpy(h_sum, d_sum, Ngpu*sizeof(int), cudaMemcpyDeviceToHost));
   //printf("In Ereignisse N_in = %d \n",h_sum[0]);
   float ergebnis=approx_pi(h_in,N);
   printf("Approx von Pi : %f \n",ergebnis);

   free(z);
   free(h_z);
   free(h_in);
   free(h_sum);
   CHECK(cudaFree(d_w));
   CHECK(cudaFree(d_z));
   CHECK(cudaFree(d_states_z));
   CHECK(cudaFree(d_states_w));
   CHECK(cudaFree(d_in));
   CHECK(cudaFree(d_sum));
}
