using namespace std;

#include <stdio.h>
#include "fstream"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <cuda.h>
#include <iostream>


__global__ void SAXS(float* X, float* q, float* F, float* Iq)
{
	int N = 2140;
	int k = blockIdx.x*blockDim.x + threadIdx.x;

	float eucdist = 0;
	float sincin = 0;

	if (k < 60){
		for (int i = 0; i<N; i++){
			for(int j = 0; j<N; j++){

				eucdist = 0;
				
				eucdist = eucdist + ((X[i]-X[j])*(X[i]-X[j]));
				eucdist = eucdist +((X[i+N]-X[j+N])*(X[i+N]-X[j+N]));
				eucdist = eucdist + ((X[i+(2*N)]-X[j+(2*N)])*(X[i+(2*N)]-X[j+(2*N)]));

				eucdist = sqrt(eucdist);

				sincin = q[k]*eucdist;
				
				if(sincin == 0)
					Iq[k] = Iq[k] + (F[i+(N*k)] * F[j+(N*k)] * 1);	
				else
					Iq[k] = Iq[k] + (F[i+(N*k)] * F[j+(N*k)] * sin(sincin)/sincin);	

			}

		}
		__syncthreads();	
	}
}


int main(int argc, char **argv)
{
  	cudaEvent_t start;
  	cudaEvent_t stop;
  	float time;
  	cudaEventCreate(&start);
  	cudaEventCreate(&stop);


	int N = 2140;
	int Q = 60;

	float *X, *q, *F, *Iq;
	float *X_d, *q_d, *F_d, *Iq_d;	

	int SIZE_X = N*3*sizeof(float);
	int SIZE_q = Q*sizeof(float);
	int SIZE_F = N*Q*sizeof(float);
	int SIZE_Iq = Q*sizeof(float);

	int DimBlock = 512;	
	int DimGrid = (int)ceil((float)SIZE_Iq/DimBlock);

	X = (float*)malloc(SIZE_X);
	q = (float*)malloc(SIZE_q);
	F = (float*)malloc(SIZE_F);
	Iq = (float*)malloc(SIZE_Iq);

	cudaMalloc(&X_d, SIZE_X);
	cudaMalloc(&q_d, SIZE_q);
	cudaMalloc(&F_d, SIZE_F);
	cudaMalloc(&Iq_d, SIZE_Iq);

	// import data files
	int i,j;

	ifstream input("X.txt");
	for(i = 0; i<N; i++){
		for(j = 0; j<3; j++){
			input>>X[i+(j*N)];
		}
	}
	input.close();

	ifstream input2("q.txt");
	for(i = 0; i<Q; i++){
		input2>>q[i];	
	}
	input2.close();

	ifstream input3("F.txt");
	for(i = 0; i<N; i++){
		for(j = 0; j<Q; j++){
			input3>>F[i+(N*j)];
		}
	}
	input3.close();


	cudaMemcpy(X_d, X, SIZE_X, cudaMemcpyHostToDevice);
	cudaMemcpy(q_d, q, SIZE_q, cudaMemcpyHostToDevice);
	cudaMemcpy(F_d, F, SIZE_F, cudaMemcpyHostToDevice);
	cudaMemcpy(Iq_d, Iq, SIZE_Iq, cudaMemcpyHostToDevice);

  	cudaEventRecord(start,0);

	SAXS<<<DimGrid,DimBlock>>>(X_d, q_d, F_d, Iq_d);

  	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
  	
	cudaMemcpy(Iq, Iq_d, SIZE_Iq, cudaMemcpyDeviceToHost);

	cudaFree(X_d);
	cudaFree(q_d);
	cudaFree(F_d);
	cudaFree(Iq_d);

	cudaEventElapsedTime(&time,start,stop);
	printf("Time for the kernel: %f ms\n", time);

}

