// MPI parallelization of SAXS code 

// compile with " mpiCC -o mympi main_MPI.cpp"

using namespace std;

#include <stdio.h>
#include "fstream"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <vector>
#include <time.h>
#include "mpi.h"


#define N 2140
#define Q 60

float X[N][3];
float q[Q];
float F[N][Q];
float Iq[Q];

int main(int argc, char **argv)
{

	int numtasks;
	int taskid;

  double start_time, end_time, exe_time;

	void SAXS(int my_proc_chunk, int my_offset);
	
	// import data files
	int i,j,k;

	ifstream input("X.txt");
	for(i = 0; i<N; i++){
		for(j = 0; j<3; j++){
			input>>X[i][j];
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
			input3>>F[i][j];
		}
	}
	input3.close();

	float eucdist = 0;
	float sincin = 0;


  	

  	MPI_Status status;
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  	start_time = MPI_Wtime();

  	int proc_chunk = (Q / numtasks);


  	int offset;
  	int dest;
  	int source;

  	if (taskid == 0)
  	{

  		offset = proc_chunk;

  		for (dest=1; dest<numtasks; dest++) 
  		{
    		MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
    		MPI_Send(&X[0][0], N*3, MPI_FLOAT, dest, 2, MPI_COMM_WORLD);
    		MPI_Send(&q[offset], proc_chunk, MPI_FLOAT, dest, 3, MPI_COMM_WORLD);
    		MPI_Send(&F[0][offset], N*proc_chunk, MPI_FLOAT, dest, 4, MPI_COMM_WORLD);
    		MPI_Send(&Iq[offset], proc_chunk, MPI_FLOAT, dest, 5, MPI_COMM_WORLD);
    		offset = offset + proc_chunk;
  		}


  		offset = 0;
  		SAXS( proc_chunk, offset);

  		for (int i = 1; i<numtasks; i++)
  		{
  			source = i;
  			MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    		MPI_Recv(&X[0][0], N*3, MPI_FLOAT, source, 2, MPI_COMM_WORLD, &status);
    		MPI_Recv(&q[offset], proc_chunk, MPI_FLOAT, source, 3, MPI_COMM_WORLD, &status);
    		MPI_Recv(&F[0][offset], N*proc_chunk, MPI_FLOAT, source, 4, MPI_COMM_WORLD, &status);
    		MPI_Recv(&Iq[offset], proc_chunk, MPI_FLOAT, source, 5, MPI_COMM_WORLD, &status);
  		}
  	
	}


	if (taskid > 0)
	{

		source = 0;

		
    	MPI_Recv(&offset, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
    	MPI_Recv(&X[0][0], N*3, MPI_FLOAT, source, 2, MPI_COMM_WORLD, &status);
    	MPI_Recv(&q[offset], proc_chunk, MPI_FLOAT, source, 3, MPI_COMM_WORLD, &status);
    	MPI_Recv(&F[0][offset], N*proc_chunk, MPI_FLOAT, source, 4, MPI_COMM_WORLD, &status);
    	MPI_Recv(&Iq[offset], proc_chunk, MPI_FLOAT, source, 5, MPI_COMM_WORLD, &status);

    	SAXS( proc_chunk, offset);

    	dest = 0;

    	MPI_Send(&offset, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
    	MPI_Send(&X[0][0], N*3, MPI_FLOAT, dest, 2, MPI_COMM_WORLD);
    	MPI_Send(&q[offset], proc_chunk, MPI_FLOAT, dest, 3, MPI_COMM_WORLD);
    	MPI_Send(&F[0][offset], N*proc_chunk, MPI_FLOAT, dest, 4, MPI_COMM_WORLD);
    	MPI_Send(&Iq[offset], proc_chunk, MPI_FLOAT, dest, 5, MPI_COMM_WORLD);

	}


  MPI_Barrier(MPI_COMM_WORLD);
  end_time = MPI_Wtime();
	MPI_Finalize();

  if (taskid ==0){
    exe_time = end_time - start_time;
    printf("Execution time = %f \n", exe_time); 
  }


}


void SAXS(int my_proc_chunk, int my_offset)
{

float eucdist = 0;
float sincin = 0;


for(int k = my_offset; k<my_offset+my_proc_chunk; k++)
{
	for(int i = 0; i<N; i++)
	{ // 214 times
		for(int j = 0; j<N; j++)
		{ // 214 times
			
			eucdist = 0;
			
			// eucdist = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
			eucdist += ((X[i][1]-X[j][1])*(X[i][1]-X[j][1]));
			eucdist += ((X[i][2]-X[j][2])*(X[i][2]-X[j][2]));
			eucdist += ((X[i][3]-X[j][3])*(X[i][3]-X[j][3]));
			

			eucdist = sqrt(eucdist);

	   		// Iq(k) = Iq(k) + F(i,k)*F(j,k)*sinc(q(k)*eucdist/pi); 
   			// We divide by PI because MATLAB's sinc function has an implicit
   			// PI factor in the sinc computation.

			sincin = q[k]*eucdist;


			if(sincin == 0)
				Iq[k] = Iq[k] + (F[i][k] * F[j][k] * 1);	
			else
				Iq[k] = Iq[k] + (F[i][k] * F[j][k] * sin(sincin)/sincin);	
				
		}
	}
}
}