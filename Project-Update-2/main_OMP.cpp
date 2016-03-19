
/*



load('INPUTS.mat')
%load('OUTPUT.mat')
Iq = SAXSIntensity(X,F,q);
plot(q,log10(Iq))
save('OUTPUT.mat')



function Iq = SAXSIntensity(X,F,q)
% INPUTS
% -X is the Nx3 matrix of Atom Coordinates
% -F is the NxQ matrix of Atom Form Factors
% -q is the Qx1 vector of q angle values
% OUTPUTS
% -Iq is the SAXS scattering.
N = length(X);
Q = length(q);
Iq = zeros(Q,1);

for k = 1:Q
    for i = 1:N
        for j = 1:N
            eucdist = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
            Iq(k) = Iq(k) + F(i,k)*F(j,k)*sinc(q(k)*eucdist/pi); 
            %We divide by PI because MATLAB's sinc function has an implicit
            %PI factor in the sinc computation.
        end
    end
end





*/

// compile with g++ -fopenmp main_OMP.cpp -o proj2


using namespace std;

#include <stdio.h>
#include "fstream"
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <omp.h>


int main(int argc, char **argv)
{

	int N = 2140;
	int Q = 60;

	double X[N][3];
	double q[Q];
	double F[N][Q];
	double Iq[Q];
	double temp;
	
	// import data files
	int h,i,j,k;

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

	

	// Begin OMP code


	double eucdist = 0;
	double sincin = 0;


  	///////
	clock_t start, end;
	double cpu_time_used;
	start = omp_get_wtime();
	///////

	int nthreads;

	#pragma omp parallel private(k,i,j,h,eucdist,sincin)
	{
  	#pragma omp for
	for(k = 0; k<Q; k++){
		for(i = 0; i<N; i++){
			
			for(j = 0; j<N; j++){
				
				for(h = 0; h<3; h++){

            		// eucdist = sqrt((X(i,:)-X(j,:))*(X(i,:)-X(j,:))');
					eucdist += ((X[i][h]-X[j][h])*(X[i][h]-X[j][h]));

				}

				eucdist = sqrt(eucdist);

        		// Iq(k) = Iq(k) + F(i,k)*F(j,k)*sinc(q(k)*eucdist/pi); 
            	// We divide by PI because MATLAB's sinc function has an implicit
            	// PI factor in the sinc computation.

				sincin = q[k]*eucdist;
				//printf("k = %d \n",k);

				if(sincin == 0)
					Iq[k] = Iq[k] + (F[i][k] * F[j][k] * 1);	
				else
					Iq[k] = Iq[k] + (F[i][k] * F[j][k] * sin(sincin)/sincin);	
				
				eucdist = 0;
				nthreads = omp_get_num_threads();

			}
		}
		
	}

	}


  ///////
	end = omp_get_wtime();
	cpu_time_used = ((double)(end-start));
	printf("Number of threads = %d\n",nthreads);
	printf("Time taken to run code = %f\n",cpu_time_used);
  ///////


}
