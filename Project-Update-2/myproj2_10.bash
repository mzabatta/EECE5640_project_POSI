#!/bin/sh
#BSUB -J myproj_omp_10
#BSUB -o omp_output_file_10
#BSUB -e omp_error_file_10
#BSUB -n 2
export OMP_NUM_THREADS=10
./proj2
./proj2
./proj2
./proj2
./proj2
./proj2
