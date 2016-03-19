#!/bin/sh
#BSUB -J myproj_omp_20
#BSUB -o omp_output_file_20
#BSUB -e omp_error_file_20
#BSUB -n 2
export OMP_NUM_THREADS=20
./proj2
./proj2
./proj2
./proj2
./proj2
./proj2
