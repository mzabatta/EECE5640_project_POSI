#!/bin/sh
#BSUB -J myproj_omp_30
#BSUB -o omp_output_file_30
#BSUB -e omp_error_file_30
#BSUB -n 2
export OMP_NUM_THREADS=30
./proj2
./proj2
./proj2
./proj2
./proj2
./proj2
