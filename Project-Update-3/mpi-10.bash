#!/bin/sh
#BSUB -J mpi_proj_10
#BSUB -o mpi_proj_out_10
#BSUB -e mpi_proj_err_10
#BSUB -n 10
#BSUB -R "span[ptile=4]"

mpirun -np 10 mympi
mpirun -np 10 mympi
mpirun -np 10 mympi