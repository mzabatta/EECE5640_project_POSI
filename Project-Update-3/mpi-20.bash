#!/bin/sh
#BSUB -J mpi_proj_20
#BSUB -o mpi_proj_out_20
#BSUB -e mpi_proj_err_20
#BSUB -n 20
#BSUB -R "span[ptile=4]"

mpirun -np 20 mympi
mpirun -np 20 mympi
mpirun -np 20 mympi