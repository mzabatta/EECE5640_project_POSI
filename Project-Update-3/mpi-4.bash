#!/bin/sh
#BSUB -J mpi_proj_4
#BSUB -o mpi_proj_out_4
#BSUB -e mpi_proj_err_4
#BSUB -n 4
#BSUB -R "span[ptile=4]"

mpirun -np 4 mympi
mpirun -np 4 mympi
mpirun -np 4 mympi