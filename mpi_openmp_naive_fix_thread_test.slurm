#!/bin/bash
#SBATCH --nodes 8
#SBATCH --ntasks=8
#SBATCH --partition snowy
#SBATCH --time 01:00:00
#SBATCH --cpus-per-task=4

set -e
set -x

module load MPICH/3.2.1-GCC-6.2.0

make clean
make nbody

export OMP_NUM_THREADS=4

printf "==========================\n"
printf "\n\n"
mpirun -np 1 ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 2 ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 4 ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 8 ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun -np 8 ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun -np 8 ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 8 ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
