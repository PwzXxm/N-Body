#!/bin/bash
#SBATCH --nodes 4
#SBATCH --ntasks=4
#SBATCH --partition snowy
#SBATCH --time 00:40:00
#SBATCH --cpus-per-task=16

set -e
set -x

module load MPICH/3.2.1-GCC-6.2.0

make clean
make nbody

printf "==========================\n"
export OMP_NUM_THREADS=2
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
export OMP_NUM_THREADS=4
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
export OMP_NUM_THREADS=8
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"

printf "==========================\n"
export OMP_NUM_THREADS=16
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/2500.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/5000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun ./nbody mpi_openmp_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
