#!/bin/bash
#SBATCH --nodes 4
#SBATCH --ntasks=4
#SBATCH --partition gpgpu
#SBATCH --gres=gpu:4
#SBATCH --time 00:15:00
#SBATCH --cpus-per-task=1

set -e
set -x

module load MPICH/3.2.1-GCC-6.2.0
module load GCC-CUDA/6.2.0-9.0.176

make clean
make nbody-with-cuda

printf "==========================\n"
printf "\n\n"
mpirun ./nbody-with-cuda cuda_single_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun ./nbody-with-cuda cuda_single_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
mpirun ./nbody-with-cuda cuda_single_naive 240 0.01 0 test_data/40000.in tmp.out
printf "\n\n"
mpirun ./nbody-with-cuda cuda_single_naive 240 0.01 0 test_data/80000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 1 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/40000.in tmp.out
printf "\n\n"
mpirun -np 1 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/80000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 2 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/40000.in tmp.out
printf "\n\n"
mpirun -np 2 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/80000.in tmp.out
printf "\n\n"

printf "==========================\n"
printf "\n\n"
mpirun -np 4 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/10000.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/20000.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/40000.in tmp.out
printf "\n\n"
mpirun -np 4 ./nbody-with-cuda cuda_mpi_naive 240 0.01 0 test_data/80000.in tmp.out
printf "\n\n"
