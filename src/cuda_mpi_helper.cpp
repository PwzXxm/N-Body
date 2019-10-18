#ifdef WITH_CUDA

#include "cuda_mpi_helper.hpp"
#include "utils.hpp"


// Code compiled by nvcc has some problem using MPI.
// All MPI related operations are moved here and compiled using mpic++.
// Both object files will be linked together.

void get_mpi_size_rank(int* size, int* rank){
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
}


void sync_particle_positions(void *buffer, int count, int root) {
    MPI_Bcast(buffer, count, mpi_vector_t, root, MPI_COMM_WORLD);
}

#endif
