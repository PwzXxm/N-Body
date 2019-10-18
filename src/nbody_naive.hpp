#pragma once

#include "utils.hpp"

void nbody_seq_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);
void nbody_mpi_openmp_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output);
