/*----------------------------------------------------------------------------
 | Authors:
 |     - Weizhi Xu   (weizhix)  752454
 |     - Zijun Chen  (zijunc3)  813190
 -----------------------------------------------------------------------------*/
#pragma once

#include "utils.hpp"

void nbody_seq_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);
void nbody_mpi_openmp_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output);
