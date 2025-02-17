/*----------------------------------------------------------------------------
 | Authors:
 |     - Weizhi Xu   (weizhix)  752454
 |     - Zijun Chen  (zijunc3)  813190
 -----------------------------------------------------------------------------*/
#pragma once

#ifdef WITH_CUDA

#include <stdio.h>
#include "utils.hpp"

void nbody_cuda_single_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);

void nbody_cuda_mpi_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);

#endif
