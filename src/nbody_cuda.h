#pragma once

#ifdef WITH_CUDA

#include <stdio.h>
#include "utils.h"

void nbody_cuda_single_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);


#endif