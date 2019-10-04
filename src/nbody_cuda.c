#ifdef WITH_CUDA

#include "nbody_cuda.h"

extern void __nbody_cuda_single_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output);

void nbody_cuda_single_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {
    __nbody_cuda_single_naive(n, m, dt, parts, grav, fp, full_output);
}


#endif
