#pragma once

#include "utils.h"

void nbody_seq_naive(int n, int m, float t, particle_t parts[], float grav, FILE* fp, bool full_output);
