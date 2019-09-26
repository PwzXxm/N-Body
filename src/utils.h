#pragma once

#include <stdio.h>
#include <stdbool.h>

typedef struct coord {
    float x, y;
} coord_t;

typedef struct particle {
    coord_t pos, v;
    float weight;
} particle_t;

// num of particles
// num of steps
// time of each step
// particle init data
// output file pointer
// whether output data for all steps
void (*algorithm_fun_ptr)(int, int, float, particle_t[], FILE*, bool);
