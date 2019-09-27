#pragma once

#include <stdio.h>
#include <stdbool.h>

typedef struct coord {
    float x, y;
} vector_t;

typedef struct particle {
    vector_t pos, v;
    float weight;
} particle_t;

// num of particles
// num of steps
// time of each step
// particle init data
// output file pointer
// whether output data for all steps
void (*algorithm_fun_ptr)(int, int, float, particle_t[], FILE*, bool);

FILE* init_output_file(const char *output_file, int n, int m, float s_time);
void output_particle_pos(int n, particle_t parts[], FILE* fp);


