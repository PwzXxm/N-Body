#pragma once

#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include <sys/time.h>
#include <inttypes.h>

#define ROOT_NODE 0

typedef struct vector {
    float x, y;
} vector_t;

MPI_Datatype mpi_vector_t;

typedef struct particle {
    vector_t pos, v;
    float mass;
} particle_t;

MPI_Datatype mpi_particle_t;

// num of particles
// num of steps
// time of each step
// particle init data
// output file pointer
// whether output data for all steps
void (*algorithm_fun_ptr)(int, int, float, particle_t[], float, FILE*, bool);

void init_MPI_datatype();
void free_MPI_datatype();


FILE* init_output_file(const char *output_file, int n, int m, float s_time);

void output_particle_pos(int n, particle_t parts[], FILE* fp);
// for cuda
#ifdef __cplusplus 
extern "C" 
#endif
void output_particle_pos_v(int n, vector_t positions[], FILE* fp);

vector_t force_between_particle(vector_t pos1, vector_t pos2, float m1, float m2, float grav);

uint64_t GetTimeStamp();

long double GetTimeSpentInSeconds(uint64_t start);