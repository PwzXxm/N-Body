#include <math.h>
#include <stddef.h>

#include "utils.h"

static const int MAGIC_NUMBER = 9036;

void init_MPI_datatype() {
    // vector_t
    {
        const int nitems = 2;
        int block_lengths[] = {1, 1};
        MPI_Datatype types[] = {MPI_FLOAT, MPI_FLOAT};
        MPI_Aint offsets[nitems];
        offsets[0] = offsetof(vector_t, x);
        offsets[1] = offsetof(vector_t, y);
        MPI_Type_create_struct(nitems, block_lengths, offsets, types, &mpi_vector_t);
        MPI_Type_commit(&mpi_vector_t);
    }
    // mpi_particle_t
    {
        const int nitems = 3;
        int block_lengths[] = {1, 1, 1};
        MPI_Datatype types[] = {mpi_vector_t, mpi_vector_t, MPI_FLOAT};
        MPI_Aint offsets[nitems];
        offsets[0] = offsetof(particle_t, pos);
        offsets[1] = offsetof(particle_t, v);
        offsets[2] = offsetof(particle_t, mass);
        MPI_Type_create_struct(nitems, block_lengths, offsets, types, &mpi_particle_t);
        MPI_Type_commit(&mpi_particle_t);
    }
}

void free_MPI_datatype() {
    MPI_Type_free(&mpi_particle_t);
    MPI_Type_free(&mpi_vector_t);
}

FILE* init_output_file(const char *output_file, int n, int m, float s_time) {
    FILE* fp = fopen(output_file, "w");
    if (!fp) return NULL;
    // identifier 
    fwrite(&MAGIC_NUMBER, sizeof(MAGIC_NUMBER), 1, fp);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(&m, sizeof(m), 1, fp);
    fwrite(&s_time, sizeof(s_time), 1, fp);
    return fp;
}

void output_particle_pos(int n, particle_t parts[], FILE* fp) {
    for (int i = 0; i < n; ++i) {
        fwrite(&parts[i].pos.x, sizeof(float), 1, fp);
        fwrite(&parts[i].pos.y, sizeof(float), 1, fp);
    }
}

void output_particle_pos_v(int n, vector_t positions[], FILE* fp) {
    for (int i = 0; i < n; ++i) {
        fwrite(&positions[i].x, sizeof(float), 1, fp);
        fwrite(&positions[i].y, sizeof(float), 1, fp);
    }
}

// return force on particle 1
vector_t force_between_particle(vector_t pos1, vector_t pos2, float m1, float m2, float grav) {
    float dis2 = powf(pos1.x - pos2.x, 2) + powf(pos1.y - pos2.y, 2);
    float dis = sqrtf(dis2);

    if (dis == 0.0f) return (vector_t) { .x = 0.0f, .y = 0.0f };

    float f = (grav * m1 * m2) / dis2;
    vector_t fv;
    fv.x = f * (pos2.x - pos1.x) / dis;
    fv.y = f * (pos2.y - pos1.y) / dis;
    return fv;
}

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

long double GetTimeSpentInSeconds(uint64_t start) {
    return (long double)((GetTimeStamp() - start) / 1000000.0f);
}
