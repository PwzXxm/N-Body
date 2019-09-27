#include <math.h>

#include "utils.h"

#define GRAV 6.6742E-11

const int magic_number = 9036;

FILE* init_output_file(const char *output_file, int n, int m, float s_time) {
    FILE* fp = fopen(output_file, "w");
    if (!fp) return NULL;
    // identifier 
    fwrite(&magic_number, sizeof(magic_number), 1, fp);
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

// return force on particle 1
vector_t force_between_particle(vector_t pos1, vector_t pos2, float m1, float m2) {
    float dis = sqrtf(powf(pos1.x - pos2.x, 2) + powf(pos1.y - pos2.y, 2));
    float f = (GRAV * m1 * m2) / powf(dis, 2);
    vector_t fv;
    fv.x = f * (pos2.x - pos1.x) / dis;
    fv.y = f * (pos2.y - pos1.y) / dis;
    return fv;
}
