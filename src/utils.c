#include "utils.h"


const int magic_number = 9036;

FILE* init_output_file(const char *output_file, int n, int m) {
    FILE* fp = fopen(output_file, "w");
    if (!fp) return NULL;
    // identifier 
    fwrite(&magic_number, sizeof(magic_number), 1, fp);
    fwrite(&n, sizeof(n), 1, fp);
    fwrite(&m, sizeof(m), 1, fp);
    return fp;
}

void output_particle_pos(int n, particle_t parts[], FILE* fp) {
    for (int i = 0; i < n; ++i) {
        fwrite(&parts[i].pos.x, sizeof(float), 1, fp);
        fwrite(&parts[i].pos.y, sizeof(float), 1, fp);
    }
}
