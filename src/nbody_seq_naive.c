#include <stdlib.h>

#include "nbody_seq_naive.h"

void nbody_seq_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {

    vector_t *acc = (vector_t *) malloc(sizeof(vector_t) * n);

    for (int m_i = 0; m_i < m; ++m_i) {

        // reset acceleration
        for (int i = 0; i < n; ++i) {
            acc[i].x = 0;
            acc[i].y = 0;
        }
        
        // compute acceleration
        // update acc pair by pair to save half calculations
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                vector_t f = force_between_particle(parts[i].pos, parts[j].pos, parts[i].mass, parts[j].mass, grav);
                // printf("fx = %f, fy = %f \n", f.x, f.y);
                acc[i].x += f.x / parts[i].mass;
                acc[i].y += f.y / parts[i].mass;
                acc[j].x -= f.x / parts[j].mass;
                acc[j].y -= f.y / parts[j].mass;
            }
        }

        // update velocity & position
        for (int i = 0; i < n; ++i) {
            // printf("%f %f\n", acc[i].x, acc[i].y);
            parts[i].v.x += acc[i].x * dt;
            parts[i].v.y += acc[i].y * dt;
            parts[i].pos.x += parts[i].v.x * dt;
            parts[i].pos.y += parts[i].v.y * dt;
        }

        // last step or full output mode
        if (m_i == (m - 1) || full_output) {
            output_particle_pos(n, parts, fp);
        }

    }
    free(acc);
}
