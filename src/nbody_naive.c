#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#include "nbody_naive.h"

void nbody_mpi_openmp_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    if (m_rank == ROOT_NODE) {
        printf("MPI_size: %d, OPENMP_threads: %d", m_size, omp_get_max_threads());
    }

    vector_t *pos_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    for (int i = 0; i < n; i++) {
        pos_arr[i] = parts[i].pos;
    }

    // distribute work
    int n_per_node = n / m_size, extra_n = n % m_size;
    int interval_of_nodes[m_size + 1];
    interval_of_nodes[0] = 0;
    for (int i = 0; i < m_size; ++i) {
        // first extra_n nodes will have an extra items`
        interval_of_nodes[i+1] = interval_of_nodes[i] + n_per_node + (i < extra_n);
    }

    // for (int i = 0; i < m_size + 1; ++i) {
    //     printf("%d\n", interval_of_nodes[i]);
    // }

    // [local_l, local_r)
    int local_l = interval_of_nodes[m_rank], local_r = interval_of_nodes[m_rank + 1];
    int local_n = local_r - local_l;
    vector_t *local_acc = (vector_t *) malloc(sizeof(vector_t) * local_n);
    vector_t *local_v = (vector_t *) malloc(sizeof(vector_t) * local_n);
    // copy v
    for (int i = 0; i < local_n; ++i) {
        local_v[i] = parts[i + local_l].v;
    }

    // printf("%d\n", local_l);
    // printf("%d\n", local_r);
    // printf("%d\n", local_n);

    // compute
    for (int m_i = 0; m_i < m; ++m_i) {
        // reset acceleration
        #pragma omp parallel for 
        for (int i = 0; i < local_n; ++i) {
            local_acc[i].x = 0;
            local_acc[i].y = 0;
        }

        // compute acceleration
        #pragma omp parallel for 
        for (int i = 0; i < local_n; ++i) {
            int global_i = i + local_l;
            for (int j = 0; j < n; ++j) {
                if (j == global_i) continue;
                vector_t f = force_between_particle(pos_arr[global_i], pos_arr[j], parts[global_i].mass, parts[j].mass, grav);
                local_acc[i].x += f.x / parts[global_i].mass;
                local_acc[i].y += f.y / parts[global_i].mass;
            }
        }

        

        // update velocity & position
        #pragma omp parallel for 
        for (int i = 0; i < local_n; ++i) {
            int global_i = i + local_l;
            // printf("%d\n", global_i);
            local_v[i].x += local_acc[i].x * dt;
            local_v[i].y += local_acc[i].y * dt;
            pos_arr[global_i].x += local_v[i].x * dt;
            pos_arr[global_i].y += local_v[i].y * dt;
        }

        // pos_arr[0].x = 0;
        // pos_arr[1].x = 0;
        // pos_arr[m_rank].x = 1;
        // // broadcast pos_arr
        // printf("r = %d, 0 = %f, 1 = %f \n", m_rank, pos_arr[0].x, pos_arr[1].x);
        for (int i = 0; i < m_size; ++i) {
            int l = interval_of_nodes[i];
            int r = interval_of_nodes[i + 1];
            MPI_Bcast(pos_arr + l, r - l, mpi_vector_t, i, MPI_COMM_WORLD);
        }
        // printf("r = %d, 0 = %f, 1 = %f \n", m_rank, pos_arr[0].x, pos_arr[1].x);
        // return;

        // last step or full output mode
        if ((m_rank == ROOT_NODE) && (m_i == (m - 1) || full_output)) {
            output_particle_pos_v(n, pos_arr, fp);
        }

    }

    free(pos_arr);
    free(local_acc);
    free(local_v);
}

void nbody_seq_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {

    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    if (m_rank != ROOT_NODE) return;

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
