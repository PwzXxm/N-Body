#include <stdio.h>
#include "utils.hpp"
#include "nbody_cuda.hpp"

__global__ void compute_acceleration_local(vector_t *d_l_acc_arr, vector_t *d_pos_arr, 
                    float *d_m_arr, int n, int local_n, int offset, float grav) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= local_n) return;
    d_l_acc_arr[i].x = 0;
    d_l_acc_arr[i].y = 0;
    int gi = i + offset;
    for (int j = 0; j < n; ++j) {
        if (j == gi) continue;
        // compute force
        float dis2 = powf(d_pos_arr[gi].x - d_pos_arr[j].x, 2) + powf(d_pos_arr[gi].y - d_pos_arr[j].y, 2);
        float dis = sqrtf(dis2);
        float fx = 0, fy = 0;
        if (dis != 0.0f) {
            float f = grav * d_m_arr[gi] * d_m_arr[j] / dis2;
            fx = f * (d_pos_arr[j].x - d_pos_arr[gi].x) / dis;
            fy = f * (d_pos_arr[j].y - d_pos_arr[gi].y) / dis;
        }
        // update acceleration
        d_l_acc_arr[i].x = fx / d_m_arr[gi];
        d_l_acc_arr[i].y = fy / d_m_arr[gi];
    }
}

__global__ void update_velocity_position_local(vector_t *d_l_acc_arr, vector_t *d_l_v_arr, 
                    vector_t *d_pos_arr, float *d_m_arr, int n, int local_n, int offset, float dt) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= local_n) return;
    int gi = i + offset;

    d_l_v_arr[i].x += d_l_acc_arr[i].x * dt;
    d_l_v_arr[i].y += d_l_acc_arr[i].y * dt;
    d_pos_arr[gi].x += d_l_v_arr[i].x * dt;
    d_pos_arr[gi].y += d_l_v_arr[i].y * dt;
}

// The code of following two methods are very similar.
// We duplicate the code instead of reusing them 
// in order to specialize and optimize in the future.
// ===============================================================================
// nbody_cuda_mpi_naive
// ===============================================================================

void nbody_cuda_mpi_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    if (m_rank == ROOT_NODE) {
        printf("MPI_size: %d\n", m_size);
    }

    // distribute work
    int n_per_node = n / m_size, extra_n = n % m_size;
    int interval_of_nodes[m_size + 1];
    interval_of_nodes[0] = 0;
    for (int i = 0; i < m_size; ++i) {
        // first extra_n nodes will have an extra items`
        interval_of_nodes[i+1] = interval_of_nodes[i] + n_per_node + (i < extra_n);
    }

    // [local_l, local_r)
    int local_l = interval_of_nodes[m_rank], local_r = interval_of_nodes[m_rank + 1];
    int local_n = local_r - local_l;
    vector_t *local_acc = (vector_t *) malloc(sizeof(vector_t) * local_n);
    vector_t *local_v = (vector_t *) malloc(sizeof(vector_t) * local_n);


    // // vector_t *acc_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    vector_t *local_v_arr = (vector_t *) malloc(sizeof(vector_t) * local_n);
    vector_t *pos_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    float *m_arr = (float *) malloc(sizeof(float) * n);

    // copy data
    for (int i = 0; i < n; ++i) {
        pos_arr[i] = parts[i].pos;
        m_arr[i] = parts[i].mass;
    }
    for (int i = 0; i < local_n; ++i) local_v_arr[i] = parts[i + local_l].v;

    vector_t *d_local_acc_arr = 0;
    vector_t *d_local_v_arr = 0;
    vector_t *d_pos_arr = 0;
    float *d_m_arr = 0;
    cudaMalloc((void**)&d_local_acc_arr, sizeof(vector_t) * local_n);
    cudaMalloc((void**)&d_local_v_arr, sizeof(vector_t) * local_n);
    cudaMalloc((void**)&d_pos_arr, sizeof(vector_t) * n);
    cudaMalloc((void**)&d_m_arr, sizeof(float) * n);
    
    cudaMemcpy(d_local_v_arr, local_v_arr, sizeof(vector_t) * local_n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_arr, pos_arr, sizeof(vector_t) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_m_arr, m_arr, sizeof(float) * n, cudaMemcpyHostToDevice);

    const int threadsPerBlock = 256;
    int threadsPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    

    for (int m_i = 0; m_i < m; ++m_i) {

        compute_acceleration_local<<<threadsPerGrid, threadsPerBlock>>>
                (d_local_acc_arr, d_pos_arr, d_m_arr, n, local_n, local_l, grav);

        update_velocity_position_local<<<threadsPerGrid, threadsPerBlock>>>
                (d_local_acc_arr, d_local_v_arr, d_pos_arr, d_m_arr, n, local_n, local_l, dt);

        // only sync the changed part
        vector_t *local_pos_arr = pos_arr + local_l;
        vector_t *d_local_pos_arr = d_pos_arr + local_l;
        cudaMemcpy(local_pos_arr, d_local_pos_arr, sizeof(vector_t) * local_n, cudaMemcpyDeviceToHost);
        // broadcast pos_arr
        for (int i = 0; i < m_size; ++i) {
            int l = interval_of_nodes[i];
            int r = interval_of_nodes[i + 1];
            MPI_Bcast(pos_arr + l, r - l, mpi_vector_t, i, MPI_COMM_WORLD);
        }
        cudaMemcpy(d_pos_arr, pos_arr, sizeof(vector_t) * n, cudaMemcpyHostToDevice);

        // last step or full output mode
        if ((m_rank == ROOT_NODE) && (m_i == (m - 1) || full_output)) {
            output_particle_pos_v(n, pos_arr, fp);
        }
    }


    cudaFree(d_local_acc_arr);
    cudaFree(d_local_v_arr);
    cudaFree(d_pos_arr);
    cudaFree(d_m_arr);

    free(local_v_arr);
    free(pos_arr);
    free(m_arr);
}


// ===============================================================================
// nbody_cuda_single_naive
// ===============================================================================

void nbody_cuda_single_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    if (m_rank != ROOT_NODE) return;

    // vector_t *acc_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    vector_t *v_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    vector_t *pos_arr = (vector_t *) malloc(sizeof(vector_t) * n);
    float *m_arr = (float *) malloc(sizeof(float) * n);

    for (int i = 0; i < n; ++i) {
        v_arr[i] = parts[i].v;
        pos_arr[i] = parts[i].pos;
        m_arr[i] = parts[i].mass;
    }

    vector_t *d_acc_arr = 0;
    vector_t *d_v_arr = 0;
    vector_t *d_pos_arr = 0;
    float *d_m_arr = 0;
    cudaMalloc((void**)&d_acc_arr, sizeof(vector_t) * n);
    cudaMalloc((void**)&d_v_arr, sizeof(vector_t) * n);
    cudaMalloc((void**)&d_pos_arr, sizeof(vector_t) * n);
    cudaMalloc((void**)&d_m_arr, sizeof(float) * n);
    
    cudaMemcpy(d_v_arr, v_arr, sizeof(vector_t) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pos_arr, pos_arr, sizeof(vector_t) * n, cudaMemcpyHostToDevice);
    cudaMemcpy(d_m_arr, m_arr, sizeof(float) * n, cudaMemcpyHostToDevice);

    const int threadsPerBlock = 256;
    int threadsPerGrid = (n + threadsPerBlock - 1) / threadsPerBlock;
    

    for (int m_i = 0; m_i < m; ++m_i) {

        compute_acceleration_local<<<threadsPerGrid, threadsPerBlock>>>
                    (d_acc_arr, d_pos_arr, d_m_arr, n, n, 0, grav);

        update_velocity_position_local<<<threadsPerGrid, threadsPerBlock>>>
                    (d_acc_arr, d_v_arr, d_pos_arr, d_m_arr, n, n, 0, dt);

        // last step or full output mode
        if ((m_rank == ROOT_NODE) && (m_i == (m - 1) || full_output)) {
            cudaMemcpy(pos_arr, d_pos_arr, sizeof(vector_t) * n, cudaMemcpyDeviceToHost);
            output_particle_pos_v(n, pos_arr, fp);
        }
    }


    cudaFree(d_acc_arr);
    cudaFree(d_v_arr);
    cudaFree(d_pos_arr);
    cudaFree(d_m_arr);

    free(v_arr);
    free(pos_arr);
    free(m_arr);
}
