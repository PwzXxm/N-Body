#include <stdio.h>
#include "utils.h"

__global__ void compute_acceleration(vector_t *d_acc_arr, vector_t *d_v_arr, 
                        vector_t *d_pos_arr, float *d_m_arr, int n, float grav) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n) return;
    d_acc_arr[i].x = 0;
    d_acc_arr[i].y = 0;
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        // compute force
        float dis2 = powf(d_pos_arr[i].x - d_pos_arr[j].x, 2) + powf(d_pos_arr[i].y - d_pos_arr[j].y, 2);
        float dis = sqrtf(dis2);
        float fx = 0, fy = 0;
        if (dis != 0.0f) {
            float f = grav * d_m_arr[i] * d_m_arr[j] / dis2;
            fx = f * (d_pos_arr[j].x - d_pos_arr[i].x) / dis;
            fy = f * (d_pos_arr[j].y - d_pos_arr[i].y) / dis;
        }
        // update acceleration
        d_acc_arr[i].x = fx / d_m_arr[i];
        d_acc_arr[i].y = fy / d_m_arr[i];
    }
}

__global__ void update_velocity_position(vector_t *d_acc_arr, vector_t *d_v_arr, 
                    vector_t *d_pos_arr, float *d_m_arr, int n, float dt) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= n) return;

    d_v_arr[i].x += d_acc_arr[i].x * dt;
    d_v_arr[i].y += d_acc_arr[i].y * dt;
    d_pos_arr[i].x += d_v_arr[i].x * dt;
    d_pos_arr[i].y += d_v_arr[i].y * dt;
}


extern "C" void __nbody_cuda_single_naive(int n, int m, float dt, particle_t parts[], float grav, FILE* fp, bool full_output) {
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

        compute_acceleration<<<threadsPerGrid, threadsPerBlock>>>
                    (d_acc_arr, d_v_arr, d_pos_arr, d_m_arr, n, grav);

        update_velocity_position<<<threadsPerGrid, threadsPerBlock>>>
                    (d_acc_arr, d_v_arr, d_pos_arr, d_m_arr, n, dt);

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
