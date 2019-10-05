#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "utils.h"
#include "seq_quad_tree.h"
#include "nbody_naive.h"

#ifdef WITH_CUDA
#include "nbody_cuda.h"

static const char CUDA_SINGLE_NAIVE[] = "cuda_single_naive";
static const char CUDA_MPI_NAIVE[] = "cuda_mpi_naive";
#endif

void print_usage(const char msg[]);
int run_task(int argc, char* argv[], int m_size, int m_rank);

static const char SEQ_QUAD_TREE[] = "seq_quad_tree";
static const char SEQ_NAIVE[] = "seq_naive";
static const char MPI_OPENMP_NAIVE[] = "mpi_openmp_naive";

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    init_MPI_datatype();
    
    int r = run_task(argc, argv, m_size, m_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    free_MPI_datatype();
    MPI_Finalize();
    return r;
}

int run_task(int argc, char* argv[], int m_size, int m_rank) {
    #ifdef WITH_CUDA
    if (m_rank == ROOT_NODE) {
        printf("CUDA Enabled Version\n");
    }
    #endif

    // parse arguments
    // Note: There seems to be some bug when using MPI to pass large input though stdin,
    //       change to file instead.
    if (argc != 7) {
        if (m_rank == ROOT_NODE) print_usage("Wrong arguments\n");
        return 1;
    }

    char *algo_name = argv[1];
    int m = atoi(argv[2]);
    float t = atof(argv[3]);
    bool full_output = atoi(argv[4]);
    const char *input_file = argv[5];
    const char *output_file = argv[6];
    if (m <= 0 || t <= 0) {
        if (m_rank == ROOT_NODE) print_usage("Invalid step or time interval\n");
        return 1;
    }

    if (strcmp(algo_name, SEQ_QUAD_TREE) == 0) {
        algorithm_fun_ptr = &qt_sim;
    } else if (strcmp(algo_name, SEQ_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_seq_naive;
    } else if (strcmp(algo_name, MPI_OPENMP_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_mpi_openmp_naive;
    #ifdef WITH_CUDA
    } else if (strcmp(algo_name, CUDA_SINGLE_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_cuda_single_naive;
    } else if (strcmp(algo_name, CUDA_MPI_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_cuda_mpi_naive;
    #endif
    } else {
        if (m_rank == ROOT_NODE) print_usage("Unsupported algorithm\n");
        return 1;
    }
    int n = -1;
    float grav;

    FILE *input_fp = NULL;
    if (m_rank == ROOT_NODE) {
        input_fp = fopen(input_file, "r");
        if (input_file == NULL) {
            printf("Unable to open input file.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        fscanf(input_fp, "%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, ROOT_NODE, MPI_COMM_WORLD);
    if (n <= 0) {
        if (m_rank == ROOT_NODE) printf("Wrong input.\n");
        return 1;
    }
    particle_t *parts = (particle_t *) malloc(sizeof(particle_t) * n);

    if (m_rank == ROOT_NODE) {        
        fscanf(input_fp, "%f", &grav);

        for (int i = 0; i < n; ++i) {
            fscanf(input_fp, "%f %f %f %f %f", &parts[i].pos.x, &parts[i].pos.y, 
                        &parts[i].v.x, &parts[i].v.y, &parts[i].mass);
        }
        fclose(input_fp);
        input_fp = NULL;
    }

    MPI_Bcast(&grav, 1, MPI_FLOAT, ROOT_NODE, MPI_COMM_WORLD);
    MPI_Bcast(parts, n, mpi_particle_t, ROOT_NODE, MPI_COMM_WORLD);

    FILE *fp = NULL;
    if (m_rank == ROOT_NODE) {
        fp = init_output_file(output_file, n, (full_output) ? m : 1, t);
        if (fp == NULL) {
            printf("Unable to open output file.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
    }

    clock_t begin = clock();
    (*algorithm_fun_ptr)(n, m, t, parts, grav, fp, full_output);
    
    if (m_rank == ROOT_NODE) {
        double time_spent = (double)(clock() - begin) / CLOCKS_PER_SEC;
        printf("Time used: %f sec\n", time_spent);
    }


    if (m_rank == ROOT_NODE) fclose(fp);
    free(parts);

    return 0;
}


void print_usage(const char msg[]) {
    puts(msg);
    printf("Usage: nbody <algorithm name> <num of steps> <time of each step> <full output> <output file>\n");
    printf("Example: nbody seq_naive 100 0.01\n");
    printf("Supported algorithms:\n");
    printf("\tseq_quad_tree\n");
    printf("\tseq_naive\n");
    printf("\tmpi_openmp_naive\n");
    printf("\tcuda_single_naive\n");
    printf("\tcuda_mpi_naive\n");
}
