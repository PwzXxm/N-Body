#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "utils.h"
// #include "seq_quad_tree.h"
#include "nbody_naive.h"

void print_usage(const char msg[]);
int run_task(int argc, char* argv[], int m_size, int m_rank);

static const char SEQ_QUAD_TREE[] = "seq_quad_tree";
static const char SEQ_NAIVE[] = "seq_naive";

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
    // parse arguments
    if (argc != 6) {
        if (m_rank == ROOT_NODE) print_usage("Wrong arguments\n");
        return 1;
    }

    char *algo_name = argv[1];
    int m = atoi(argv[2]);
    float t = atof(argv[3]);
    bool full_output = atoi(argv[4]);
    char *output_file = argv[5];
    if (m <= 0 || t <= 0) {
        if (m_rank == ROOT_NODE) print_usage("Invalid step or time interval\n");
        return 1;
    }

    if (strcmp(algo_name, SEQ_QUAD_TREE) == 0) {
        // algorithm_fun_ptr = &qt_sim;
    } else if (strcmp(algo_name, SEQ_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_seq_naive;
    } else {
        if (m_rank == ROOT_NODE) print_usage("Unsupported algorithm\n");
        return 1;
    }
    int n = -1;
    float grav;
    

    if (m_rank == ROOT_NODE) scanf("%d", &n);
    MPI_Bcast(&n, 1, MPI_INT, ROOT_NODE, MPI_COMM_WORLD);
    if (n <= 0) {
        if (m_rank == ROOT_NODE) printf("Wrong input.\n");
        return 1;
    }
    particle_t *parts = (particle_t *) malloc(sizeof(particle_t) * n);

    if (m_rank == ROOT_NODE) {        
        scanf("%f", &grav);

        for (int i = 0; i < n; ++i) {
            scanf("%f %f %f %f %f", &parts[i].pos.x, &parts[i].pos.y, 
                        &parts[i].v.x, &parts[i].v.y, &parts[i].mass);
        }
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

    (*algorithm_fun_ptr)(n, m, t, parts, grav, fp, full_output);

    if (m_rank == ROOT_NODE) fclose(fp);
    free(parts);

    return 0;
}


void print_usage(const char msg[]) {
    printf(msg);
    printf("Usage: nbody <algorithm name> <num of steps> <time of each step> <full output> <output file>\n");
    printf("Example: nbody seq_naive 100 0.01\n");
    printf("Supported algorithms:\n");
    printf("\tseq_quad_tree\n");
    printf("\tseq_naive\n");
}
