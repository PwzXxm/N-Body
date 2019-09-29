#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
// #include "seq_quad_tree.h"
#include "nbody_seq_naive.h"

void print_usage();

static const char SEQ_QUAD_TREE[] = "seq_quad_tree";
static const char SEQ_NAIVE[] = "seq_naive";

// Unit:
//      10^5 kg
//      10^10 m

int main(int argc, char* argv[]) {
    if (argc != 6) {
        printf("Wrong arguments\n");
        print_usage();
        return 1;
    }
    char *algo_name = argv[1];
    int m = atoi(argv[2]);
    float t = atof(argv[3]);
    bool full_output = atoi(argv[4]);
    char *output_file = argv[5];
    if (m <= 0 || t <= 0) {
        printf("Invalid step or time interval\n");
        print_usage();
        return 1;
    }

    int n = -1;
    float grav;
    scanf("%d", &n);
    scanf("%f", &grav);
    particle_t *parts = (particle_t *) malloc(sizeof(particle_t) * n);
    
    for (int i = 0; i < n; ++i) {
        scanf("%f %f %f %f %f", &parts[i].pos.x, &parts[i].pos.y, 
                    &parts[i].v.x, &parts[i].v.y, &parts[i].mass);
    }
    if (n <= 0) {
        printf("Wrong input.\n");
        return 1;
    }

    FILE *fp = init_output_file(output_file, n, (full_output) ? m : 1, t);
    if (fp == NULL) {
        printf("Unable to open output file.\n");
        return 1;
    }

    if (strcmp(algo_name, SEQ_QUAD_TREE) == 0) {
        // algorithm_fun_ptr = &qt_sim;
    } else if (strcmp(algo_name, SEQ_NAIVE) == 0) {
        algorithm_fun_ptr = &nbody_seq_naive;
    } else {
        printf("Unsupported algorithm\n");
        print_usage();
        return 1;
    }

    (*algorithm_fun_ptr)(n, m, t, parts, grav, fp, full_output);

    output_particle_pos(n, parts, fp);


    fclose(fp);
    free(parts);
    return 0;
}



void print_usage() {
    printf("Usage: nbody <algorithm name> <num of steps> <time of each step> <full output> <output file>\n");
    printf("Example: nbody seq_naive 100 0.01\n");
    printf("Supported algorithms:\n");
    printf("\tseq_quad_tree\n");
    printf("\tseq_naive\n");
}
