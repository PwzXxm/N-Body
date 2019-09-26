#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

void print_usage();

int main(int argc, char* argv[]) {
    if (argc != 6) {
        print_usage();
        return 1;
    }
    char *algo_name = argv[1];
    int m = atoi(argv[2]);
    float t = atof(argv[3]);
    bool full_output = atoi(argv[4]);
    char *output_file = argv[5];
    if (m <= 0 && t <= 0) {
        print_usage();
        return 1;
    }

    int n = -1;
    scanf("%d", &n);
    particle_t *parts = (particle_t *) malloc(sizeof(particle_t) * n);
    
    for (int i = 0; i < n; ++i) {
        scanf("%f %f %f %f %f", &parts[i].pos.x, &parts[i].pos.y, 
                    &parts[i].v.x, &parts[i].v.y, &parts[i].weight);
    }
    if (n <= 0) {
        printf("Wrong input.\n");
        return 1;
    }

    FILE *fp = init_output_file(output_file, n, (full_output) ? 1 : m, t);
    if (fp == NULL) {
        printf("Unable to open output file.\n");
        return 1;
    }

    output_particle_pos(n, parts, fp);



    fclose(fp);
    free(parts);
    return 0;
}



void print_usage() {
    printf("Wrong arguments\n");
    printf("Usage: nbody <algorithm name> <num of steps> <time of each step> <full output> <output file>\n");
    printf("Example: nbody seq_naive 100 0.01\n");
}
