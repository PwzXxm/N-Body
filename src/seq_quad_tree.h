#pragma once
#include "stdlib.h"
#include "float.h"
#include "stdio.h"
#include "utils.h"
#include "math.h"

// #define DEBUG

// theta = width/distance

static const int CHILDREN_CNT = 4;
static const float SCALE_FACTOR = 1.5f;
static const int dx[] = {0, 1, 0, 1};
static const int dy[] = {0, 0, -1, -1};

#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct qt_node {
    vector_t mass_pos;
    float mass;

    vector_t min_pos;
    float width;

    particle_t *particle;

    // top-left, top-right, bot-left, bot-right
    struct qt_node **children;
} qt_node_t;

void qt_sim(int n_particle, int n_steps, float time_step, particle_t *particles, FILE *f_out, bool is_full_out);

void qt_init(qt_node_t *root);

void qt_insert(particle_t *particle, qt_node_t *node);

size_t qt_find_ind(particle_t *p, qt_node_t *node);

qt_node_t *qt_new_node(vector_t pos, float width);

float qt_find_boundary(int n_particle, particle_t *particles);

void qt_print_tree(qt_node_t *root, int level);


