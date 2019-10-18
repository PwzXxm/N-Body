#pragma once
#include "stdlib.h"
#include "float.h"
#include "stdio.h"
#include "utils.hpp"
#include "math.h"

// #define DEBUG

static const int CHILDREN_CNT = 4;
static const float SCALE_FACTOR = 1.5f;
static const int DX[] = {0, 1, 0, 1};
static const int DY[] = {0, 0, -1, -1};
static const float THETA = 0.6f;

#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct qt_mass {
    vector_t pos;
    float mass;
} qt_mass_t;

typedef struct qt_node {
    qt_mass_t mass_info;

    vector_t min_pos;
    float width;

    particle_t *particle;

    // top-left, top-right, bot-left, bot-right
    struct qt_node **children;
} qt_node_t;

void qt_sim(int n_particle, int n_steps, float time_step, particle_t *particles, float grav, FILE *f_out, bool is_full_out);

void qt_init(qt_node_t *root);

void qt_insert(particle_t *particle, qt_node_t *node);

size_t qt_find_ind(particle_t *p, qt_node_t *node);

qt_node_t *qt_new_node(vector_t pos, float width);

float qt_find_boundary(int n_particle, particle_t *particles);

void qt_print_tree(qt_node_t *root, int level);

qt_mass_t qt_compute_mass(qt_node_t *root);

vector_t qt_compute_force(particle_t *particle, qt_node_t *root, float grav);

float qt_dist(vector_t a, vector_t b);

void qt_free_tree(qt_node_t *root);

bool qt_is_out_of_boundary(particle_t p, float boundary);
