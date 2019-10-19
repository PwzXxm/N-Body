#pragma once
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "seq_quad_tree.hpp"
#include "util.h"

#define QT_SWAP_INT(a,b) {int t;t=(a);(a)=(b);(b)=t;}

typedef struct qt_ORB_node {
    vector_t min_pos;
    vector_t len;

    bool end_node;
    int work_rank;
    int l, r;

    struct qt_ORB_node *left, *right;

    tree_vec_t *tree_vec;
} qt_ORB_node_t;

void qt_p_sim(int n_particle, int n_steps, float time_step, particle_t *ps, float grav, FILE *f_out, bool is_full_out);

void qt_ORB_with_level(qt_ORB_node_t *node, particle_t *ps, int *idx, int l, int r, int d, int lvl, int *w_rank, int size);

float qt_quick_select(particle_t *ps, int *idx, int n, bool is_horizon, int k);

float qt_ps_idx_to_xy(particle_t *ps, int i, size_t offset);

bool qt_offset_gt(particle_t *ps, int *idx, int l, int r, size_t offset);

qt_ORB_node_t *qt_new_ORB_node(float x, float y, float x_len, float y_len);

void qt_free_ORB_tree(qt_ORB_node_t *root);

void qt_print_ORB_tree(qt_ORB_node_t *root, int d, particle_t *ps, int *idx);

void qt_test_find_medium(particle_t *ps, int n);

void qt_p_construct_BH(particle_t *ps, int *idx, qt_ORB_node_t *orb_root, int rank);


// void qt_serialize(qt_array_t *arr, qt_node_t *root);
// qt_node_t *qt_deserialize(qt_array_t *a);
// void qt_serialize_int(qt_array_t *a, int *x);
// void qt_serialize_float(qt_array_t *a, float *x);
// void qt_serialize_vector_t(qt_array_t *a, vector_t * v);


void qt_p_bcast(qt_ORB_node_t *node, int rank);

void qt_p_compute_force(qt_ORB_node_t *orb_root, particle_t *ps, int *idx, float dt, float grav, int rank);

void qt_p_gather_particle(qt_ORB_node_t *node, particle_t *ps, int *idx, int rank);