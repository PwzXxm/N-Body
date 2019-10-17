#pragma once
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>

#include "seq_quad_tree.h"
#include "util.h"

#define QT_SWAP_INT(a,b) {int t;t=(a);(a)=(b);(b)=t;}

typedef struct qt_ORB_node {
    vector_t min_pos;
    vector_t len;

    bool end_node;
    int work_rank;
    int l, r;

    struct qt_ORB_node *left, *right;
} qt_ORB_node_t;

void qt_p_sim(int n_particle, int n_steps, float time_step, particle_t *ps, float grav, FILE *f_out, bool is_full_out);

void qt_ORB_with_level(qt_ORB_node_t *node, particle_t *ps, int *idx, int l, int r, int d, int lvl, int *w_rank, int size);

float qt_quick_select(particle_t *ps, int *idx, int l, int r, size_t offset, int k);

float qt_ps_idx_to_xy(particle_t *ps, int i, size_t offset);

bool qt_offset_gt(particle_t *ps, int *idx, int l, int r, size_t offset);

qt_ORB_node_t *qt_new_ORB_node(float x, float y, float x_len, float y_len);

void qt_free_ORB_tree(qt_ORB_node_t *root);

void qt_print_ORB_tree(qt_ORB_node_t *root, int d);

void qt_test_find_medium(particle_t *ps, int n);