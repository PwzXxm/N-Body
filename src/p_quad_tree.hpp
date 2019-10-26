/*----------------------------------------------------------------------------
 | Authors:
 |     - Weizhi Xu   (weizhix)  752454
 |     - Zijun Chen  (zijunc3)  813190
 -----------------------------------------------------------------------------*/
#pragma once
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "seq_quad_tree.hpp"
#include "utils.hpp"

#define QT_SWAP_INT(a,b) {int t;t=(a);(a)=(b);(b)=t;}

/*
 * The structure of Orthogonal Recursive Bisection (ORB) tree node
 */
typedef struct qt_ORB_node {
    vector_t min_pos; // top-left position that defines the area
    vector_t len; // length (width, height) of the area

    bool end_node; // if this is a leaf node
    int work_rank; // what node in MPI should work on this node
    int l, r; // the left and right indices in the array of indices of particles

    struct qt_ORB_node *left, *right; // the left and right children

    tree_vec_t *tree_vec; // if it is a leaf node, a BH tree is attached to the node
} qt_ORB_node_t;

/*
 * Simulate the process over timestamps, the main entrance of the algorithm
 */ 
void qt_p_sim(int n_particle, int n_steps, float time_step, particle_t *ps, float grav, FILE *f_out, bool is_full_out);

/*
 * Construt a ORB tree with a desired level as input
 * so the process stops at a particular level
 */
void qt_ORB_with_level(qt_ORB_node_t *node, particle_t *ps, int *idx, int l, int r, int d, int lvl, int *w_rank, int size);

/*
 * Use QuickSelect algorithm to find the median particle according to its x/y coordinate
 */ 
float qt_quick_select(particle_t *ps, int *idx, int n, bool is_horizon, int k);

/*
 * Get x/y coordinate from the array using index
 */
float qt_ps_idx_to_xy(particle_t *ps, int i, size_t offset);

/*
 * check if the particles position in x/y coordinate(depending on the offset)
 */
bool qt_offset_gt(particle_t *ps, int *idx, int l, int r, size_t offset);

/*
 * Initialize an ORB node
 */
qt_ORB_node_t *qt_new_ORB_node(float x, float y, float x_len, float y_len);

/*
 * Free an ORB tree
 */
void qt_free_ORB_tree(qt_ORB_node_t *root);

/*
 * print an ORB tree
 */
void qt_print_ORB_tree(qt_ORB_node_t *root, int d, particle_t *ps, int *idx);

// void qt_test_find_medium(particle_t *ps, int n);

/*
 * Construct a BH tree according to the particles in ORB nodes
 */
void qt_p_construct_BH(particle_t *ps, int *idx, qt_ORB_node_t *orb_root, int rank);


// void qt_serialize(qt_array_t *arr, qt_node_t *root);
// qt_node_t *qt_deserialize(qt_array_t *a);
// void qt_serialize_int(qt_array_t *a, int *x);
// void qt_serialize_float(qt_array_t *a, float *x);
// void qt_serialize_vector_t(qt_array_t *a, vector_t * v);


/*
 * Broadcast the particles in the ORB node to all other MPI nodes
 */
void qt_p_bcast(qt_ORB_node_t *node, int rank);

/*
 * Compute the forces for all particles in an ORB node
 */
void qt_p_compute_force(qt_ORB_node_t *orb_root, particle_t *ps, int *idx, float dt, float grav, int rank);

/*
 * gather all particles information to the MPI root node
 */
void qt_p_gather_particle(qt_ORB_node_t *node, particle_t *ps, int *idx, int rank);
