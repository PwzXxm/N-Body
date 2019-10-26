/*----------------------------------------------------------------------------
 | Authors:
 |     - Weizhi Xu   (weizhix)  752454
 |     - Zijun Chen  (zijunc3)  813190
 -----------------------------------------------------------------------------*/
#pragma once

#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "utils.hpp"

// #define QT_SEQ_DEBUG

static const int QT_CHILDREN_CNT = 4;
static const float QT_SCALE_FACTOR = 1.5f;
static const float QT_THETA = 0.6f;
static const int INIT_CAPACITY = 1024;

static const int QT_DX[] = {0, 1, 0, 1};
static const int QT_DY[] = {0, 0, -1, -1};

#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*
 * The mass information about a node representing a cluster of particles
 */
typedef struct qt_mass {
    vector_t pos;
    float mass;
} qt_mass_t;

/*
 * The data structure for a BH node.
 */
typedef struct qt_node {
    qt_mass_t mass_info; // informatino about its mass and central of mass position

    vector_t min_pos; // top-left of the area
    vector_t len; // lenfth(width, height)

    int particle_idx; // the index of the particle in the node, if absent, default to -1

    // top-left, top-right, bot-left, bot-right
    int child_idx[QT_CHILDREN_CNT];
} qt_node_t;

extern MPI_Datatype mpi_qt_mass_t;
extern MPI_Datatype mpi_qt_node_t;

void init_MPI_datatype_qt();
void free_MPI_datatype_qt();

// typedef struct {
//     qt_node_t *arr;
//     size_t size;
//     size_t cap;
// } qt_array_t;
typedef std::vector<qt_node_t> tree_vec_t;

/*
 * The main entrance of the algorithm, the simualtion over the timestamps
 */
void qt_sim(int n_particle, int n_steps, float time_step, particle_t *particles, float grav, FILE *f_out, bool is_full_out);

/*
 * Insert a particle into the BH tree
 */
void qt_insert(particle_t *ps, int p_idx, tree_vec_t &tree_vec, int node_idx);

/*
 * Find the index in the children node of a particle accroding to its position
 */
size_t qt_find_ind(particle_t *ps, int p_idx, tree_vec_t &tree_vec, int node_idx);

/*
 * Initialize a BH tree node
 */
qt_node_t *qt_new_node(vector_t pos, vector_t len);

/*
 * Find the boundary of the simulation, using a factor times the maximum abosoluate coordinates
 */
float qt_find_boundary(int n_particle, particle_t *particles);

/*
 * Print the whole BH tree
 */
void qt_print_tree(particle_t *ps, tree_vec_t &tree_vec, int node_idx, int level);

/*
 * Compute the mass for the particular particle in the BH tree
 */
qt_mass_t qt_compute_mass(particle_t *ps, tree_vec_t &tree_vec, int node_idx);

/*
 * Compute the force for the partcular particle in the BH tree
 */ 
vector_t qt_compute_force(particle_t *ps, int idx, tree_vec_t &tree_vec, int node_idx, float grav);

/*
 * Calculate the distance between two position
 */
float qt_dist(vector_t a, vector_t b);

// void qt_reset_tree(qt_array_t *bh_tree);

/* 
 * Check if an particle is out of boundary or not
 */
bool qt_is_out_of_boundary(particle_t *p, float boundary);

// void qt_array_init(qt_array_t *a, int init_cap);
/*
 * Initialize a BH tree node 
 */
int qt_vec_append(tree_vec_t &v, vector_t pos, vector_t len);
// void qt_array_reserve(qt_array_t *a, int size);
// void qt_array_free(qt_array_t *a);
