#include "seq_quad_tree.h"

void qt_sim(int n_particle, int n_steps, float time_step, particle_t *particles, FILE *f_out, bool is_full_out) {
    float boundary = qt_find_boundary(n_particle, particles);

    // TODO: add min x, y, width
    qt_node_t *root = qt_new_node((vector_t){.x = -boundary, .y = -boundary}, boundary*2);

    for (int i = 0; i < n_particle; i++) {
        qt_insert(&particles[i], root);
    }
}

void qt_init(qt_node_t *root) {}

void qt_insert(particle_t *p, qt_node_t *node) {
    if (node->particle == NULL) {
        if (node->children == NULL) {
            // leaf node with no particle
            // add particle in this leaf
            node->particle = p;
        } else {
            // intermediate node, traverse down
            qt_insert(p, node->children[qt_find_ind(p, node)]);
        }
    } else {
        // already has a particle in the node, move the particle into children
        node->children = (qt_node_t **)malloc(sizeof(qt_node_t *) * CHILDREN_CNT);

        float w_2 = node->width / 2;

        for (int i = 0; i < CHILDREN_CNT; i++) {
            node->children[i] = qt_new_node(
                (vector_t){
                    .x = node->min_pos.x + dx[i] * w_2,
                    .y = node->min_pos.y + dy[i] * w_2,
                },
                w_2);
        }

        qt_insert(node->particle, node->children[qt_find_ind(node->particle, node)]);
        node->particle = NULL;

        qt_insert(p, node->children[qt_find_ind(p, node)]);
    }
}

size_t qt_find_ind(particle_t *p, qt_node_t *node) {
    if (node->children == NULL) {
        fprintf(stderr, "The node does not have any child\n");
        exit(1);
    }

    float w_2 = node->width / 2;
    float mid_x = node->min_pos.x + w_2;
    float mid_y = node->min_pos.y + w_2;

    if (p->pos.y < mid_y) {
        if (p->pos.x < mid_x) {
            // top left
            return 0;
        } else {
            // top right
            return 1;
        }
    } else {
        if (p->pos.x < mid_x) {
            // bot left
            return 2;
        } else {
            // bot right
            return 3;
        }
    }

    if (node->children == NULL) {
        fprintf(stderr, "The node does not belong to any of child\n");
        exit(1);
    }
}

qt_node_t *qt_new_node(vector_t pos, float w) {
    qt_node_t *node = (qt_node_t *)malloc(sizeof(qt_node_t));

    node->min_pos = pos;
    node->width = w;

    node->particle = NULL;
    node->children = NULL;

    return node;
}

float qt_find_boundary(int n_particle, particle_t *particles) {
    float maxi = FLT_MIN;

    for (int i = 0; i < n_particle; i++) {
        float x = fabs(particles[i].pos.x);
        float y = fabs(particles[i].pos.y);
        maxi = MAX(maxi, MAX(x, y));
    }

    return maxi*SCALE_FACTOR;
}