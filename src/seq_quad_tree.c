#include "seq_quad_tree.h"

void qt_sim(int n_particle, int n_steps, float time_step, particle_t *particles, float grav, FILE *f_out, bool is_full_out) {
    float boundary = qt_find_boundary(n_particle, particles);

    qt_node_t *root = qt_new_node((vector_t){.x = -boundary, .y = boundary}, boundary*2);

    for (int i = 0; i < n_particle; i++) {
        qt_insert(&particles[i], root);

#ifdef DEBUG
        printf("Iteration: %d\n", i);
        qt_print_tree(root, 0);
#endif
    }
}

void qt_init(qt_node_t *root) {}

void qt_insert(particle_t *p, qt_node_t *node) {

    if (node == NULL) return ;

#ifdef DEBUG
    printf("... trying to insert p at (%f, %f) into node with min_pos (%f, %f) ...\n", p->pos.x, p->pos.y, node->min_pos.x, node->min_pos.y);
#endif

    if (node->particle == NULL) {
        if (node->children == NULL) {
            // leaf node with no particle
            // add particle in this leaf
            node->particle = p;
        } else {
            // intermediate node, traverse down
            int index = qt_find_ind(p, node);
            qt_insert(p, index == CHILDREN_CNT ? NULL : node->children[index]);
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

        particle_t *tmp_particle = node->particle;
        node->particle = NULL;
        int index = qt_find_ind(tmp_particle, node);
        qt_insert(tmp_particle, index == CHILDREN_CNT ? NULL : node->children[index]);

        index = qt_find_ind(p, node);
        qt_insert(p, index == CHILDREN_CNT ? NULL : node->children[index]);
    }
}

size_t qt_find_ind(particle_t *p, qt_node_t *node) {
    if (node->children == NULL) {
        fprintf(stderr, "The node does not have any child\n");
        exit(1);
    }

    float w_2 = node->width / 2;
    float mid_x = node->min_pos.x + w_2;
    float mid_y = node->min_pos.y - w_2;

    if (p->pos.y > mid_y && p->pos.y <= node->min_pos.y) {
        if (p->pos.x < mid_x && p->pos.x >= node->min_pos.x) {
            // top left
            return 0;
        } else {
            // top right
            return 1;
        }
    } else {
        if (p->pos.x < mid_x && p->pos.x >= node->min_pos.x) {
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

void qt_print_tree(qt_node_t *root, int level) {
    const char* const BLOCK_DIR[] = {"tl", "tr", "bl", "br"};

    if (root != NULL) {
        printf("Node level: %d\n", level);
        printf("\tmin_pos x, y, width: %f, %f, %f\n", root->min_pos.x, root->min_pos.y, root->width);
        printf("\tmass: %f\n", root->mass);
        printf("\tmass_pos x, y: %f, %f\n", root->mass_pos.x, root->mass_pos.y);

        if (root->particle != NULL) {
            printf("\tparticle x, y: %f, %f\n", root->particle->pos.x, root->particle->pos.y);
        }

        if (root->children != NULL) {
            printf("\tchildren:\n\t\t");
            for (int i = 0; i < CHILDREN_CNT; i++) {
                printf("%s: %d; ", BLOCK_DIR[i], (root->children[i]==NULL ? 0 : 1));
            }
            printf("\n");
            printf("==============================\n");

            for (int i = 0; i < CHILDREN_CNT; i++) {
                if (root->children[i] != NULL) {
                    printf("Go down to %s\n", BLOCK_DIR[i]);
                    qt_print_tree(root->children[i], level+1);
                }
            }
        } else {
            printf("==============================\n");
        }
    }
}