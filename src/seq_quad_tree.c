#include "seq_quad_tree.h"

void qt_sim(int n_particle, int n_steps, float dt, particle_t *particles, float grav, FILE *f_out, bool is_full_out) {
    float boundary = qt_find_boundary(n_particle, particles);

    qt_node_t *root;
    vector_t *acc = (vector_t *)malloc(sizeof(vector_t) * n_particle);

    // tree construction
    for (int step = 0; step < n_steps; step++) {
        root = qt_new_node((vector_t){.x = -boundary, .y = boundary}, boundary * 2);

#ifdef DEBUG
        printf("step: %d\n", step);
#endif

        // initialization
        for (int i = 0; i < n_particle; i++) {
            acc[i].x = 0.0f;
            acc[i].y = 0.0f;
        }

#ifdef DEBUG
        printf("******************************\n");
        printf("* Construct Tree             *\n");
        printf("******************************\n");
#endif
        for (int i = 0; i < n_particle; i++) {
            qt_insert(&particles[i], root);

#ifdef DEBUG
            // printf("Iteration: %d\n", i);
            // qt_print_tree(root, 0);
#endif
        }

        // compute centre mass and total mass at each node
        qt_compute_mass(root);

#ifdef DEBUG
        printf("******************************\n");
        printf("* Compute mass               *\n");
        printf("******************************\n");
        // qt_print_tree(root, 0);
#endif

        for (int i = 0; i < n_particle; i++) {
            vector_t forces = qt_compute_force(&particles[i], root, grav);
            acc[i].x += forces.x / particles[i].mass;
            acc[i].y += forces.y / particles[i].mass;
        }

        for (int i = 0; i < n_particle; i++) {
            particles[i].v.x += acc[i].x * dt;
            particles[i].v.y += acc[i].y * dt;
            particles[i].pos.x += particles[i].v.x * dt;
            particles[i].pos.y += particles[i].v.y * dt;

            boundary = MAX(boundary, MAX(fabs(particles[i].pos.x), fabs(particles[i].pos.y)));
        }
        boundary += 1.0f;

#ifdef DEBUG
        printf("Particles:\n");
        for (int i = 0; i < n_particle; i++) {
            printf("%d: %f\n", i, particles[i].pos.x);
        }

        // qt_print_tree(root, 0);
        // break;
#endif

        qt_free_tree(root);

        if (step == (n_steps - 1) || is_full_out) {
            output_particle_pos(n_particle, particles, f_out);
        }
    }

    free(acc);
}

void qt_init(qt_node_t *root) {}

void qt_insert(particle_t *p, qt_node_t *node) {
    if (node == NULL) return;

#ifdef DEBUG
    printf("... trying to insert p at (%f, %f) into node with min_pos (%f, %f) ...\n",
           p->pos.x,
           p->pos.y,
           node->min_pos.x,
           node->min_pos.y);
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
                    .x = node->min_pos.x + DX[i] * w_2,
                    .y = node->min_pos.y + DY[i] * w_2,
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

    return maxi * SCALE_FACTOR + 1.0f;
}

void qt_print_tree(qt_node_t *root, int level) {
    const char *const BLOCK_DIR[] = {"tl", "tr", "bl", "br"};

    if (root != NULL) {
        printf("Node level: %d\n", level);
        printf("\tmin_pos x, y, width: %f, %f, %f\n", root->min_pos.x, root->min_pos.y, root->width);
        printf("\tmass: %f\n", root->mass_info.mass);
        printf("\tmass_pos x, y: %f, %f\n", root->mass_info.pos.x, root->mass_info.pos.y);

        if (root->particle != NULL) {
            printf("\tparticle x, y: %f, %f\n", root->particle->pos.x, root->particle->pos.y);
        }

        if (root->children != NULL) {
            printf("\tchildren:\n\t\t");
            for (int i = 0; i < CHILDREN_CNT; i++) {
                printf("%s: %d; ", BLOCK_DIR[i], (root->children[i] == NULL ? 0 : 1));
            }
            printf("\n");
            printf("==============================\n");

            for (int i = 0; i < CHILDREN_CNT; i++) {
                if (root->children[i] != NULL) {
                    printf("Go down to %s\n", BLOCK_DIR[i]);
                    qt_print_tree(root->children[i], level + 1);
                }
            }
        } else {
            printf("==============================\n");
        }
    }
}

qt_mass_t qt_compute_mass(qt_node_t *root) {
    float sum_mass = 0.0f;
    vector_t cm_pos = {.x = 0.0f, .y = 0.0f};

    if (root->particle != NULL) {
        // contains one particle
        sum_mass = root->particle->mass;
        cm_pos = root->particle->pos;
    } else {
        if (root->children != NULL) {
            // intermediate node
            qt_mass_t *masses = (qt_mass_t *)malloc(sizeof(qt_mass_t) * CHILDREN_CNT);

            for (int i = 0; i < CHILDREN_CNT; i++) {
                masses[i] = qt_compute_mass(root->children[i]);
                sum_mass += masses[i].mass;
                cm_pos.x += (masses[i].mass * masses[i].pos.x);
                cm_pos.y += (masses[i].mass * masses[i].pos.y);
            }

            cm_pos.x /= sum_mass;
            cm_pos.y /= sum_mass;

            free(masses);
        }
    }

    root->mass_info.mass = sum_mass;
    root->mass_info.pos = cm_pos;

    return (qt_mass_t){.mass = sum_mass, .pos = cm_pos};
}

vector_t qt_compute_force(particle_t *p, qt_node_t *root, float grav) {
    vector_t f = {.x = 0.0f, .y = 0.0f};

    if (root->particle != NULL) {
        // leaf node
#ifndef DEBUG
        return force_between_particle(p->pos, root->particle->pos, p->mass, root->particle->mass, grav);
#else
        vector_t temp = force_between_particle(p->pos, root->particle->pos, p->mass, root->particle->mass, grav);
        printf("Compute forces between ((pos), mass): ((%f, %f), %f) and ((%f, %f), %f); result: (%f, %f)\n",
               p->pos.x,
               p->pos.y,
               p->mass,
               root->particle->pos.x,
               root->particle->pos.y,
               root->particle->mass,
               temp.x,
               temp.y);
        return temp;
#endif
    } else {
        // intermediate node

        if (root->children != NULL) {
            if ((root->width / qt_dist(root->mass_info.pos, p->pos)) < THETA) {
                return force_between_particle(p->pos, root->mass_info.pos, p->mass, root->mass_info.mass, grav);
            } else {
                for (int i = 0; i < CHILDREN_CNT; i++) {
                    vector_t temp = qt_compute_force(p, root->children[i], grav);
                    f.x += temp.x;
                    f.y += temp.y;
                }
            }
        }
    }

    return f;
}

inline float qt_dist(vector_t a, vector_t b) { return sqrtf(powf(a.x - b.x, 2) + powf(a.y - b.y, 2)); }

void qt_free_tree(qt_node_t *root) {
    if (root == NULL) return;

    if (root->children != NULL) {
        for (int i = 0; i < CHILDREN_CNT; i++) {
            qt_free_tree(root->children[i]);
        }

        free(root->children);
    }

    free(root);
}

inline bool qt_is_out_of_boundary(particle_t p, float boundary) {
    return (p.pos.x < -boundary || p.pos.x > boundary || p.pos.y < -boundary || p.pos.y > boundary);
}