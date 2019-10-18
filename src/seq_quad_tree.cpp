#include "seq_quad_tree.hpp"

void qt_sim(int n_particle, int n_steps, float dt, particle_t *particles, float grav, FILE *f_out, bool is_full_out) {
    float boundary = qt_find_boundary(n_particle, particles);

    int root_idx;
    vector_t *acc = (vector_t *)malloc(sizeof(vector_t) * n_particle);
    bool *boundary_flags = (bool *)malloc(sizeof(bool) * n_particle);

    std::vector<qt_node_t> tree_vec;

    // tree construction
    for (int step = 0; step < n_steps; step++) {
        //root = qt_new_node((vector_t){.x = -boundary, .y = boundary}, (vector_t){.x = boundary*2, .y = boundary*2});
        root_idx = qt_vec_append(tree_vec, (vector_t){.x = -boundary, .y = boundary}, (vector_t){.x = boundary*2, .y = boundary*2});

#ifdef QT_SEQ_DEBUG
        printf("step: %d\n", step);
#endif
        
        // initialization
        for(int i = 0; i < n_particle; i++) {
            acc[i].x = 0.0f;
            acc[i].y = 0.0f;
            boundary_flags[i] = false;
        }

#ifdef QT_SEQ_DEBUG
        printf("******************************\n");
        printf("* Construct Tree             *\n");
        printf("******************************\n");
#endif
        for (int i = 0; i < n_particle; i++) {
            if (qt_is_out_of_boundary(&(particles[i]), boundary)) {
                boundary_flags[i] = true;
            }

            if (!boundary_flags[i]) {
                qt_insert(particles, i, tree_vec, root_idx);
            }

#ifdef QT_SEQ_DEBUG
            // printf("Iteration: %d\n", i);
            // qt_print_tree(root, 0);
#endif
        }

        // compute centre mass and total mass at each node
        qt_compute_mass(particles, tree_vec, root_idx);

#ifdef QT_SEQ_DEBUG
        printf("******************************\n");
        printf("* Compute mass               *\n");
        printf("******************************\n");
        // qt_print_tree(root, 0);
#endif

        for (int i = 0; i < n_particle; i++) {
            if (!boundary_flags[i]) {
                vector_t forces = qt_compute_force(particles, i, tree_vec, root_idx, grav);
                acc[i].x += forces.x / particles[i].mass;
                acc[i].y += forces.y / particles[i].mass;
            }
        }

        for (int i = 0; i < n_particle; i++) {
            if (!boundary_flags[i]) {
                particles[i].v.x += acc[i].x * dt;
                particles[i].v.y += acc[i].y * dt;
                particles[i].pos.x += particles[i].v.x * dt;
                particles[i].pos.y += particles[i].v.y * dt;
            }
        }



#ifdef QT_SEQ_DEBUG
        printf("Particles:\n");
        for (int i = 0; i < n_particle; i++) {
            printf("%d: %f\n", i, particles[i].pos.x);
        }

        // qt_print_tree(root, 0);
        // break;
#endif

        if (step == (n_steps - 1) || is_full_out) {
            output_particle_pos(n_particle, particles, f_out);
        }

        tree_vec.clear();
    }

    free(acc);
    free(boundary_flags);
}

void qt_init(qt_node_t *root) {}

void qt_insert(particle_t *ps, int p_idx, tree_vec_t &tree_vec, int node_idx) {
    if (node_idx < 0) return;

#ifdef QT_SEQ_DEBUG
    printf("... trying to insert p at (%f, %f) into node with min_pos (%f, %f) ...\n",
           p->pos.x,
           p->pos.y,
           tree_vec.at(node_idx).min_pos.x,
           tree_vec.at(node_idx).min_pos.y);
#endif
    tree_vec.at(node_idx);

    if (tree_vec.at(node_idx).particle_idx == -1) {
        if (tree_vec.at(node_idx).child_idx[0] == -1) {
            // leaf node with no particle
            // add particle in this leaf
            tree_vec.at(node_idx).particle_idx = p_idx;
        } else {
            // intermediate node, traverse down
            int index = qt_find_ind(ps, p_idx, tree_vec, node_idx);
            qt_insert(ps, p_idx, tree_vec, tree_vec.at(node_idx).child_idx[index]);
        }
    } else {
        // already has a particle in the node, move the particle into children
        float xl_2 = (tree_vec.at(node_idx).len.x) / 2;
        float yl_2 = (tree_vec.at(node_idx).len.y) / 2;

        for (int i = 0; i < QT_CHILDREN_CNT; i++) {
            tree_vec.at(node_idx).child_idx[i] = qt_vec_append(tree_vec,
                (vector_t){
                    .x = tree_vec.at(node_idx).min_pos.x + QT_DX[i] * xl_2,
                    .y = tree_vec.at(node_idx).min_pos.y + QT_DY[i] * yl_2,
                },
                (vector_t) {
                    .x = xl_2,
                    .y = yl_2,
                }
                );
        }

        int tmp_particle_idx = tree_vec.at(node_idx).particle_idx;
        tree_vec.at(node_idx).particle_idx = -1;
        int index = qt_find_ind(ps, tmp_particle_idx, tree_vec, node_idx);
        qt_insert(ps, tmp_particle_idx, tree_vec, tree_vec.at(node_idx).child_idx[index]);

        index = qt_find_ind(ps, p_idx, tree_vec, node_idx);
        qt_insert(ps, p_idx, tree_vec, tree_vec.at(node_idx).child_idx[index]);
    }
}

size_t qt_find_ind(particle_t *ps, int p_idx, tree_vec_t &tree_vec, int node_idx) {
    if (tree_vec.at(node_idx).child_idx[0] == -1) {
        fprintf(stderr, "The node does not have any child\n");
        exit(1);
    }

    float mid_x = tree_vec.at(node_idx).min_pos.x + ((tree_vec.at(node_idx).len.x)/2);
    float mid_y = tree_vec.at(node_idx).min_pos.y - ((tree_vec.at(node_idx).len.y)/2);

    if (ps[p_idx].pos.y > mid_y && ps[p_idx].pos.y <= tree_vec.at(node_idx).min_pos.y) {
        if (ps[p_idx].pos.x < mid_x && ps[p_idx].pos.x >= tree_vec.at(node_idx).min_pos.x) {
            // top left
            return 0;
        } else {
            // top right
            return 1;
        }
    } else {
        if (ps[p_idx].pos.x < mid_x && ps[p_idx].pos.x >= tree_vec.at(node_idx).min_pos.x) {
            // bot left
            return 2;
        } else {
            // bot right
            return 3;
        }
    }

    fprintf(stderr, "The node does not belong to any of child\n");
    exit(1);
}

float qt_find_boundary(int n_particle, particle_t *particles) {
    float maxi = FLT_MIN;

    for (int i = 0; i < n_particle; i++) {
        float x = fabs(particles[i].pos.x);
        float y = fabs(particles[i].pos.y);
        maxi = MAX(maxi, MAX(x, y));
    }

    return maxi * QT_SCALE_FACTOR;
}

void qt_print_tree(particle_t *ps, tree_vec_t &tree_vec, int node_idx, int level) {
    const char *const BLOCK_DIR[] = {"tl", "tr", "bl", "br"};

    printf("Node level: %d\n", level);
    printf("\tmin_pos x, y, x_len, y_len: %f, %f, %f, %f\n", tree_vec.at(node_idx).min_pos.x, tree_vec.at(node_idx).min_pos.y, tree_vec.at(node_idx).len.x, tree_vec.at(node_idx).len.y);
    printf("\tmass: %f\n", tree_vec.at(node_idx).mass_info.mass);
    printf("\tmass_pos x, y: %f, %f\n", tree_vec.at(node_idx).mass_info.pos.x, tree_vec.at(node_idx).mass_info.pos.y);

    if (tree_vec.at(node_idx).particle_idx != -1) {
        printf("\tparticle x, y: %f, %f\n", ps[tree_vec.at(node_idx).particle_idx].pos.x, ps[tree_vec.at(node_idx).particle_idx].pos.y);
    }

    if (tree_vec.at(node_idx).child_idx[0] == -1) {
        printf("\tchildren:\n\t\t");
        for (int i = 0; i < QT_CHILDREN_CNT; i++) {
            printf("%s: %d; ", BLOCK_DIR[i], (tree_vec.at(node_idx).child_idx[i] == -1 ? 0 : 1));
        }
        printf("\n");
        printf("==============================\n");

        for (int i = 0; i < QT_CHILDREN_CNT; i++) {
            if (tree_vec.at(node_idx).child_idx[i] != -1) {
                printf("Go down to %s\n", BLOCK_DIR[i]);
                qt_print_tree(ps, tree_vec, tree_vec.at(node_idx).child_idx[i], level + 1);
            }
        }
    } else {
        printf("==============================\n");
    }
}

qt_mass_t qt_compute_mass(particle_t *ps, tree_vec_t &tree_vec, int node_idx) {
    float sum_mass = 0.0f;
    vector_t cm_pos = {.x = 0.0f, .y = 0.0f};

    if (tree_vec.at(node_idx).particle_idx != -1) {
        // contains one particle
        sum_mass = ps[tree_vec.at(node_idx).particle_idx].mass;
        cm_pos = ps[tree_vec.at(node_idx).particle_idx].pos;
    } else {
        if (tree_vec.at(node_idx).child_idx[0] != -1) {
            // intermediate node
            qt_mass_t *masses = (qt_mass_t *)malloc(sizeof(qt_mass_t) * QT_CHILDREN_CNT);

            for (int i = 0; i < QT_CHILDREN_CNT; i++) {
                masses[i] = qt_compute_mass(ps, tree_vec, tree_vec.at(node_idx).child_idx[i]);
                sum_mass += masses[i].mass;
                cm_pos.x += (masses[i].mass * masses[i].pos.x);
                cm_pos.y += (masses[i].mass * masses[i].pos.y);
            }

            cm_pos.x /= sum_mass;
            cm_pos.y /= sum_mass;

            free(masses);
        }
    }

    tree_vec.at(node_idx).mass_info.mass = sum_mass;
    tree_vec.at(node_idx).mass_info.pos = cm_pos;

    return (qt_mass_t){.pos = cm_pos, .mass = sum_mass};
}

vector_t qt_compute_force(particle_t *ps, int p_idx, tree_vec_t &tree_vec, int node_idx, float grav) {
    vector_t f = {.x = 0.0f, .y = 0.0f};

    if (tree_vec.at(node_idx).particle_idx != -1) {
        // leaf node
#ifndef QT_SEQ_DEBUG
        return force_between_particle(ps[p_idx].pos, ps[tree_vec.at(node_idx).particle_idx].pos, ps[p_idx].mass, ps[tree_vec.at(node_idx).particle_idx].mass, grav);
#else
        vector_t temp = force_between_particle(p->pos, tree_vec.at(node_idx).particle->pos, p->mass, tree_vec.at(node_idx).particle->mass, grav);
        printf("Compute forces between ((pos), mass): ((%f, %f), %f) and ((%f, %f), %f); result: (%f, %f)\n",
            p->pos.x, p->pos.y, p->mass,
            tree_vec.at(node_idx).particle->pos.x, tree_vec.at(node_idx).particle->pos.y, tree_vec.at(node_idx).particle->mass,
            temp.x, temp.y
        );
        return temp;
#endif
    } else {
        // intermediate node

        if (tree_vec.at(node_idx).child_idx[0] != -1) {
            if ((((tree_vec.at(node_idx).len.x + tree_vec.at(node_idx).len.y)/2.0f) / qt_dist(tree_vec.at(node_idx).mass_info.pos, ps[p_idx].pos)) < QT_THETA) {
                return force_between_particle(ps[p_idx].pos, tree_vec.at(node_idx).mass_info.pos, ps[p_idx].mass, tree_vec.at(node_idx).mass_info.mass, grav);
            } else {
                for (int i = 0; i < QT_CHILDREN_CNT; i++) {
                    vector_t temp = qt_compute_force(ps, p_idx, tree_vec, tree_vec.at(node_idx).child_idx[i], grav);
                    f.x += temp.x;
                    f.y += temp.y;
                }
            }
        }
    }

    return f;
}

inline float qt_dist(vector_t a, vector_t b) { return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2)); }

// void qt_reset_tree(qt_array_t *bh_tree) {
//     bh_tree->size = 0;
// }

bool qt_is_out_of_boundary(particle_t *p, float boundary) {
    // printf("--- %f < %f\n", p->pos.x, -boundary);
    // printf("--- %f >= %f\n", p->pos.x, boundary);
    // printf("--- %f < %f\n", p->pos.y, -boundary);
    // printf("--- %f >= %f\n", p->pos.y, boundary);
    return (
        p->pos.x < -boundary ||
        p->pos.x >= boundary ||
        p->pos.y < -boundary ||
        p->pos.y >= boundary
    );
}

/******************************/
/*       Dynamic array        */
/******************************/
// void qt_array_init(qt_array_t *a, int init_cap) {
//     a->arr = (qt_node_t *)malloc(sizeof(qt_node_t) * init_cap);

//     if (a->arr == null) {
//         fprintf(stderr, "cannot allocate memory for bh_tree\n");
//         exit(1);
//     }

//     a->size = 0;
//     a->cap = init_cap;
// }

int qt_vec_append(std::vector<qt_node_t> &v, vector_t pos, vector_t len) {
    qt_node_t tmp;

    tmp.min_pos = pos;
    tmp.len = len;
    for (int i = 0; i < QT_CHILDREN_CNT; i++) {
        tmp.child_idx[i] = -1;
    }
    tmp.particle_idx = -1;

    tmp.mass_info.mass = 0.0f;
    tmp.mass_info.pos.x = 0.0f;
    tmp.mass_info.pos.y = 0.0f;

    v.push_back(tmp);

    // printf("size: %d, cap: %d\n", a->size, a->cap);
    // printf("%f \t%f\n", tmp.min_pos.x, pos.x);
    // printf("%f \t%f\n", tmp.min_pos.y, pos.y);
    // printf("%f \t%f\n", tmp.len.x, len.x);
    // printf("%f \t%f\n", tmp.len.y, len.y);

    return v.size()-1;
}

// void qt_array_reserve(qt_array_t *a, int size) {
//     while ((a->cap - a->size) < size) {
//         a->cap *= 2;
//         a->arr = (qt_node_t *)realloc(a->arr, sizeof(qt_node_t) * a->cap);
//     }
// }

// void qt_array_free(qt_array_t *a) {
//     free(a->arr);
//     a->arr = null;
// }
// /******************************/
