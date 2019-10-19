#include "p_quad_tree.hpp"

void qt_p_sim(int n_particle, int n_steps, float dt, particle_t *ps, float grav, FILE *f_out, bool is_full_out) {
    int m_size, m_rank;
    uint64_t start;

    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    // the number of particles is small, switch to sequential
    // if (log2(n_particle) < m_size) {
    //     qt_sim(n_particle, n_steps, time_step, ps, grav, f_out, is_full_out);
    //     return ;
    // }

    const float boundary = qt_find_boundary(n_particle, ps);

    // init an array of integers representing the index of the particle in 'ps'
    // as the order in the original particle array needs to be preserved
    int *ps_idx = (int *)malloc(sizeof(int) * n_particle);
    if (ps_idx == NULL) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int step = 0; step < n_steps; step++) {

        /* balance load using ORB */
        // start = GetTimeStamp();
        int work_rank_assign = 0;

        int orb_lvl = (int)pow(2, ceil(log2(m_size)))-1;
        // if (m_rank == ROOT_NODE) {
        //     printf("Using ORB level: %d\n", orb_lvl);
        // }

        // counter for how many particles within the boundary
        int p_cnt = 0;
        for (int i = 0; i < n_particle; i++) {
            if (!qt_is_out_of_boundary(&(ps[i]), boundary)) {
                // printf("%d: \t%f\t%f\n", p_cnt, ps[i].pos.x, ps[i].pos.y);
                ps_idx[p_cnt++] = i;
            }
        }

        // for (int i = 0; i < n_particle; i++) {
        //     printf("%d: %f\n", ps_idx[i], ps[ps_idx[i]].pos.x);
        // }

        // qt_test_find_medium(ps, n_particle);

        qt_ORB_node_t *orb_root = qt_new_ORB_node(-boundary, boundary, boundary*2, boundary*2);
        orb_root->l = 0;
        orb_root->r = p_cnt-1;
        qt_ORB_with_level(orb_root, ps, ps_idx, 0, p_cnt-1, 0, orb_lvl, &work_rank_assign, m_size);

        // qt_print_ORB_tree(orb_root, 0);

        // printf("rank %d: Load balancing: %f sec\n", m_rank, GetTimeSpentInSeconds(start));

        /* construct quad tree (BH) on each node */
        // start = GetTimeStamp();
        qt_p_construct_BH(ps, ps_idx, orb_root, m_rank);
        // printf("rank %d: constructing tree: %f sec\n", m_rank, GetTimeSpentInSeconds(start));

        // send/recv tree
        qt_p_bcast(orb_root, m_rank);

        // compute force
        qt_p_compute_force(orb_root, ps, ps_idx, dt, grav, m_rank);

        if (m_rank == ROOT_NODE) {
            if (step == (n_steps - 1) || is_full_out) {
                output_particle_pos(n_particle, ps, f_out);
            }
        }
    }
    free(ps_idx);
}

void qt_ORB_with_level(qt_ORB_node_t *node, particle_t *ps, int *idx, int l, int r, int d, int lvl, int *w_r, int size) {
    // printf("lvl %d: %d, %d\n", d, l, r);

    if (d == lvl || r == l) {
        node->end_node = true;
        node->work_rank = ((*w_r)++ % size);
        return ;
    }

    int n = r-l+1;
    int k = ((n & 1) ? (n+1)/2 : n/2);

    bool is_horizon = d&1;
    size_t offset = is_horizon ? offsetof(vector_t, y) : offsetof(vector_t, x);
    // size_t offset = offsetof(vector_t, x);
    float median = qt_quick_select(ps, idx, l, r, offset, k);

    if (median == FLT_MAX) {
        fprintf(stderr, "Unable to find median\n");
        exit(1);
    }

    // printf("median: %f with %d\n", median, k);

    d++;

    qt_ORB_node_t *new_node;

    // TODO: openmp
    if (l <= l+k-1) {
        // divide the space horizontally, thus, the height, y_len shrinks
        if (is_horizon) {
            new_node = qt_new_ORB_node(node->min_pos.x, node->min_pos.y,
                                       node->len.x, (node->min_pos.y)-median);
        } else {
            new_node = qt_new_ORB_node(node->min_pos.x, node->min_pos.y,
                                       median-(node->min_pos.x), node->len.y);
        }
        new_node->l = l;
        new_node->r = l+k-1;
        node->left = new_node;
        qt_ORB_with_level(new_node, ps, idx, l, l+k-1, d, lvl, w_r, size);
    }
    if (l+k <= r) {
        if (is_horizon) {
            new_node = qt_new_ORB_node(node->min_pos.x, median,
                                       node->len.x, median - (node->min_pos.y - node->len.y));
        } else {
            new_node = qt_new_ORB_node(median, node->min_pos.y,
                                       (node->min_pos.x + node->len.x) - median, node->len.y);
        }
        new_node->l = l+k;
        new_node->r = r;
        node->right = new_node;
        qt_ORB_with_level(new_node, ps, idx, l+k, r, d, lvl, w_r, size);
    }
}

/*
 * find the k-th smallest elem in the array
 * offset is used to switch between x and y coordinate
 * 
 * quickselect algorithm reference: https://www.stat.cmu.edu/~ryantibs/median/quickselect.c
 */
float qt_quick_select(particle_t *ps, int *idx, int l, int r, size_t offset, int k) {
    // printf("Start finding median\n");
    // for (int i = 0; i < 10; i++) {
    //     printf("%d ", idx[i]);
    // }
    // printf("\n");
    // printf("l, r: %d, %d\n", l, r);
    // printf("offset: %d\n", offset);
    // printf("k: %d\n", k);
    if (k > (r-l)) return FLT_MAX;

    for (;;) {
        if (r <= l+1) {
            if (r == l+1 && qt_offset_gt(ps, idx, l, r, offset)) {
                QT_SWAP_INT(idx[l], idx[r]);
            }
            return qt_ps_idx_to_xy(ps, idx[k], offset);
        }

        // choose middle item as the pivot
        int mid = (l+r) >> 2, ll = l+1, rr=r;
        QT_SWAP_INT(idx[mid], idx[l+1]);
        
        // sort element in l, l+1, r in ascending order
        if (qt_offset_gt(ps, idx, l, r, offset)) { QT_SWAP_INT(idx[l], idx[r]); }
        if (qt_offset_gt(ps, idx, l+1, r, offset)) { QT_SWAP_INT(idx[l+1], idx[r]); }
        if (qt_offset_gt(ps, idx, l, l+1, offset)) { QT_SWAP_INT(idx[l], idx[l+1]); }

        // put all elems smaller than mid in the left, others in the right
        for (;;) {
            // these two for loops are independent thus can be parallelized
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    for (ll++;qt_offset_gt(ps, idx, l+1, ll, offset);ll++);
                }
                #pragma omp section
                {
                    for (rr--;qt_offset_gt(ps, idx, rr, l+1, offset);rr--);
                }
            }

            if (rr < ll) break;
            QT_SWAP_INT(idx[ll], idx[rr]);
        }
        // swap the pivot back
        QT_SWAP_INT(idx[l+1], idx[rr]);
        if (rr >= k) r = rr-1;
        if (rr <= k) l = ll;
    }
}

/*
 * Get x/y coordinate from the array using index
 */
inline float qt_ps_idx_to_xy(particle_t *ps, int i, size_t offset) {
    return *(float *)(((char *)&(ps[i].pos))+offset);
}

/*
 * check if the particles position in x/y coordinate(depending on the offset)
 */
bool qt_offset_gt(particle_t *ps, int *idx, int l, int r, size_t offset) {
    if (qt_ps_idx_to_xy(ps, idx[l], offset) > qt_ps_idx_to_xy(ps, idx[r], offset))  {
        return true;
    } else {
        return false;
    }
}

/*
 * Construct a new ORB node
 */
qt_ORB_node_t *qt_new_ORB_node(float x, float y, float x_len, float y_len) {
    qt_ORB_node_t *node = (qt_ORB_node_t *)malloc(sizeof(qt_ORB_node_t));

    node->min_pos.x = x;
    node->min_pos.y = y;
    node->len.x = x_len;
    node->len.y = y_len;

    node->end_node = false;

    node->work_rank = 0;

    node->l = node->r = 0;

    node->tree_vec = new tree_vec_t();

    node->left = NULL;
    node->right = NULL;

    return node;
}

/*
 * free the whole ORB tree
 */
void qt_free_ORB_tree(qt_ORB_node_t *root) {
    if (root == NULL) return;
    if (root->left != NULL) qt_free_ORB_tree(root->left);
    if (root->right != NULL) qt_free_ORB_tree(root->right);

    root->tree_vec->clear();
    if (root->tree_vec != NULL) delete root->tree_vec;
    root->tree_vec = NULL;

    free (root);
    root = NULL;
}

void qt_print_ORB_tree(qt_ORB_node_t *root, int d) {
    if (root == NULL) return ;

    printf("Level: %d:\nmin_pos: %f, %f\nlen:%f, %f\nend_node: %d\nl, r: %d, %d\nwork_rank: %d\n",
            d, root->min_pos.x, root->min_pos.y, root->len.x, root->len.y,
            root->end_node, root->l, root->r, root->work_rank);

    qt_print_ORB_tree(root->left, d+1);
    qt_print_ORB_tree(root->right, d+1);
}

void qt_test_find_medium(particle_t *ps, int n) {
    int idx[] = {2, 6, 1, 5, 0, 9, 3, 7, 8, 4};
    
    n = 5;

    printf("offset of y: %lu\n", offsetof(vector_t, y));
    vector_t v = {.x = 1.6f, .y = -2.4f};
    printf("====== %f\n", *(float *)(((char *)&(v))+4));


    float ans = qt_quick_select(ps, idx, 0, 4, offsetof(vector_t, y), (n&1)?(n+1)/2:n/2);

    printf("med: %f\n", ans);
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", idx[i], ps[idx[i]].pos.x);
    }
}

void qt_p_construct_BH(particle_t *ps, int *idx, qt_ORB_node_t *node, int rank) {
    if (node == NULL) return ;
    if (node->end_node && node->work_rank == rank) {
        // construct local BH

        // printf("rk: %d, pos: %f, %f, len: %f, %f, v_size: %d\n", rank, node->min_pos.x, node->min_pos.y, node->len.x, node->len.y, node->tree_vec->size());
        int root_node = qt_vec_append(*(node->tree_vec), node->min_pos, node->len);

        for (int i = node->l; i <= node->r; i++) {
            qt_insert(ps, idx[i], *(node->tree_vec), root_node);
        }

        // compute mass for the tree
        qt_compute_mass(ps, *(node->tree_vec), root_node);
    } else {
        qt_p_construct_BH(ps, idx, node->left, rank);
        qt_p_construct_BH(ps, idx, node->right, rank);
    }
}

// /******************************/
// /*       Serialization        */
// /******************************/
// void qt_serialize(qt_array_t *a, qt_node_t *root) {
//     int x = 111111;
//     qt_serialize_int(a, &x);

//     printf("size: %d, cap: %d\n", a->size, a->cap);
//     for (int i = 0; i < a->size; i++) {
//         printf("%X ", a->arr[i]);
//     }
//     printf("\n");
// }

// qt_node_t *qt_deserialize(qt_array_t *a) {
// }

// void qt_serialize_int(qt_array_t *a, int *x) {
//     qt_array_reserve(a, sizeof(*x));
//     memcpy(a->arr+a->size, x, sizeof(*x));
//     a->size += sizeof(*x);
// }


// void qt_serialize_float(qt_array_t *a, float *x) {
//     qt_array_reserve(a, sizeof(*x));
//     memcpy(a->arr+a->size, x, sizeof(*x));
//     a->size += sizeof(*x);
// }

// void qt_serialize_vector_t(qt_array_t *a, vector_t * v) {
// }
/******************************/

void qt_p_bcast(qt_ORB_node_t *node, int rank) {
    if (node == NULL) return;

    if (node->end_node) {
        int n = node->tree_vec->size();
        MPI_Bcast(&n, 1, MPI_INT, node->work_rank, MPI_COMM_WORLD);
        node->tree_vec->resize(n);
        MPI_Bcast(node->tree_vec->data(), n, mpi_qt_node_t, node->work_rank, MPI_COMM_WORLD);
    } else {
        qt_p_bcast(node->left, rank);
        qt_p_bcast(node->right, rank);
    }
}

void qt_p_compute_force(qt_ORB_node_t *node, particle_t *ps, int *idx, float dt, float grav, int rank) {
    if (node == NULL) return;
    if (node->end_node && node->work_rank == rank) {
        for (int i = node->l; i <= node->r; i++) {
            vector_t forces = qt_compute_force(ps, idx[i], *(node->tree_vec), 0, grav);
            float acc_x = 0.0f, acc_y = 0.0f;

            particle_t *p = &ps[idx[i]];

            acc_x += forces.x / p->mass;
            acc_y += forces.y / p->mass;

            p->v.x += acc_x * dt;
            p->v.y += acc_y * dt;

            p->pos.x += p->v.x * dt;
            p->pos.y += p->v.y * dt;
        }
    } else {
        qt_p_compute_force(node->right, ps, idx, dt, grav, rank);
        qt_p_compute_force(node->left, ps, idx, dt, grav, rank);
    }
}

void gather_particle(qt_ORB_node_t *node, particle_t *ps, int *idx, int rank) {
    if (node == NULL) return;
    if (node->end_node) {
        if (rank == ROOT_NODE) {
            if (node->work_rank != rank) {
                for (int i = node->l; i <= node->r; i++) {
                    MPI_Recv(&ps[idx[i]], 1, mpi_particle_t, node->work_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        } else {
            if (node->work_rank == rank) {
                for (int i = node->l; i <= node->r; i++) {
                    MPI_Send(&ps[idx[i]], 1, mpi_particle_t, ROOT_NODE, 0, MPI_COMM_WORLD);
                }
            }
        }
    } else {
        gather_particle(node->left, ps, idx, rank);
        gather_particle(node->right, ps, idx, rank);
    }
}