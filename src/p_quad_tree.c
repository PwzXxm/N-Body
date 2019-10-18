#include "p_quad_tree.h"

void qt_p_sim(int n_particle, int n_steps, float time_step, particle_t *ps, float grav, FILE *f_out, bool is_full_out) {
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    // the number of particles is small, switch to sequential
    // if (log2(n_particle) < m_size) {
    //     qt_sim(n_particle, n_steps, time_step, ps, grav, f_out, is_full_out);
    //     return ;
    // }

    float boundary = qt_find_boundary(n_particle, ps);

    /* balance load using ORB */
    uint64_t start;
    start = GetTimeStamp();
    int work_rank_assign = 0;

    // int orb_lvl = 0;
    int orb_lvl = (int)pow(2, ceil(log2(m_size)))-1;
    if (m_rank == ROOT_NODE) {
        printf("Using ORB level: %d\n", orb_lvl);
    }

    // init an array of integers representing the index of the particle in 'ps'
    // as the order in the original particle array needs to be preserved
    int *ps_idx = (int *)malloc(sizeof(int) * n_particle);
    if (ps_idx == NULL) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < n_particle; i++) {
        ps_idx[i] = i;
        printf("%d: \t%f\t%f\n", i, ps[i].pos.x, ps[i].pos.y);
    }

    // for (int i = 0; i < n_particle; i++) {
    //     printf("%d: %f\n", ps_idx[i], ps[ps_idx[i]].pos.x);
    // }

    // qt_test_find_medium(ps, n_particle);

    qt_ORB_node_t *orb_root = qt_new_ORB_node(-boundary, boundary, boundary*2, boundary*2);
    orb_root->l = 0;
    orb_root->r = n_particle-1;
    qt_ORB_with_level(orb_root, ps, ps_idx, 0, n_particle-1, 0, orb_lvl, &work_rank_assign, m_size);

    qt_print_ORB_tree(orb_root, 0);

    printf("rank %d: Load balancing: %f sec\n", m_rank, GetTimeSpentInSeconds(start));

    /* construct quad tree (BH) on each node */
    start = GetTimeStamp();


    printf("rank %d: constructing tree: %f sec\n", m_rank, GetTimeSpentInSeconds(start));
    // send/recv tree

    // compute force

    if (m_rank == ROOT_NODE) {
        // if (step == (n_steps - 1) || is_full_out) {
            output_particle_pos(n_particle, ps, f_out);
        // }
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

    printf("median: %f with %d\n", median, k);

    d++;

    qt_ORB_node_t *new_node;
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
    return *(float *)(((void *)&(ps[i].pos))+offset);
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
    free (root);
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

    printf("offset of y: %d\n", offsetof(vector_t, y));
    vector_t v = {.x = 1.6f, .y = -2.4f};
    printf("====== %f\n", *(float *)(((void *)&(v))+4));


    float ans = qt_quick_select(ps, idx, 0, 4, offsetof(vector_t, y), (n&1)?(n+1)/2:n/2);

    printf("med: %f\n", ans);
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", idx[i], ps[idx[i]].pos.x);
    }
}