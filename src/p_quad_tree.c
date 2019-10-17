#include "p_quad_tree.h"

void qt_p_sim(int n_particle, int n_steps, float time_step, particle_t *ps, float grav, FILE *f_out, bool is_full_out) {
    int m_size, m_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

    float boundary = qt_find_boundary(n_particle, ps);

    int orb_lvl = (int)pow(2, ceil(log2(m_size)));

    if (m_rank == ROOT_NODE) {
        printf("Using ORB level: %d\n", orb_lvl);
    }

    uint64_t start, end;

    /* balance load using ORB */
    // start = GetTimeStamp();

    // init an array of integers representing the index of the particle in 'ps'
    // as the order in the original particle array needs to be preserved
    int *ps_idx = (int *)malloc(sizeof(int) * n_particle);
    if (ps_idx == NULL) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < n_particle; i++) {
        ps_idx[i] = i;
        printf("%d: %f\n", i, ps[i].pos.x);
    }
    // printf("rank %d: load balancing: %f\n", GetTimeSpentInSeconds(start));

    for (int i = 0; i < n_particle; i++) {
        printf("%d: %f\n", ps_idx[i], ps[ps_idx[i]].pos.x);
    }

    // qt_ORB_node_t *orb_root;
    // // qt_new_ORB_node(orb_root, boundary);
    // qt_ORB_with_level(orb_root, ps, ps_idx, 0, n_particle-1, orb_lvl);

    // construct quad tree - BH

    // send/recv tree

    // compute force

        // if (step == (n_steps - 1) || is_full_out) {
            output_particle_pos(n_particle, ps, f_out);
        // }

    free(ps_idx);
}

void qt_ORB_with_level(qt_ORB_node_t *node, particle_t *ps, int *idx, int l, int r, int lvl) {
    int n = r-l;
    int k = ((n & 1) ? (n+1)/2 : n/2);
    printf("median: %f with %d\n", qt_quick_select(ps, idx, l, r, offsetof(vector_t, x), k), k);
}

/*
 * find the k-th smallest elem in the array
 * offset is used to switch between x and y coordinate
 * 
 * quickselect algorithm reference: https://www.stat.cmu.edu/~ryantibs/median/quickselect.c
 */
float qt_quick_select(particle_t *ps, int *idx, int l, int r, size_t offset, int k) {
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
    return *(float *)(&(ps[i].pos)+offset);
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