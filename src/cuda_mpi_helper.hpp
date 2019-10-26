/*----------------------------------------------------------------------------
 | Authors:
 |     - Weizhi Xu   (weizhix)  752454
 |     - Zijun Chen  (zijunc3)  813190
 -----------------------------------------------------------------------------*/
#pragma once

void get_mpi_size_rank(int* size, int* rank);
void sync_particle_positions(void *buffer, int count, int root);
