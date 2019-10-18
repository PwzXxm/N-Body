#pragma once

void get_mpi_size_rank(int* size, int* rank);
void sync_particle_positions(void *buffer, int count, int root);
