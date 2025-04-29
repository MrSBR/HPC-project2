// GSparallel.h
#ifndef GSPARALLEL_H // Use a standard include guard name
#define GSPARALLEL_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> // For memcpy, if used
#include <math.h>   // For sqrt
#include <mpi.h>

// --- Function Declarations (Prototypes) ---

// Function defined in allocatearray3D.c
void allocate_array3D(int kmax, int jmax, int imax, double ****array);

// Function defined in euclideandistance.c
double euclidean_distance(int kmax, int jmax, int imax, double ***arr1, double ***arr2);

// Function defined in GSiteration2chunks.c (the *serial* version for validation)
void GS_iteration_2_chunks(int kmax, int jmax, int imax, double ***phi);

// Function defined in gs_iteration_2_chunks_mpi.c (or wherever you implement it)
// Make sure the implementation file exists!
void GS_iteration_2_chunks_mpi(int my_rank, int kmax, int my_jmax, int imax, double ***my_phi);

void GS_iteration_2_chunks_mpi_test(int my_rank, int kmax, int my_jmax, int imax, double ***my_phi);

void handle_input(int nargs, char **args, int *num_iters, int *kmax, int *jmax, int *imax, bool *print_verbose);

#endif // GSPARALLEL_H