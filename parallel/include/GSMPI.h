#ifndef GSPARALLEL_H 
#define GSPARALLEL_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> // For memcpy, if used
#include <math.h>   // For sqrt
#include <mpi.h>


// Function defined in allocatearray3D.c
void allocate_array3D(int kmax, int jmax, int imax, double ****array);

// Function defined in euclideandistance.c
double euclidean_distance(int kmax, int jmax, int imax, double ***arr1, double ***arr2);

// Function defined in GSiteration2chunksserial.c 
void GS_iteration_2_chunks(int kmax, int jmax, int imax, double ***phi);

// Function defined in GSiteration2chunksmpi.c
void GS_iteration_2_chunks_mpi(int my_rank, int kmax, int my_jmax, int imax, double ***my_phi);


#endif // GSPARALLEL_H