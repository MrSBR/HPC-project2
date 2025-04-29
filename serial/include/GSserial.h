// gs_serial.h
#ifndef GS_SERIAL_H
#define GS_SERIAL_H

#include <stdio.h> // For size_t, NULL, fprintf, stderr, printf
#include <stdlib.h> // For malloc, free, exit, atoi, srand, rand
#include <math.h>   // For sqrt
#include <time.h>   // For time()

// Allocates a 3D array dynamically.
// Note: Uses **** because the address of the 3D array pointer (***) is passed.
void allocate_array3D(int kmax, int jmax, int imax, double ****array);

// Performs one iteration of Gauss-Seidel using the "normal" traversal.
void GS_iteration_normal(int kmax, int jmax, int imax, double ***phi);

// Performs one iteration of Gauss-Seidel using the two-chunk wavefront traversal.
void GS_iteration_2_chunks(int kmax, int jmax, int imax, double ***phi);

// Calculates the Euclidean distance between the interior points of two 3D arrays.
double euclidean_distance(int kmax, int jmax, int imax, double ***arr1, double ***arr2);

#endif // GS_SERIAL_H