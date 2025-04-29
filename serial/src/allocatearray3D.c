#include "../include/GSserial.h"

void allocate_array3D(int kmax, int jmax, int imax, double ****array) {
    int k, j;

    // Allocate the first dimension (array of pointers to 2D arrays)
    *array = (double ***)malloc(kmax * sizeof(double **));
    if (*array == NULL) {
        perror("Failed to allocate memory for the first dimension");
        exit(EXIT_FAILURE);
    }

    // Allocate the second dimension (array of pointers to 1D arrays) for each k
    for (k = 0; k < kmax; k++) {
        (*array)[k] = (double **)malloc(jmax * sizeof(double *));
        if ((*array)[k] == NULL) {
            perror("Failed to allocate memory for the second dimension");
            // Clean up previously allocated memory before exiting
            for (int k_cleanup = 0; k_cleanup < k; k_cleanup++) {
                free((*array)[k_cleanup]);
            }
            free(*array);
            exit(EXIT_FAILURE);
        }

        // Allocate the third dimension (actual data array) for each k, j
        for (j = 0; j < jmax; j++) {
            (*array)[k][j] = (double *)malloc(imax * sizeof(double));
            if ((*array)[k][j] == NULL) {
                perror("Failed to allocate memory for the third dimension");
                // Clean up previously allocated memory before exiting
                 for (int k_cleanup = 0; k_cleanup <= k; k_cleanup++) {
                    for (int j_cleanup = 0; j_cleanup < (k_cleanup == k ? j : jmax) ; j_cleanup++) {
                         free((*array)[k_cleanup][j_cleanup]);
                    }
                    free((*array)[k_cleanup]);
                }
                free(*array);
                exit(EXIT_FAILURE);
            }
        }
    }
}


