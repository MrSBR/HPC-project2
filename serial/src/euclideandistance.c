#include "math.h"

double euclidean_distance(int kmax, int jmax, int imax, double ***arr1, double ***arr2) {
    double sum_sq_diff = 0.0;
    int k, j, i;

    // Iterate only over interior points
    for (k = 0; k < kmax; k++) {
        for (j = 0; j < jmax; j++) {
            for (i = 0; i < imax; i++) {
                double diff = arr1[k][j][i] - arr2[k][j][i];
                sum_sq_diff += diff * diff;
            }
        }
    }

    return sqrt(sum_sq_diff);
}
