double euclidean_distance(int kmax, int jmax, int imax, double ***arr1, double ***arr2) {
    double sum_sq_diff = 0.0;
    int k, j, i;

    // Iterate only over interior points
    for (k = 1; k < kmax - 1; k++) {
        for (j = 1; j < jmax - 1; j++) {
            for (i = 1; i < imax - 1; i++) {
                double diff = arr1[k][j][i] - arr2[k][j][i];
                sum_sq_diff += diff * diff;
            }
        }
    }

    return sqrt(sum_sq_diff);
}