void GS_iteration_2_chunks (int kmax, int jmax, int imax, double ***phi)
{
    int k, j, i;
    int jmid = jmax / 2; // Midpoint for chunk division

    // --- First Wavefront ---
    // Only computation on the left chunk at level k=1
    k = 1;
    for (j = 1; j < jmid; j++) { // Left chunk: j=1 to jmid-1
        for (i = 1; i < imax - 1; i++) {
            phi[k][j][i] = (phi[k-1][j][i] + phi[k][j-1][i]
                          + phi[k][j][i-1] + phi[k][j][i+1]
                          + phi[k][j+1][i] + phi[k+1][j][i]) / 6.0;
        }
    }

    for (k = 2; k <= kmax - 2; k++) {
        // Computation on left chunk at level k
        for (j = 1; j < jmid; j++) { // Left chunk: j=1 to jmid-1
            for (i = 1; i < imax - 1; i++) {
                 phi[k][j][i] = (phi[k-1][j][i] + phi[k][j-1][i]
                               + phi[k][j][i-1] + phi[k][j][i+1]
                               + phi[k][j+1][i] + phi[k+1][j][i]) / 6.0;
            }
        }

        // Computation on right chunk at level k-1
        int k_right = k - 1; // Level for the right chunk
        for (j = jmid; j < jmax - 1; j++) { // Right chunk: j=jmid to jmax-2
            for (i = 1; i < imax - 1; i++) {
                 phi[k_right][j][i] = (phi[k_right-1][j][i] + phi[k_right][j-1][i]
                                     + phi[k_right][j][i-1] + phi[k_right][j][i+1]
                                     + phi[k_right][j+1][i] + phi[k_right+1][j][i]) / 6.0;
            }
        }
    }

    // --- Last Wavefront ---
    // Only computation on the right chunk at level k = kmax - 2
    k = kmax - 2;
    for (j = jmid; j < jmax - 1; j++) { // Right chunk: j=jmid to jmax-2
        for (i = 1; i < imax - 1; i++) {
            phi[k][j][i] = (phi[k-1][j][i] + phi[k][j-1][i]
                          + phi[k][j][i-1] + phi[k][j][i+1]
                          + phi[k][j+1][i] + phi[k+1][j][i]) / 6.0;
        }
    }
}