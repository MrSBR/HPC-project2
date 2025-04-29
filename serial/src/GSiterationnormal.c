#include <stdio.h>

void GS_iteration_normal (int kmax, int jmax, int imax, double ***phi)
{
    int k,j,i;
    // Iterate only over interior points
    for (k=1; k<kmax-1; k++) {
        for (j=1; j<jmax-1; j++) {
            for (i=1; i<imax-1; i++) {

                phi[k][j][i] = (phi[k-1][j][i] + phi[k][j-1][i]
                              + phi[k][j][i-1] + phi[k][j][i+1]
                              + phi[k][j+1][i] + phi[k+1][j][i]) / 6.0;
            }
        }
    }
}
