#include "../include/GSMPI.h"

void GS_iteration_2_chunks_mpi(int my_rank, int kmax, int my_jmax, int imax, double ***my_phi) {
    int k, j, i;
    int count = imax - 2;
    int neighbor_rank = 1 - my_rank;
    MPI_Status status;

    double *send_col = (double *)malloc(count * sizeof(double));
    double *recv_col = (double *)malloc(count * sizeof(double));
    if (!send_col || !recv_col) {
        perror("Failed to allocate communication buffers");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // First Wavefront (k=1)
    if (my_rank == 0) {
        k = 1;
        for (j = 1; j < my_jmax - 1; j++) {
            for (i = 1; i < imax - 1; i++) {
                my_phi[k][j][i] = (my_phi[k - 1][j][i] + my_phi[k][j - 1][i] +
                                   my_phi[k][j][i - 1] + my_phi[k][j][i + 1] +
                                   my_phi[k][j + 1][i] + my_phi[k + 1][j][i]) / 6.0;
            }
        }
        // Send boundary to Rank 1
        for (i = 0; i < count; ++i) send_col[i] = my_phi[k][my_jmax - 2][i + 1];
        MPI_Send(send_col, count, MPI_DOUBLE, neighbor_rank, 10, MPI_COMM_WORLD);
        // Receive from Rank 1 (initial boundary, possibly not needed if k=1 not used)
        MPI_Recv(recv_col, count, MPI_DOUBLE, neighbor_rank, 11, MPI_COMM_WORLD, &status);
        for (i = 0; i < count; ++i) recv_col[i] = my_phi[k][my_jmax - 1][i + 1];
    } else {
        k = 1;
        // Send initial boundary (not computed, use initialized value)
        for (i = 0; i < count; ++i) send_col[i] = my_phi[k][1][i + 1];
        MPI_Send(send_col, count, MPI_DOUBLE, neighbor_rank, 11, MPI_COMM_WORLD);
        // Receive from Rank 0
        MPI_Recv(recv_col, count, MPI_DOUBLE, neighbor_rank, 10, MPI_COMM_WORLD, &status);
        for (i = 0; i < count; ++i) recv_col[i] = my_phi[k][0][i + 1];
    }

    // Intermediate Wavefronts
    for (k = 2; k <= kmax - 2; k++) {
        int k_right = k - 1;
        if (my_rank == 0) {
            for (j = 1; j < my_jmax - 1; j++) {
                for (i = 1; i < imax - 1; i++) {
                    my_phi[k][j][i] = (my_phi[k - 1][j][i] + my_phi[k][j - 1][i] +
                                       my_phi[k][j][i - 1] + my_phi[k][j][i + 1] +
                                       my_phi[k][j + 1][i] + my_phi[k + 1][j][i]) / 6.0;
                }
            }
        } else {
            for (j = 1; j < my_jmax - 1; j++) {
                for (i = 1; i < imax - 1; i++) {
                    my_phi[k_right][j][i] = (my_phi[k_right - 1][j][i] + my_phi[k_right][j - 1][i] +
                                             my_phi[k_right][j][i - 1] + my_phi[k_right][j][i + 1] +
                                             my_phi[k_right][j + 1][i] + my_phi[k_right + 1][j][i]) / 6.0;
                }
            }
        }
        // Communication
        int tag_k0 = k * 10 + 0;
        int tag_k1 = k * 10 + 1;
        if (my_rank == 0) {
            for (i = 0; i < count; ++i) send_col[i] = my_phi[k][my_jmax - 2][i + 1];
            MPI_Sendrecv(send_col, count, MPI_DOUBLE, neighbor_rank, tag_k0,
                         recv_col, count, MPI_DOUBLE, neighbor_rank, tag_k1,
                         MPI_COMM_WORLD, &status);
            for (i = 0; i < count; ++i) my_phi[k_right][my_jmax - 1][i + 1] = recv_col[i];
        } else {
            for (i = 0; i < count; ++i) send_col[i] = my_phi[k_right][1][i + 1];
            MPI_Sendrecv(send_col, count, MPI_DOUBLE, neighbor_rank, tag_k1,
                         recv_col, count, MPI_DOUBLE, neighbor_rank, tag_k0,
                         MPI_COMM_WORLD, &status);
            for (i = 0; i < count; ++i) my_phi[k][0][i + 1] = recv_col[i];
        }
    }

    // Last Wavefront
    if (my_rank == 1) {
        k = kmax - 2;
        for (j = 1; j < my_jmax - 1; j++) {
            for (i = 1; i < imax - 1; i++) {
                my_phi[k][j][i] = (my_phi[k - 1][j][i] + my_phi[k][j - 1][i] +
                                   my_phi[k][j][i - 1] + my_phi[k][j][i + 1] +
                                   my_phi[k][j + 1][i] + my_phi[k + 1][j][i]) / 6.0;
            }
        }
        // Send final boundary
        for (i = 0; i < count; ++i) send_col[i] = my_phi[k][1][i + 1];
        MPI_Send(send_col, count, MPI_DOUBLE, neighbor_rank, k * 10 + 1, MPI_COMM_WORLD);
    }
    if (my_rank == 0) {
        k = kmax - 2;
        for (i = 0; i < count; ++i) send_col[i] = my_phi[k][my_jmax - 2][i + 1];
        MPI_Recv(recv_col, count, MPI_DOUBLE, neighbor_rank, k * 10 + 1, MPI_COMM_WORLD, &status);
        for (i = 0; i < count; ++i) my_phi[k][my_jmax - 1][i + 1] = recv_col[i];
        MPI_Send(send_col, count, MPI_DOUBLE, neighbor_rank, k * 10 + 0, MPI_COMM_WORLD);
    } else {
        k = kmax - 2;
        MPI_Recv(recv_col, count, MPI_DOUBLE, neighbor_rank, k * 10 + 0, MPI_COMM_WORLD, &status);
        for (i = 0; i < count; ++i) my_phi[k][0][i + 1] = recv_col[i];
    }

    free(send_col);
    free(recv_col);
}