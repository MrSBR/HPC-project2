#include "../include/GSMPI.h"

void GS_iteration_2_chunks_mpi(int my_rank, int kmax, int my_jmax, int imax, double ***my_phi) {
    int count = imax; // Number of elements to send/recieve
    int tag;          // Links the correct send/recieve statements to eachother
    int destination;  // Destination rank
    int source;       // Source rank
    
    // First wavefront: Computation left chunk level k=1
    if (my_rank == 0) {
        destination = 1;
        source = 1;

        for (int j=1; j<my_jmax-1; j++) {
            for (int i=1; i<imax-1; i++) {
                my_phi[1][j][i] = (my_phi[0][j][i]   + my_phi[1][j-1][i]
                                  +my_phi[1][j][i-1] + my_phi[1][j][i+1]
                                  +my_phi[1][j+1][i] + my_phi[2][j][i])/6.0;
            }
        }
        // Sending the second to last row at k = 1 to process 1 
        tag = 3; // Tag for level k=1
        MPI_Send(my_phi[1][my_jmax-2], count, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
        
        // Computation left chunk level k
        for (int k=2; k<=kmax-2; k++) {
            // Receiving the first row at k - 1 from process 1
            tag = 4; // Tag for level k-1
            source = 1;
            MPI_Recv(my_phi[k-1][my_jmax-1], count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int j=1; j<my_jmax-1; j++) {
                for (int i=1; i<imax-1; i++) {
                    my_phi[k][j][i] = (my_phi[k-1][j][i] + my_phi[k][j-1][i]
                                      +my_phi[k][j][i-1] + my_phi[k][j][i+1]
                                      +my_phi[k][j+1][i] + my_phi[k+1][j][i])/6.0;
                }
            }

            tag = 5; // Tag for level k
            MPI_Send(my_phi[k][my_jmax-2], count, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
        }

        // Receiving the last row at k = kmax-2 from process 1
        tag = 6; // Tag for level k=kmax-2
        source = 1;
        MPI_Recv(my_phi[kmax-2][my_jmax-1], count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
        
    if (my_rank == 1) {
        destination = 0;
        source = 0;

        // Recieving the second to last row at k = 1 to process 1
        tag = 3; // Tag for level k=1
        MPI_Recv(my_phi[1][0], count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Rewriting the k-loop to seperate the two chunks
        // Computation right chunk level k-1
        for (int k=2; k<=kmax-2; k++) {
            for (int j=1; j<my_jmax-1; j++) {
                for (int i=1; i<imax-1; i++) {
                    my_phi[k-1][j][i] = (my_phi[k-2][j][i]   + my_phi[k-1][j-1][i]
                                        +my_phi[k-1][j][i-1] + my_phi[k-1][j][i+1]
                                        +my_phi[k-1][j+1][i] + my_phi[k][j][i])/6.0;
                    }                
                }
            // Sending the first row at k - 1 to process 0
            tag = 4; // Tag for level k-1
            MPI_Send(my_phi[k-1][1], count, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);

            tag = 5; // Tag for level k
            MPI_Recv(my_phi[k][0], count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Last wavefront: Computation right chunk level k=kmax-2
        for (int j=1; j<my_jmax-1; j++) {
            for (int i=1; i<imax-1; i++) {
                my_phi[kmax-2][j][i] = (my_phi[kmax-3][j][i]   + my_phi[kmax-2][j-1][i]
                                       +my_phi[kmax-2][j][i-1] + my_phi[kmax-2][j][i+1]
                                       +my_phi[kmax-2][j+1][i] + my_phi[kmax-1][j][i])/6.0;
                }
        }

        // Sending the last row at k = kmax-2 to process 0
        tag = 6; // Tag for level k=kmax-2
        MPI_Send(my_phi[kmax-2][1], count, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
    }
}