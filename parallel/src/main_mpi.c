#include "../include/GSMPI.h"

int main(int nargs, char **args) {
    int num_iters, kmax, jmax, imax;

    // Initialize MPI and make sure there are 2 processes
    MPI_Init(&nargs, &args);
    
    int my_rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (num_procs != 2) {
        fprintf(stderr, "This program requires exactly 2 processes. Read the README.md to see how this is done\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int num_values;                         // How many values to send/recieve with MPI 
    MPI_Datatype data_type;                 // Data type to send/recieve 
    int destination;                        // Rank of the recieving process
    int source;                             // Rank of the sending process
    int tag;                                // Tag for the message, which works as an ID to link specific send and recieve. 
    MPI_Comm comn = MPI_COMM_WORLD;         // Communicator for the message; world communicator
    MPI_Status *status = MPI_STATUS_IGNORE; // Pointer to status object containing info on the recieved message. We ignore status since it's not used
    
    // Dividing the main file into rank 0 and rank 1, where rank 0 is the root node
    if (my_rank == 0) {
        //Initializing input from terminal
        num_iters = atoi(args[1]);
        kmax = atoi(args[2]);
        jmax = atoi(args[3]);
        imax = atoi(args[4]);

        printf("Running Gauss-Seidel Comparison:\n");
        printf("  Iterations: %d\n", num_iters);
        printf("  Dimensions: kmax=%d, jmax=%d, imax=%d\n", kmax, jmax, imax);
        printf("  Interior points: %d x %d x %d\n", kmax-2, jmax-2, imax-2);

        source = 1; 
        destination = 1; 

        int num_values = 1;
        MPI_Datatype data_type = MPI_INT;
        int tag = 1; 
        // Send the input values to process 1
        printf("Sending input values to process 1\n");
        MPI_Send(&num_iters, num_values, data_type, destination, tag, comn);
        MPI_Send(&kmax, num_values, data_type, destination, tag, comn);
        MPI_Send(&jmax, num_values, data_type, destination, tag, comn);
        MPI_Send(&imax, num_values, data_type, destination, tag, comn);

        // Allocate memory for the first 3D array on process 0
        double ***array_mpi; 
        double ***array_serial;
        int jmid = jmax/2 + 1;

        printf("Allocating and initializing first half of the 3D array \n");
        allocate_array3D(kmax, jmid, imax, &array_mpi);
        
        // Initialize the first half of the 3D array 
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                for (int i = 0; i < imax; i++) {
                    array_mpi[k][j][i] = pow(k*jmax*imax + 
                                            j*imax + 
                                            i, 2); 
                }
            }
        }
        
        allocate_array3D(kmax, jmax, imax, &array_serial);
        printf("Allocating and initializing first half of the 3D array complete\n");


        // Compute the first half of the 3D array
        printf("Performing %d iterations on first half of the 3D array\n", num_iters);
        for (int _ = 0; _ < num_iters; _++) {
            GS_iteration_2_chunks_mpi(my_rank, kmax, jmid, imax, array_mpi);
            GS_iteration_2_chunks(kmax, jmax, imax, array_serial);
        }
        printf("Iterations complete.\n");

        // Receive the results from process 1 and store in temp array 
        printf("Receiving results from process 1.\n");
        double ***recieved_array;
        allocate_array3D(kmax, jmid, imax, &recieved_array);

        num_values = imax;
        data_type = MPI_DOUBLE;
        tag = 2; 
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                MPI_Recv(recieved_array[k][j], num_values, data_type, source, tag, comn, status);
            }
        }

        // Construct the global array from the two halves
        printf("Creating and filling the global array\n");
        double ***global_array;
        allocate_array3D(kmax, jmax, imax, &global_array);
        for (int k = 1; k < kmax-1; k++) {
            for (int j = 1; j < jmax-1; j++) {
                for (int i = 1; i < imax-1; i++) {
                if (j < jmid-2) {
                        // Copy the left half from process 0
                        global_array[k][j][i] = array_mpi[k][j][i];
                    } else {
                        // Copy the right half from process 1
                        global_array[k][j][i] = recieved_array[k][j-jmid+2][i];
                    }
                }
            }
        }
        printf("Creating and filling the global array complete\n");
        printf("\n--- Final Comparison ---\n");
        printf("Num Iters=%d, kmax=%d, jmax=%d, imax=%d \n", num_iters, kmax, jmax, imax);
        printf("EuclideanDistance (Serial vs Normal): %g \n", euclidean_distance(kmax, jmax, imax, array_serial, global_array));

        // Free the allocated arrays
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                free(array_mpi[k][j]);
                free(array_serial[k][j]);
                free(global_array[k][j]);
                free(recieved_array[k][j]);
            }
            free(array_mpi[k]);
            free(array_serial[k]);
            free(global_array[k]);
            free(recieved_array[k]);
        }
        free(array_mpi);
        free(array_serial);
        free(global_array);
        free(recieved_array);
    }
 
    if (my_rank == 1) {
        // Receive the input values from process 0
        source = 0;
        destination = 0;
        
        num_values = 1;
        data_type = MPI_INT;
        tag = 1;
        MPI_Recv(&num_iters, num_values, data_type, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&kmax, num_values, data_type, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&jmax, num_values, data_type, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&imax, num_values, data_type, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Allocate and initialize the second 3D array on process 1
        double ***array_mpi;
        int jmid = jmax/2 + 1;
        allocate_array3D(kmax, jmid, imax, &array_mpi);
        
        // Initialize the second half of the 3D array 
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                for (int i = 0; i < imax; i++) {
                    array_mpi[k][j][i] = pow(k*jmax*imax + 
                                            j*imax + 
                                            i + (jmid-2)*imax, 2); 
                }
            }
        }
        
        // Compute the second half of the 3D array
        for (int _ = 0; _ < num_iters; _++) {
            GS_iteration_2_chunks_mpi(my_rank, kmax, jmid, imax, array_mpi);
        }

        // Send the results back to process 0
        num_values = imax;
        data_type = MPI_DOUBLE;
        tag = 2;
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                MPI_Send(array_mpi[k][j], num_values, data_type, destination, tag, comn);
            }
        }

        // Free the allocated arrays
        for (int k = 0; k < kmax; k++) {
            for (int j = 0; j < jmid; j++) {
                free(array_mpi[k][j]);
            }
            free(array_mpi[k]);
        }
        free(array_mpi);
    }
    MPI_Finalize();
    return 0;
}