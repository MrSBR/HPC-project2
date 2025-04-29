#include "../include/GSserial.h"

int main (int nargs, char **args)
{
    double ***arr1 = NULL, ***arr2 = NULL;
    int n, k, j, i, num_iters, kmax, jmax, imax;

    // --- Read command line arguments ---
    if (nargs != 5) {
        fprintf(stderr, "Usage: %s <num_iters> <kmax> <jmax> <imax>\n", args[0]);
        return EXIT_FAILURE;
    }

    num_iters = atoi(args[1]);
    kmax      = atoi(args[2]);
    jmax      = atoi(args[3]);
    imax      = atoi(args[4]);

    if (num_iters < 0 || kmax <= 2 || jmax <= 2 || imax <= 2) {
         fprintf(stderr, "Error: num_iters must be >= 0, dimensions (kmax, jmax, imax) must be > 2.\n");
         return EXIT_FAILURE;
    }
     // Ensure jmax is suitable for division (even number often simpler for jmax/2)
     if (jmax % 2 != 0) {
         fprintf(stderr, "Warning: jmax is odd (%d). Integer division jmax/2 = %d will be used.\n", jmax, jmax/2);
     }


    printf("Running Gauss-Seidel Comparison:\n");
    printf("  Iterations: %d\n", num_iters);
    printf("  Dimensions: kmax=%d, jmax=%d, imax=%d\n", kmax, jmax, imax);
    printf("  Interior points: %d x %d x %d\n", kmax-2, jmax-2, imax-2);


    // --- Allocate memory ---
    printf("Allocating arrays...\n");
    allocate_array3D(kmax, jmax, imax, &arr1);
    allocate_array3D(kmax, jmax, imax, &arr2);
    printf("Allocation complete.\n");

    // --- Initialize arrays ---
    printf("Initializing arrays with identical values...\n");

    // Seed random number generator ONCE
    srand(12345);

    for (k = 0; k < kmax; k++) {
        for (j = 0; j < jmax; j++) {
            for (i = 0; i < imax; i++) {
                // Assign some non-constant value, e.g., based on indices or random
                double val = (double)rand() / RAND_MAX * 100.0; // Example: random values 0-100
                // Or use index-based value: double val = k * 1.0 + j * 0.1 + i * 0.01;
                arr1[k][j][i] = val;
                arr2[k][j][i] = val; // Ensure both start identically
            }
        }
    }   
    printf("Initialization complete.\n");


    // --- Perform iterations ---
    printf("Performing %d iterations...\n", num_iters);
    for (n = 0; n < num_iters; n++) {
        GS_iteration_normal(kmax, jmax, imax, arr1);
        GS_iteration_2_chunks(kmax, jmax, imax, arr2);
         // Optional: Print progress
         // if ((n + 1) % 10 == 0) {
         //     printf("  Iteration %d completed.\n", n + 1);
         // }
    }
    printf("Iterations complete.\n");

    // --- Calculate and print difference ---
    printf("Comparing Serial result with normal result...\n");
    double diff = euclidean_distance(kmax, jmax, imax, arr1, arr2);

    printf("\n--- Final Comparison ---\n");
    printf("Num Iters: %d, kmax: %d, jmax: %d, imax: %d\n", num_iters, kmax, jmax, imax);
    printf("Euclidean Distance (Serial vs Normal): %g\n", diff);

    // Free global arrays (only allocated on Rank 0)
        for (int k=0; k<kmax; ++k) {
            for (int j=0; j<jmax; ++j) { free(arr1[k][j]); free(arr2[k][j]); }
            free(arr1[k]); free(arr2[k]);
        }
        free(arr1); free(arr2);


    return 0; 
}