#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/GSMPI.h"

void handle_input(int nargs, char **args, int *num_iters, int *kmax, int *jmax, int *imax, bool *print_verbose) {
    // Default values
    *num_iters = 1000;
    *kmax = 4;
    *jmax = 6;
    *imax = 8;

    // Check number of arguments 
    if (nargs == 5 || nargs == 6) {
        // Read from command line the values of num_iters, kmax, jmax, and imax
        *num_iters = atoi(args[1]);
        *kmax = atoi(args[2]);
        *jmax = atoi(args[3]);
        *imax = atoi(args[4]);
        *print_verbose = false;
        
        // Check if the values are valid
        if (*num_iters <= 0 || *kmax <= 2 || *jmax <= 3 || *imax <= 2) {
            fprintf(stderr, "***********************************************************\n");
            fprintf(stderr, "* Invalid input values. Valid values are:                 *\n");
            fprintf(stderr, "* num_iters > 0                                           *\n");
            fprintf(stderr, "* kmax > 2                                                *\n");
            fprintf(stderr, "* jmax > 3 (and divisible by 2)                           *\n");
            fprintf(stderr, "* imax > 2                                                *\n");
            fprintf(stderr, "* Usage: %s [num_iters] [kmax] [jmax] [imax] *\n", args[0]);
            fprintf(stderr, "* Received: num_iters=%d, kmax=%d, jmax=%d, imax=%d         *\n", *num_iters, *kmax, *jmax, *imax);
            fprintf(stderr, "***********************************************************\n\n");
            exit(EXIT_FAILURE);
        }

        if (nargs == 6) {
            // Check if the verbose flag is set
            if (args[5][0] == '-' && args[5][1] == 'v') {
                *print_verbose = true;
            } else {
                fprintf(stderr, "Invalid last argument flag. Use -v for verbose output.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    // If no arguments or just the verbose flag, is provided, use defualt values
    else if (nargs == 1 || nargs == 2) {
        char *message;
        int len_message = asprintf(&message, "* No args provided. Using default values: num_iters=%d, kmax=%d, jmax=%d, imax=%d *", *num_iters, *kmax, *jmax, *imax);

        for (int i = 0; i < len_message; i++) {printf("*");}
        printf("\n"); 
        printf("*");
        for (int i = 0; i < len_message-2; i++) {printf(" ");}
        printf("*\n");
        printf("%s\n", message);
        printf("*");
        for (int i = 0; i < len_message-2; i++) {printf(" ");}
        printf("*\n");
        for (int i = 0; i < len_message; i++) {printf("*");}
        printf("\n\n");


        if (nargs == 2) {
            // Check if the verbose flag is set
            if (args[1][0] == '-' && args[1][1] == 'v') {
                *print_verbose = true;
            } else {
                fprintf(stderr, "Invalid last argument flag. Use -v for verbose output.\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    // If the number of arguments is not 1, 2, 5, or 6, print an error message
    else {
        char *message;
        int len_message = asprintf(&message, "* Invalid number of arguments. Enter no args for default values (num_iters=%d, kmax=%d, jmax=%d, imax=%d) or 4 args for custom values. *", *num_iters, *kmax, *jmax, *imax);
        for (int i = 0; i < len_message; i++) {fprintf(stderr, "*");}
        fprintf(stderr, "\n");
        fprintf(stderr, "*");
        for (int i = 0; i < len_message-2; i++) {fprintf(stderr, " ");}
        fprintf(stderr, "*\n");
        fprintf(stderr, "%s\n", message);
        fprintf(stderr, "*");
        for (int i = 0; i < len_message-2; i++) {fprintf(stderr, " ");}
        fprintf(stderr, "*\n");
        for (int i = 0; i < len_message; i++) {fprintf(stderr, "*");}
        fprintf(stderr, "\n\n");
        exit(EXIT_FAILURE);
    }
}