# IN3200/IN4200 Oblig 2 - Gauss-Seidel Parallelization

This submission contains the code for Obligatory Assignment 2.

## Directory Structure

The submission is organized as follows:

-   `README.md`: This file.
-   `serial/`: Contains the source code for the serial implementation.
    -   `*.c`: C source files (one per function requested in task 3.1).
    -   `gs_serial.h`: Header file for the serial code.
    -   `Makefile`: To compile the serial code.
-   `parallel/`: (To be added later) Contains the source code for the MPI parallel implementation.
    -   `*.c`: C source files for the parallel version.
    -   `gs_parallel.h`: Header file for the parallel code.
    -   `Makefile`: To compile the parallel code.

## Serial Implementation (`serial/`)

This folder contains a serial C implementation of the Gauss-Seidel algorithm for a 3D grid, strictly following the functions requested in task section 3.1. It includes:
1.  The standard iteration method (`GS_iteration_normal`).
2.  A wavefront-based iteration using two chunks (`GS_iteration_2_chunks`).
3.  Helper functions for 3D array allocation (`allocate_array3D`) and calculating the Euclidean distance (`euclidean_distance`) between two arrays.
4.  A main program (`main_serial.c`) that runs both iteration methods for a specified number of iterations and compares their results.

**Important Note:** Following the strict list of required functions in the task description, memory allocated for the 3D arrays is *not* explicitly freed upon program termination. This will result in memory leaks, which is generally discouraged in production code but adheres to the specific constraints of this assignment section.

### Compilation

To compile the serial code:
1.  Navigate to the `serial` directory:
    ```bash
    cd IN3200_Oblig2_xxx/serial
    ```
    or
    ```bash
    cd IN4200_Oblig2_xxx/serial
    ```
2.  Run the `make` command:
    ```bash
    make
    ```
    This will create an executable file named `gs_serial`.

### Running

To run the serial executable:
```bash
./main_serial <num_iters> <kmax> <jmax> <imax>

forexample:
```bash
./main_serial 100 3 4 3  

## Parallel Implementation (`parallel/`)

This folder contains the MPI parallel implementation using 2 processes and domain decomposition in the j-direction.

### Compilation

Requires an MPI implementation (like OpenMPI, MPICH) to be installed and the `mpicc` compiler wrapper available in your PATH.
1.  Navigate to the `parallel` directory:
    ```bash
    cd IN3200_Oblig2_xxx/parallel
    ```
    or
    ```bash
    cd IN4200_Oblig2_xxx/parallel
    ```
2.  Run the `make` command from the parallelGS folder:
    ```bash
    make
    ```
    This will use `mpicc` to compile the code and create an executable file named `gs_parallel` inside the `exe` sub-directory (`parallel/exe/gs_parallel`).

    Run the 'clean' command from the exe folder:
    ``bash
    make clean
    ```
    This will remove the object and dependency files in the exe folder
### Running

Use the `mpirun` or `mpiexec` command (provided by your MPI installation) to launch the program with exactly 2 processes.
```bash
# General format
mpirun -np 2 main_parallel <num_iters> <kmax> <jmax> <imax>

# Example:
mpirun -np 2 main_parallel 100 3 4 3