# IN3200/IN4200 Oblig 2 - Gauss-Seidel Parallelization

This submission contains the code for Obligatory Assignment 2.

## Directory Structure

The submission is organized as follows:

-   `README.md`: This file.
-   `serial/`: Contains the source code for the serial implementation.
    -   `include/`: contains header file
        - `GSserial.h`
    -   `obj/`: contains object files
        ...
    -   `src/`: contains c files
        -   `allocatearray3D.c`
        -   `euclideandistance.c`: 
        -   `GSiteration2chuncks.c`: 
        -   `GSiterationnormal.c`:
        -   `main.c`:
    -   `main_serial`:
    -   `Makefile`:
-   `parallel/`: (To be added later) Contains the source code for the MPI parallel implementation.
    -   `include/`: contains header file
        - `GSMPI.h`
    -   `obj/`: contains object files
        ...
    -   `src/`: contains c files
        -   `allocatearray3D.c`
        -   `euclideandistance.c`: 
        -   `GSiteration2chunksmpi.c`: 
        -   `GSiteration2chunksserial.c`:
        -   `main_mpi.c`:
    -   `main_parallel`:
    -   `Makefile`:

## Implementation

This folder contains a serial C implementation of the Gauss-Seidel algorithm for a 3D grid. It includes:
1.  The standard iteration method (`GS_iteration_normal`).
2.  A wavefront-based iteration using two chunks (`GS_iteration_2_chunks`).
3.  Helper functions for 3D array allocation (`allocate_array3D`) and calculating the Euclidean distance (`euclidean_distance`) between two arrays.
4.  A main program (`main_serial.c`) that runs both iteration methods for a specified number of iterations and compares their results.

This folder also contains the MPI parallel implementation using 2 processes doing the same as the serial implementation. It contains the same as the serial implementation but with MPI modifications.

### Compilation

To compile the serial code:
1.  Navigate to the `serial` directory:
    ```bash
    cd IN4200_Oblig2_xxx/serial
    ```
2.  Run the `make` command:
    ```bash
    make
    ```
    This will create an executable file named `main_serial`.


The MPI implementation requires installation of OpenMPI or MPICH such that the `mpicc` compiler wrapper can be made available in your PATH.
1.  Navigate to the `parallel` directory:
    ```bash
    cd IN4200_Oblig2_xxx/parallel
    ```
2.  Run the `make` command from the parallel folder:
    ```bash
    make
    ```
    This will use `mpicc` to compile the code and create an executable file named `main_parallel` inside the `parallel` directory (`parallel/main_parallel`).

### Running

To run the serial and parallel executable respectivley:
```bash
./main_serial <num_iters> <kmax> <jmax> <imax>

forexample:
```bash
./main_serial 100 3 4 3  

```bash
# General format for 2 processes
mpirun -np 2 main_parallel <num_iters> <kmax> <jmax> <imax>

# Example:
mpirun -np 2 main_parallel 100 3 4 3