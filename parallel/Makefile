# Makefile for Parallel MPI Gauss-Seidel

# MPI Compiler Wrapper
MPICC = mpicc
# Use MPICC as the C compiler
CC = $(MPICC)
CFLAGS = -Wall -Wextra -I./include 
LDFLAGS = -lm

# Directories 
SRC_PARALLEL_DIR = src
INCLUDE_PARALLEL_DIR = include
OBJ_DIR_PARALLEL = obj

# Source files
PARALLEL_SRCS = $(wildcard $(SRC_PARALLEL_DIR)/*.c)

# Object files
PARALLEL_OBJS = $(patsubst $(SRC_PARALLEL_DIR)/%.c,$(OBJ_DIR_PARALLEL)/%.o,$(PARALLEL_SRCS))

# Target
TARGET_PARALLEL = main_parallel

# Default target
all: parallel

# Parallel version
parallel: $(OBJ_DIR_PARALLEL) $(PARALLEL_OBJS)
	$(CC) $(CFLAGS) -fopenmp $(PARALLEL_OBJS) -o $(TARGET_PARALLEL) $(LDFLAGS)

# Create object directories
$(OBJ_DIR_PARALLEL):
	mkdir -p $(OBJ_DIR_PARALLEL)

# Compile parallel source files
$(OBJ_DIR_PARALLEL)/%.o: $(SRC_PARALLEL_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE_PARALLEL_DIR) -fopenmp -c $< -o $@

# Debug target
debug: CFLAGS += -g -DDEBUG
debug: clean all

clean:
	rm -rf $(OBJ_DIR_PARALLEL)
	rm -f $(TARGET_PARALLEL)

.PHONY: all parallel debug clean