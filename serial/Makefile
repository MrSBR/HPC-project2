# Makefile for Serial Gauss-Seidel

# Compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -I./include 
LDFLAGS = -lm

# Directories 
SRC_SERIAL_DIR = src
INCLUDE_SERIAL_DIR = include
OBJ_DIR_SERIAL = obj

# Source files
SERIAL_SRCS = $(wildcard $(SRC_SERIAL_DIR)/*.c)

# Object files
SERIAL_OBJS = $(patsubst $(SRC_SERIAL_DIR)/%.c,$(OBJ_DIR_SERIAL)/%.o,$(SERIAL_SRCS))

# Target
TARGET_SERIAL = main_serial

# Default target
all: serial 

# Serial version
serial: $(OBJ_DIR_SERIAL) $(SERIAL_OBJS)
	$(CC) $(CFLAGS) $(SERIAL_OBJS) -o $(TARGET_SERIAL) $(LDFLAGS)

# Create object directories
$(OBJ_DIR_SERIAL):
	mkdir -p $(OBJ_DIR_SERIAL)

# Compile source files
$(OBJ_DIR_SERIAL)/%.o: $(SRC_SERIAL_DIR)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE_SERIAL_DIR) -c $< -o $@

# Debug target
debug: CFLAGS += -g -DDEBUG
debug: clean all

clean:
	rm -rf $(OBJ_DIR_SERIAL)
	rm -f $(TARGET_SERIAL)

.PHONY: all serial debug clean