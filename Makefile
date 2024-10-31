# Compiler to be used
CC = gcc

# Compiler flags
CFLAGS = -g -fopenmp

# Libraries to link against
LIBS = -lm -lgsl -lgslcblas

# The target executable
TARGET = rec_plots_chaotic_systems

# Source files
SRCS = rec_plots_chaotic_systems.c \
       auxiliary_and_recurrence_functions.c

# Object files (derived from source files)
OBJS = $(SRCS:.c=.o)

# File to signify build time
TIMESTAMP_FILE = .buildtimestamp

# Default target to build and run
run: $(TARGET)
	@./$(TARGET)

# Rule to link object files to create the executable
$(TARGET): $(OBJS) $(TIMESTAMP_FILE)
	@$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Rule to compile .c files into .o files
%.o: %.c
	@$(CC) $(CFLAGS) -c $< -o $@

# Create timestamp file to signify build time
$(TIMESTAMP_FILE):
	@touch $(TIMESTAMP_FILE)

# Clean up the build (remove object files, the executable, and the timestamp file)
clean:
	@rm -f $(OBJS) $(TARGET) $(TIMESTAMP_FILE)

# Phony targets to avoid conflicts with files of the same name
.PHONY: run clean
