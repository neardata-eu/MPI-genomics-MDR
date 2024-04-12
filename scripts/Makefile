# Makefile for compiling mdr-mpi.c with mpicc

# Compiler
CC = mpicc

# Compiler flags
CFLAGS = -Wall -Wextra -O3 

# Libraries
LDLIBS = -lz

# Source files
SRCS = mdr-mpi.c

# Object files
OBJS = $(SRCS:.c=.o)

# Executable name
EXEC = mdr-mpi-c

# Main target
all: $(EXEC)

# Compile rule for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link rule for the executable
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LDLIBS)

# Clean rule
clean:
	rm -f $(EXEC) $(OBJS)

