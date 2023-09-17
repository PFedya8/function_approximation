# makefile for approximation.f90

# Compiler
FC = gfortran-11

# Flags
FFLAGS = -Wall -Wextra -pedantic -fcheck=all -fbacktrace -fimplicit-none -Og -g -std=f2018 -Wconversion-extra -fmax-errors=10

# Objects
OBJS = approximation.o

# Executable
EXEC = approximation

# File 
FILE = approximation.f90

# Default rule (first one is default)
all: $(EXEC)

# Compile the object file
$(OBJS): $(FILE)
	$(FC) $(FFLAGS) -c $(FILE)

# Build the executable
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(EXEC)

# Clean
clean:
	rm -f $(OBJS) $(EXEC)

# Run
run: $(EXEC)
	./$(EXEC)

# Mark rules as phony
.PHONY: all clean run
