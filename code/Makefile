########################################
# Makefile for mpcSim		       #
#			               #
# Leonor Garcia-Gurierrez	       #
# May 2013			       #
########################################
# C compiler and flags; compiler for MPI: mpicc (wrapper, recall NAG course). U
# GNU compiler: gcc. Load modules for Intel (icc) and PGI
CC = gcc
CFLAGS = -O3 -Wall -Wextra -Wundef -Wshadow -Wstrict-overflow=5 -Winline #-fno-inline -pg #-Wfloat-equal
INC = -I.

# Command to use for linking and executable
# LD=mpicc   # use this for MPI
# LD is the linker. EXE is the name of the executable that will be created when linking the files 
LD = gcc 
LDFLAGS = -lm -lgsl -lgslcblas #-pg # Include the maths library. For MPI, possibly: -lpng -nofor-main
# Executable name
EXE = mpcSim
# List of objects
OBJECTS = mpc.o stream.o collide.o

# Input data file
INPUT = input.dat

# Ignore any file named clean and check_mem, if any
.PHONY: clean
.PHONY: memory
.PHONY: doc
.PHONY: profile
.PHONY: rmbackup

# Rules for building object files
# Each object if built from a source file of the same name (dependency), using the compiler
# specified in the macro CC, with compiler flags CFLAGS, is assigned the name of this target (-o $@)
# The source files are compiled or assembled but not linked (-c)
# and $< stands for the .c file used to build the target .o
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $(INC) $<
	
# Build the executable by linking the object files.
$(EXE) : $(OBJECTS) 
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) 
	#@echo "LDFLAGS="$(LDFLAGS) 
	#@echo "CFLAGS="$(CFLAGS)

#------------------------------------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------------------------------------
# Purge object files and executable
clean :
	rm *.o ./$(EXE)

# Use Valgrind to track down double free or corruption or the like
memory :
	make CFLAGS='-O1 -g'
	valgrind --leak-check=yes ./$(EXE) $(INPUT)

# Generate documentation with doxygen
doc :
	doxygen Doxyfile
	
# Profile the code: force no inlining and remove optimization flags
profile:
	make CFLAGS='-O0 -Wall -Wextra -Wundef -Wshadow -Wstrict-overflow=5 -fno-inline -pg' LDFLAGS='-pg -lm -lgsl -lgslcblas'

# Remove system and emacs backup files
rmbackup:
	rm *~
	rm \#*
