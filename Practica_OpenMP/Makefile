#
# Exact genetic sequence alignment
#
# Parallel computing (Degree in Computer Engineering)
# 2023/2024
#
# (c) 2024 Arturo Gonzalez-Escribano
# Grupo Trasgo, Universidad de Valladolid (Spain)
#

# Compilers
CC=gcc
OMPFLAG=-fopenmp

# Flags for optimization and external libs
FLAGS=-O3
LIBS=-lm

# Targets to build
OBJS=R2D2

# Rules. By default show help
help:
	@echo
	@echo "Exact genetic sequence alignment"
	@echo
	@echo "Group Trasgo, Universidad de Valladolid (Spain)"
	@echo
	@echo "make align_seq	Build only the sequential version"
	@echo "make align_omp	Build only the OpenMP version"
	@echo
	@echo "make all	Build all versions (Sequential, OpenMP)"
	@echo "make debug	Build all version with demo output for small surfaces"
	@echo "make clean	Remove targets"
	@echo

all: $(OBJS)

R2D2: R2D2.c
	$(CC) $(FLAGS) $(DEBUG) $(OMPFLAG) $< $(LIBS) -o $@

# Remove the target files
clean:
	rm -rf $(OBJS)

# Compile in debug mode
debug:
	make DEBUG="-DDEBUG -g" FLAGS= all
