.DEFAULT_GOAL := LIB

C_ROOT = /opt/gcc-5.2.0
MPI_ROOT = /opt/openmpi/1.8.7/gnu520
MPICC = $(MPI_ROOT)/bin/mpicxx
CC = $(C_ROOT)/bin/g++
CFLAGS = -std=c++14 -O3 -Wall -Wextra -DMKL_ILP64 -m64 -I/opt/intel/Compiler/11.0/074/mkl/include
MPI_CFLAGS = --std=c++14 -O3 -Wextra -DMKL_ILP64 -DPARALLEL -m64 -I/opt/intel/Compiler/11.0/074/mkl/include
MPI_LDFLAGS = -lpotentialPack -L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_openmpi_ilp64 -lpthread -lm -Wl,--rpath=/opt/openmpi/1.10.1/gcc520/lib,--rpath=/opt/intel/mkl/lib/intel64,--rpath=$(C_ROOT)/lib64

OBJDIR = .
LIBDIR = .

CPPFILES = KernelMatrix.cpp \

#MPI_CPPFILES = KernelMatrix.cpp
#OMP_CPPFILES = KernelMatrix.cpp

OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPPFILES))
MPI_OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(MPI_CPPFILES))
OMP_OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(OMP_CPPFILES))


$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c -fopenmp $< -o $@

#$(OBJDIR)/%.o: %.cpp
#	$(MPICC) $(MPI_CFLAGS) $(MPI_LDFLAGS) -c -I$(MPI_ROOT) $< -o $@

#$(OBJDIR)/%.o: %OMP.cpp
#	$(CC) $(OMP_CFLAGS) -c -fopenmp $< -o $@

clean:
	rm -f $(LIBDIR)/libVRTransition.a
	rm -rf $(OBJDIR)/*.o
	rm -f $(OBJDIR)/console

LIB: $(OBJS)
	ar rvs $(LIBDIR)/libVRTransition.a $(OBJS)

mpi: $(OBJS) $(MPI_OBJS) 
	ar rvs $(LIBDIR)/libVRTransition.a $(OBJS) $(MPI_OBJS) 

omp: $(OBJS) $(OMP_OBJS)
	ar rvs $(LIBDIR)/libVRTransition.a $(OBJS) $(OMP_OBJS)

test: $(OBJDIR)/main.o $(OBJS)
	$(CC) $(OBJDIR)/main.o $(OBJS) -Wl,--rpath=$(C_ROOT)/lib64,--rpath=$(MPI_ROOT)/lib -o $(OBJDIR)/console

