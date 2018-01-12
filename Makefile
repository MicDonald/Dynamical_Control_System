.DEFAULT_GOAL := LIB

#C_ROOT = /opt/gcc-5.2.0
C_ROOT = /opt
#MPI_ROOT = /opt/openmpi/gcc-7.2.0
MPI_ROOT = /opt
MPICC = $(MPI_ROOT)/bin/mpicxx
CC = $(C_ROOT)/bin/g++
CFLAGS = -std=c++14 -O3 -Wall -Wextra -Wno-sign-compare -DMKL_ILP64 -m64 -I/opt/include/ -I/opt/intel/mkl/include/ -I/nas/home/mkdonald/lammps/lammps-14May16/src/USER-DYNAMICAL_CTRLSYSTEM
MPI_CFLAGS = --std=c++14 -O3 -Wextra -DMKL_ILP64 -DPARALLEL -m64 -I/opt/include -I/opt/intel/mkl/include/ -I/nas/home/mkdonald/lammps/lammps-14May16/src/USER-DYNAMICAL_CTRLSYSTEM
MPI_LDFLAGS = -lpotentialPack -L/opt/intel/mkl/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_openmpi_ilp64 -lpthread -lm -Wl,--rpath=/opt/openmpi/gcc-7.2.0/lib,--rpath=/opt/lib/gcc/x86_64-pc-linux-gnu7.2.0,--rpath=$(C_ROOT)/lib64

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
	rm -f $(LIBDIR)/libDynamicalCtrlSystem.a
	rm -rf $(OBJDIR)/*.o
	rm -f $(OBJDIR)/console

LIB: $(OBJS)
	ar rvs $(LIBDIR)/libDynamicalCtrlSystem.a $(OBJS)

mpi: $(OBJS) $(MPI_OBJS) 
	ar rvs $(LIBDIR)/libDynamicalCtrlSystem.a $(OBJS) $(MPI_OBJS) 

omp: $(OBJS) $(OMP_OBJS)
	ar rvs $(LIBDIR)/libDynamicalCtrlSystem.a $(OBJS) $(OMP_OBJS)

test: $(OBJDIR)/main.o $(OBJS)
	$(CC) $(OBJDIR)/main.o $(OBJS) -Wl,--rpath=$(C_ROOT)/lib64,--rpath=$(MPI_ROOT)/lib -o $(OBJDIR)/console

