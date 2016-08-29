SHELL = /bin/sh
#F90 = ifort
#OPT = -O3 -xhost -implicitnone -ftrapuv -traceback -check bounds -warn all -openmp-stubs
#OPT = -implicitnone -O3 -xhost -mkl -openmp-stubs -pad
#OPT = -implicitnone -O3 -xhost -pad
#OPT = -O0 -implicitnone -ftrapuv -traceback -check bounds -openmp -warn all
#OPT = -implicitnone -O3 -xhost -openmp -pad -mkl
#OPT = -implicitnone -O3 -xhost -openmp -pad -ipo -mkl

# Code currently doesn't compile under gfortran because of our use of [] for initializing
# two-dimensional array constants
F90 = gfortran
#OPT = -ffree-form -ffree-line-length-none -std=f2008ts -O3 -fopenmp -ffpe-trap=invalid -fbacktrace
OPT = -ffree-form -ffree-line-length-none -std=f2008ts -O3 -march=native -fopenmp

all:ptm

rksuite.o: rksuite.f90
	$(F90) $(OPT) -c -o rksuite.o rksuite.f90

global.o: global.f90
	$(F90) $(OPT) -c -o global.o global.f90

finite_differences.o: finite_differences.f90
	$(F90) $(OPT) -c -o finite_differences.o finite_differences.f90

interpolation.o: interpolation.f90
	$(F90) $(OPT) -c -o interpolation.o interpolation.f90

fields.o: fields.f90
	$(F90) $(OPT) -c -o fields.o fields.f90

particles.o: particles.f90
	$(F90) $(OPT) -c -o particles.o particles.f90

fileio.o: fileio.f90
	$(F90) $(OPT) -c -o fileio.o fileio.f90

stepper.o: stepper.f90
	$(F90) $(OPT) -c -o stepper.o stepper.f90

ptm.o: ptm.f90
	$(F90) $(OPT) -c -o ptm.o ptm.f90

ptm: rksuite.o global.o finite_differences.o interpolation.o fields.o particles.o fileio.o stepper.o ptm.o
	$(F90) $(OPT) -o ptm rksuite.o global.o finite_differences.o interpolation.o fields.o particles.o fileio.o stepper.o ptm.o

clean:
	rm -f *.o
	rm -f *.mod

