FC      := gfortran 
CC      := gcc
LD      := gfortran
F_OPT   :=-Ofast -Wall -w -fopenmp  -ftree-vectorize -ffast-math -mavx2 -mavx -fbacktrace -fno-trapping-math -freal-4-real-8 -Wno-argument-mismatch -fno-range-check
#F_OPT   :=-Ofast -Wall -w -fopenmp  -ftree-vectorize -ffast-math -mavx2 -mavx -fopt-info-vec-all -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero,overflow -freal-4-real-8 -Wno-argument-mismatch -fno-range-check
F_PREC  :=-fdefault-real-8

FFLAGS  := $(F_OPT)
CPPFLAGS := -fpp
INCLUDE :=-I./

TEST_MAIN:= test_main
OBJS_BASIC := machine.o 
MAIN := main.o nms.o mod_legendre.o mod_adaptiveInterpolation.o $(OBJS_BASIC)
TUTORIAL := tutorial.o mod_adaptiveInterpolation.o

all : $(TEST_MAIN)


main.o        : nms.o machine.o mod_legendre.o mod_adaptiveInterpolation.o
tutorial.o    : mod_adaptiveInterpolation.o
mod_legendre.o : $(OBJS_BASIC)


clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat


$(TEST_MAIN) : $(MAIN) $(TUTORIAL)
	$(LD) $(FFLAGS) $(MAIN) -o a.out
	$(LD) $(FFLAGS) $(TUTORIAL) -o tutorial
	#./a.out &> out

#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@
