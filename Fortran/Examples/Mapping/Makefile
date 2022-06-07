FC      := gfortran 
CC      := gcc
LD      := gfortran
F_OPT   :=-O -Wall -w  -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero,overflow -freal-4-real-8 -Wno-argument-mismatch -fno-range-check
F_PREC  :=-fdefault-real-8

FFLAGS  := $(F_OPT)
CPPFLAGS := -fpp
INCLUDE :=-I./

TEST_INTERP:= test_interp
OBJS_BASIC := machine.o 
INTERP := testInterpolation.o nms.o mod_legendre.o mod_adaptiveInterpolation.o $(OBJS_BASIC)

all : $(TEST_INTERP)


testInterpolation.o        : nms.o machine.o mod_legendre.o mod_adaptiveInterpolation.o
mod_legendre.o : $(OBJS_BASIC)


clean :
	rm -f $(OBJS) *.mod *.out *.o  *.dat


$(TEST_INTERP) : $(INTERP)
	$(LD) $(FFLAGS) $(INTERP) -o a.out
	#./a.out &> out

#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE)  $< -o $@
