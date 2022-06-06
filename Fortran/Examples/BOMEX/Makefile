#FC     := ifort
FC     := gfortran
#CC     := icc
CC     := gcc
#LD     := ifort
LD     := gfortran

#F_OPT   :=-O2 -g -traceback -fp-model strict -ftz -fpe0 -align array64byte
F_OPT   :=-O -Wall -w  -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero,overflow -freal-4-real-8 -Wno-argument-mismatch -fno-range-check#-ftz -fpe0 -align array64byte
#F_PREC  :=-r8
F_PREC  :=-fdefault-real-8

FFLAGS  :=$(F_PREC) $(F_OPT)
INCLUDE :=-I./

MAIN    := main
PROGRAM := a.out

OBJS_BASIC := machine.o physcons.o mod_advection.o
OBJS := main.o mod_bomex.o mod_adaptiveInterpolation.o nms.o mfpbl.o moninedmf.o scm_physics.o kessler.o mod_legendre.o $(OBJS_BASIC)

all : $(MAIN)


main.o         : mod_bomex.o moninedmf.o scm_physics.o mod_adaptiveInterpolation.o mod_legendre.o $(OBJS_BASIC)
mod_bomex.o    : $(OBJS_BASIC)
mod_adaptiveInterpolation.o  : mod_legendre.o $(OBJS_BASIC)
moninedmf.o    : mfpbl.o $(OBJS_BASIC)
mfpbl.o        : $(OBJS_BASIC)
mod_legendre.o : $(OBJS_BASIC)
scm_physics.o  : kessler.o moninedmf.o $(OBJS_BASIC)


clean :
	rm -f $(OBJS) *.mod a.out

$(MAIN) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE) $(CPPFLAGS) $< -o $@
