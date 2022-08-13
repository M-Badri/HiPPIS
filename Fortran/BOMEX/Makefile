## intel compiler configurations
#FC      := ifort
#F_OPT   :=-O3 -xMIC-AVX512 -fp-model strict  -ftz -fpe0 -traceback -qopenmp-simd -align array64byte  
#F_PREC  :=-r8

## gfortran compiler configurations
FC     := gfortran
F_OPT   :=-O2 -Wall -fbacktrace -w -freal-4-real-8 -fno-range-check 
F_PREC  :=-fdefault-real-8

FFLAGS  :=$(F_PREC) $(F_OPT)
INCLUDE :=-I./

MAIN    := main

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
	rm -f $(OBJS) *.mod main

$(MAIN) : $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) -o $(MAIN)
	mkdir -p bomex_data
	./$(MAIN)
        

#------------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------------
%.o : %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $< -o $@

%.o : %.F90
	$(FC) -c $(FFLAGS) $(INCLUDE) $(CPPFLAGS) $< -o $@
